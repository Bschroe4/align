PROGRAM align

  implicit none

  double precision, parameter :: factr=0.529177240859d0
  double precision, parameter :: factm=1822.888484
  character(*),     parameter :: cvers="XVERSION"


  logical :: lexist

  integer :: nat,iat
  integer :: ixyz
  integer :: irep
  integer :: ncmd,icmd,ifil
  integer, allocatable, dimension(:)   :: mask

  character(len=100) :: cbuf
  character(len=32)  :: ccmd,clist
  character(len=3)   :: catlab,catom
  character(len=3), allocatable, dimension(:)   :: ztag

  double precision :: xmass
  double precision, allocatable, dimension(:)   :: atmass
  double precision, allocatable, dimension(:,:) :: xyzref
  double precision, allocatable, dimension(:,:) :: xyz
  double precision, allocatable, dimension(:,:) :: tmp
  double precision, dimension(3,3) :: u
  double precision, dimension(3,3) :: I,Xref
  double precision, dimension(3)   :: ie,be 

  integer, dimension(3,6) :: ireptab_abc2xyz = RESHAPE((/ 3, 1, 2,  & !   IR
                                                          2, 3, 1,  & !  IIR
                                                          1, 2, 3,  & ! IIIR
                                                          3, 2, 1,  & !   IL
                                                          1, 3, 2,  & !  IIL
                                                          2, 1, 3/),& ! IIIL
                                                          shape(ireptab_abc2xyz))
!
  !1) read first argument i.e. xyzref unto which shall be aligned
  clist=""
  ncmd=COMMAND_ARGUMENT_COUNT()
  ! add KILL/CATCH if ncmd zero
  icmd=0
  ifil=0
  do 
    icmd=icmd+1
    if(icmd.gt.ncmd) exit
    call get_command_argument(icmd,ccmd)
    select case (trim(ccmd))
      case ("-h","--help")
        ! print help
      case ("-f","--fragment") ! align to a specific fragment 
        icmd=icmd+1
        call get_command_argument(icmd,clist)
      case default ! file designations
        ifil=ifil+1
        select case (ifil)
          case (1)
!           call get_command_argument(1,cfile)
            inquire(file=ccmd,exist=lexist)
            if(.not.lexist) call file_error(1,ccmd)
            open(unit=11,file=ccmd,status='OLD')
            read(11,'(I4)') nat
            allocate(xyzref(3,nat))
            allocate(atmass(nat))
            allocate(ztag(nat))
            read(11,*) cbuf
            DO IAT=1,NAT
              read(11,*) catom,(xyzref(ixyz,iat),ixyz=1,3)
              catlab=StrUpCase(catom,len(catom),1)
              atmass(iat)=convert_atlab_to_atmass(catlab)*factm
              ztag(iat)=StrUpCase(catom,len(catom),0)
              call dscal(3,1.d0/factr,xyzref(1,iat),1)
            ENDDO
            close(unit=11)
          case(2)
            !2) read second argument i.e. xyz which shall be aligned
            inquire(file=ccmd,exist=lexist)
            if(.not.lexist) call file_error(1,ccmd)
            open(unit=11,file=ccmd,status='OLD')
            read(11,*) iat
            if(iat.ne.nat) call nat_error(nat,iat)
            allocate(xyz(3,nat))
            read(11,*) cbuf
            DO IAT=1,NAT
              read(11,*) catom,(xyz(ixyz,iat),ixyz=1,3)
              catlab=StrUpCase(catom,len(catom),1)
              xmass=convert_atlab_to_atmass(catlab)*factm
              if(xmass.ne.atmass(iat)) call mass_error(iat,xmass/factm,atmass(iat)/factm)
              call dscal(3,1.d0/factr,xyz(1,iat),1)
            ENDDO
            close(unit=11)
          case(3)
            open(unit=11,file=ccmd,status='REPLACE')
        end select
    end select
  end do
!
  allocate(mask(nat))
  mask=1
  call convert_list2mask(nat,clist,mask)
  call transform_com(nat,mask,atmass,xyzref)
  allocate(tmp(3,nat))
  tmp=xyz
  call transform_com(nat,mask,atmass,tmp)
  call quaternion(nat,tmp,xyzref,mask,atmass,u)
  do iat=1,nat
    call dgemm('N','N',3,1,3,1.d0,u,3,tmp(1,iat),3,0.d0,xyz(1,iat),3)
  enddo
  write(11,'(I4/)') nat
  DO IAT=1,NAT
    call dscal(3,factr,xyz(1,iat),1)
    write(11,'(A3,3F15.10)') ztag(iat),(xyz(ixyz,iat),ixyz=1,3)
  ENDDO

CONTAINS

FUNCTION StrUpCase ( Input_String, ilen, iopt ) RESULT ( Output_String )
  ! -- Parameters
  CHARACTER(LEN=26), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(LEN=26), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  ! -- Argument and result
  INTEGER :: ilen,iopt
  CHARACTER(LEN=ilen), INTENT( IN ) :: Input_String
  CHARACTER(LEN=ilen)               :: Output_String
  ! -- Local variables
  INTEGER :: i, n, m

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
    ! -- Find location of letter in lower case constant string
    n = INDEX( LOWER_CASE, Output_String( i:i ) )
    m = INDEX( UPPER_CASE, Output_String( i:i ) )
    ! -- If current substring is a lower case letter, make it upper case
    IF ( n /= 0 ) THEN
      Output_String( i:i ) = UPPER_CASE( n:n )
    ELSE IF ( m == 0 .and. iopt == 1 ) then
      Output_String( i:i ) = ' '
    ENDIF
  END DO
END FUNCTION StrUpCase

FUNCTION convert_atlab_to_atmass ( catom ) RESULT ( xmass )

  IMPLICIT NONE
  ! -- Parameters
  INTEGER,PARAMETER ::I4B =SELECTED_INT_KIND(9)
  INTEGER,PARAMETER ::DP =KIND(1.0D0)
  real(dp), PARAMETER :: me=5.48579909065d-4
  
  integer(i4b), parameter :: iat = 20
  character(len=2), dimension(iat) :: cpertable = (/"H ","D ","T ","HE","LI","BE", &
                                                    "B ","C ","N ","O ","F ","NE", &
                                                    "NA","MG","AL","SI","P ","S ", &
                                                    "CL","AR" &
!                                                            ,"K ","CA","SC","TI", &
!                                                   "V ","CR","MN","FE","CO","NI", &
!                                                   "CU","ZN","GA","GE","AS","SE", &
!                                                   "BR","KR","RB","SR","Y ","ZR", &
!                                                   "NB","MO","TC","RU","RH","PD", &
!                                                   "AG","CD","IN","SN","SB","TE", &
!                                                   "I ","XE","CS","BA","HF","TA", &
!                                                   "W ","RE","OS","IR","PT","AU", &
!                                                   "HG","TL","PB","BI","PO","AT", &
!                                                   "RN" &
                                                    /)
  real(dp), dimension(iat) :: pertabmas = (/  1.00782503207d0,  2.014101778d0 ,  3.0160492777d0, & ! H,D,T
                                              4.00260325315d0,  7.01600455d0  ,  9.0121822d0   , & ! He,Li,Be
                                             11.0093054d0    , 12.d0          , 14.0030740048d0, & ! B,C,N
                                             15.99491461956d0, 18.99840322d0  , 19.9924401754d0, & ! O,F,Ne
                                             22.9897692809d0 , 23.985041700d0 , 26.98153863d0  , & ! Na,Mg,Al
                                             27.9769265325d0 , 30.97376163d0  , 31.97207100d0  , & ! Si,P,S
                                             34.96885268d0   , 39.9623831225d0                   & ! Cl,Ar,K
                                            /)
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: catom
  real(dp) :: xmass
  ! -- Local variables
  integer(i4b) :: i, n
 
  n = 0 
  do i=1,iat
    if(trim(catom) == trim(cpertable(i))) n=i
    if(n.ne.0) exit
  enddo
  if(n.eq.0) then
    write(6,'(A,A2)') "!!! ERROR !!! ELEMENT SYMBOL NOT RECOGNIZED: ",CATOM
    stop 
  endif
  xmass=pertabmas(n)
  return

END FUNCTION convert_atlab_to_atmass

SUBROUTINE transform_com(nat,mask,atmass,xyzeq)

  IMPLICIT NONE

  integer,          intent(in)    :: nat
  integer,          intent(in)    :: mask(nat)
  double precision, intent(in)    :: atmass(nat)
  double precision, intent(inout) :: xyzeq(3,nat)

  integer :: iat,ixyz

  double precision :: totmas
  double precision :: com(3)

  totmas=0.d0
  DO iat=1,nat
    if(mask(iat).ne.1) cycle
    totmas=totmas+atmass(iat)
  enddo
  com=0.d0
  do ixyz=1,3
    do iat=1,nat
      if(mask(iat).ne.1) cycle
      com(ixyz)=com(ixyz)+atmass(iat)*xyzeq(ixyz,iat)
    enddo
    com(ixyz)=com(ixyz)/totmas
  enddo
  do iat=1,NAT
    do ixyz=1,3
      xyzeq(ixyz,iat)=xyzeq(ixyz,iat)-com(ixyz)
    enddo
  enddo
  return

END SUBROUTINE transform_com

SUBROUTINE transform_stdorient(nat,coord,xabc,irep,map)
  implicit none

  integer :: nat
  integer :: irep
  integer, dimension(3,6) :: map
  double precision, dimension(3,nat) :: coord
  double precision, dimension(3,3)   :: xabc
  double precision, dimension(3,nat) :: tmp

  integer :: iat,iabc,ixyz
  double precision, dimension(3) :: xyz

  tmp=0.d0
  do iat=1,nat
    xyz=0.d0
    do iabc=1,3
      do ixyz=1,3
        xyz(map(iabc,irep))=xyz(map(iabc,irep))+coord(ixyz,iat)*xabc(ixyz,iabc)
      enddo
    enddo
    tmp(1:3,iat)=xyz(1:3)
  enddo
  coord=tmp
  return

end subroutine

SUBROUTINE calc_inertia(na,am,co,tm)
  implicit none

  integer :: na
  double precision, dimension(na)   :: am
  double precision, dimension(3,na) :: co
  double precision, dimension(3,3)  :: tm

  integer :: ia,ic,jc
  double precision xm
  double precision, dimension(3)    :: ac

  tm=0.d0
  do ia=1,na
    xm=am(ia)
    do ic=1,3
      ac(ic)=co(ic,ia)
    enddo
    tm(1,1) = tm(1,1) + (ac(2)*ac(2) + ac(3)*ac(3))*xm
    tm(2,1) = tm(2,1) - ac(1)*ac(2)*xm
    tm(2,2) = tm(2,2) + (ac(1)*ac(1) + ac(3)*ac(3))*xm
    tm(3,1) = tm(3,1) - ac(1)*ac(3)*xm
    tm(3,2) = tm(3,2) - ac(2)*ac(3)*xm
    tm(3,3) = tm(3,3) + (ac(1)*ac(1) + ac(2)*ac(2))*xm
  end do
  do ic=1,3
    do jc=1,ic
      tm(jc,ic)=tm(ic,jc)
    enddo
  enddo
  return
end subroutine 

SUBROUTINE calc_be(ie,be)

  implicit none

  double precision, PARAMETER :: FACTE=219474.63137054d0   

  integer :: na

  double precision, dimension(3) :: ie
  double precision, dimension(3) :: be

  integer :: iabc

  do iabc=1,3 
    be(iabc)=facte/(2.d0*ie(iabc))
  enddo
  return

end subroutine 

SUBROUTINE get_top(abc,itop)

  implicit none

  integer :: isig,itop
  double precision, dimension(3) :: abc

  double precision :: x
  character(len=7), dimension(2) :: ctop=(/"PROLATE"," OBLATE"/)

  x=(2*abc(2)-abc(1)-abc(3))/(abc(1)-abc(3))
  if(abs(x).gt.0.5d0) then
    isig=sign(1.d0,x)
    select case (isig)
      case (-1)
        itop=1  ! PROLATE, USING Ir
!       write(6,'(/1x,"NEAR ",A7," TOP   x =",F10.6)') ctop(1),x
      case (+1)
        itop=3  ! OBLATE   USING IIIR
!       write(6,'(/1x,"NEAR ",A7," TOP   x =",F10.6)') ctop(2),x
    end select
  else
    itop=2      !          USING IIR
!   write(6,'(/1x,"STRONG ASYMMETRY  x =",F10.6)') x
  endif
  return

END SUBROUTINE

SUBROUTINE rotaxis_phase(xabc)

  implicit none
  
  double precision, intent(inout) :: xabc(3,3)

  integer :: iabc
  integer :: ixyz,jxyz

  double precision :: su

  ! PHASE CONVENTION FOR ROTATIONAL AXES
  abcloop: do iabc=1,3
    ixyz=0
    do 
      ixyz=ixyz+1
      su=xabc(ixyz,iabc)
      if(abs(su).gt.1.d-10) then
        if(su.gt.0.d0) cycle abcloop 
        do jxyz=1,3
          su=xabc(jxyz,iabc)
          if(abs(su).gt.1.d-10) xabc(jxyz,iabc)=-su
        enddo
        cycle abcloop
      endif
    enddo
  enddo abcloop

END SUBROUTINE

subroutine quaternion(nat,xyz,xyzref,mask,atmass,u)

  implicit none

  integer, intent(in) :: nat

  double precision, intent(in)  :: xyz(3,nat)
  double precision, intent(in)  :: xyzref(3,nat)
  integer,          intent(in)  :: mask(nat)
  double precision, intent(in)  :: atmass(nat)

  double precision, intent(inout) :: u(3,3)

  double precision :: xmass
  double precision, dimension(nat) :: xm,xp,ym,yp,zm,zp

  double precision :: c(4,4),q(4),xk(4,4)

  integer :: iat,iq,jq

  xm=0.d0
  xp=0.d0
  ym=0.d0
  yp=0.d0
  zm=0.d0
  zp=0.d0
  do iat=1,nat
    if(mask(iat).ne.1) cycle
    xm(iat)=xyzref(1,iat)-xyz(1,iat)
    xp(iat)=xyzref(1,iat)+xyz(1,iat)
    ym(iat)=xyzref(2,iat)-xyz(2,iat)
    yp(iat)=xyzref(2,iat)+xyz(2,iat)
    zm(iat)=xyzref(3,iat)-xyz(3,iat)
    zp(iat)=xyzref(3,iat)+xyz(3,iat)
  enddo
  c=0.d0
  do iat=1,nat
    xmass=atmass(iat)
    c(1,1) = c(1,1) + xmass*(xm(iat)*xm(iat)+ym(iat)*ym(iat)+zm(iat)*zm(iat))
    c(2,1) = c(2,1) + xmass*(yp(iat)*zm(iat)-ym(iat)*zp(iat))
    c(3,1) = c(3,1) + xmass*(xm(iat)*zp(iat)-xp(iat)*zm(iat))
    c(4,1) = c(4,1) + xmass*(xp(iat)*ym(iat)-xm(iat)*yp(iat))
    c(2,2) = c(2,2) + xmass*(xm(iat)*xm(iat)+yp(iat)*yp(iat)+zp(iat)*zp(iat))
    c(3,2) = c(3,2) + xmass*(xm(iat)*ym(iat)-xp(iat)*yp(iat))
    c(4,2) = c(4,2) + xmass*(xm(iat)*zm(iat)-xp(iat)*zp(iat))
    c(3,3) = c(3,3) + xmass*(xp(iat)*xp(iat)+ym(iat)*ym(iat)+zp(iat)*zp(iat))
    c(4,3) = c(4,3) + xmass*(ym(iat)*zm(iat)-yp(iat)*zp(iat))
    c(4,4) = c(4,4) + xmass*(xp(iat)*xp(iat)+yp(iat)*yp(iat)+zm(iat)*zm(iat))
  enddo
  do iq=1,4
    do jq=iq,4
      c(iq,jq)=c(jq,iq)
    enddo
  enddo
!
  call DIAG(4,4,c,q,xk)
  call dcopy(4,xk(1,1),1,q,1)
!
  u(1,1) = q(1)**2+q(2)**2-q(3)**2-q(4)**2
  u(2,1) = 2.d0 * (q(2)*q(3)-q(1)*q(4))
  u(3,1) = 2.d0 * (q(2)*q(4)+q(1)*q(3)) 
  u(1,2) = 2.d0 * (q(2)*q(3)+q(1)*q(4))
  u(2,2) = q(1)**2-q(2)**2+q(3)**2-q(4)**2
  u(3,2) = 2.d0 * (q(3)*q(4)-q(1)*q(2))!  
  u(1,3) = 2.d0 * (q(2)*q(4)-q(1)*q(3))
  u(2,3) = 2.d0 * (q(3)*q(4)+q(1)*q(2))!  
  u(3,3) = q(1)**2-q(2)**2-q(3)**2+q(4)**2
!  
  return

end subroutine

SUBROUTINE DIAG(M,N,A,D,X)
  
  implicit none

  integer, intent(in) :: m,n
  double precision, intent(inout) :: A(m,n)
  double precision, intent(inout) :: D(n)
  double precision, intent(inout) :: X(m,n)

  integer :: info
  integer :: lwork,liwork
  double precision, allocatable :: work(:)
  integer, allocatable :: iwork(:)

  X=A
  lwork=-1
  liwork=-1
  allocate(iwork(1),work(1))
  call dsyevd('V','U',n,X,n,D,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(0,*) "ERROR IN DIAG dsyev work INFO=",info
    stop
  endif
  lwork=int(work(1))
  liwork=iwork(1)
  deallocate(work,iwork)
  allocate(work(lwork),iwork(liwork));
  call dsyevd('V','U',n,X,n,D,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(0,*) "ERROR IN DIAG dsyev diag INFO=",info
    stop
  endif
  deallocate(work,iwork)
  return

END SUBROUTINE DIAG

SUBROUTINE convert_list2mask(nat,list,mask)

  implicit none

  integer,      intent(in)  :: nat
  character(*), intent(in)  :: list
  integer,      intent(out) :: mask(nat)

  character(len=5) :: cnum

  integer :: nlen,ipos
  integer :: ic,ir
  integer :: iat,jat,kat
  integer :: inum
  integer :: istat

  mask=0
  nlen=len(trim(list))
  write(0,*) list, nlen
  cnum=""
  do ipos=1,nlen
    select case(list(ipos:ipos))
      case("-")
        ! collected lower bound
        ir=1
        read(cnum,'(I5)') iat
        cnum=""
      case(",")
        if(ir.eq.1) then
          !collected upper bound
          ir=0
          ic=0
          read(cnum,'(I5)') jat
          if(iat.gt.jat) call mask_error(1,iat,jat) 
          do kat=iat,jat
            mask(kat)=1
          enddo
          cnum=""
        else
          ! single index
          ic=0
          read(cnum,'(I5)') iat
          mask(iat)=1
          cnum=""
        endif
      case default 
        read(list(ipos:ipos),'(I1)',iostat=istat) inum
        if(istat.ne.0) call mask_error(2,ipos,0)
        if(inum.eq.0) cycle 
        ic=1
        write(cnum,'(A,I1)') trim(adjustl(cnum)),inum
    end select
  enddo
  if(ic.eq.1) then
    if(ir.eq.1) then
      read(cnum,'(I5)') jat
      if(iat.gt.jat) call mask_error(1,iat,jat) 
      do kat=iat,jat
        mask(kat)=1
      enddo
    else
      read(cnum,'(I5)') iat
      mask(iat)=1
      cnum=""
    endif
  endif

END SUBROUTINE

subroutine file_error(ifile,cfile)

  implicit none

  integer,      intent(in) :: ifile
  character(*), intent(in) :: cfile

  write(0,'("!!!ERROR!!! File ",I1," not found! Filename: ",a)') ifile,trim(cfile)
  stop

end subroutine

subroutine nat_error(iat,jat)

  implicit none

  integer, intent(in) :: iat,jat

  write(0,'("!!!ERROR!!! Number of atoms not consistent!",2I4)') iat,jat
  stop

end subroutine

subroutine mass_error(iat,ami,amj)

  implicit none

  integer,          intent(in) :: iat
  double precision, intent(in) :: ami,amj
  write(0,'("!!!ERROR!!! Inconsistent ordering! IAT: ",I4," Mi:",F10.0," Mj:",F10.0)') iat,ami,amj
  stop

end subroutine

subroutine mask_error(ierr,i,j)

  implicit none

  integer, intent(in) :: ierr
  integer, intent(in) :: i,j
  select case (ierr)
    case (1)
      write(0,'("!!!ERROR!!! Inconsistent index range! LOWER: ",I5," UPPER:",I5)') i,j
    case (2)
      write(0,'("!!!ERROR!!! Illegal character in framgment definition! POS:",I5)') i
  end select
  stop

end subroutine

END
