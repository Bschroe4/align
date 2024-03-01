# align
Alignment of xyz coordinates according to Eckart conditions.

Compile the code with `ifort` using the following commandline:
```
ifort -o align.x align.F90 -mkl
```

The program requires 3 arguments (file names) that are passed through the command line call:
```
./align.x struct1.xyz struct2.xyz2 struct3.xyz
```
where the given file names are placeholders and can be given however the user wants. The `struct1.xyz` and `struct2.xyz` files are required to be existent, otherwise the program aborts. Struct3.xyz is generated by the program and any existing file with the corresponding name is overwritten.

Possible extensions: 
- only use a fragment of the molecule for alignment
