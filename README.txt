PROJECT: sbakker3
FILES: pattern.c
AUTHOR: Sasha Bakker
RELEASE DATE: 10/14/21
COMPILER: MinGW GCC
IDE: ECLIPSE
OS: Windows 10

DESCRIPTION
-----------
The purpose of this project is to solve a system of two discretized coupled partial differential equations on a square grid for a fixed amount of time.



RUNNING THE PROGRAM
-------------------
The program requires the header files `stdio.h`, `math.h`, and `string.h` from the standard C library. The file `pattern.c` must be complied and run in the terminal.

How open the terminal in Eclipse (make sure version is newest):

1. Right click on the name of the project in the project explorer
2. Left click "show in local terminal"
3. Left click "terminal"

An alternative method is to hold Ctrl+Alt+T in Windows.

How to compile the program: "gcc -o pattern pattern.c".

How to run program: "filename" for reflective boundary conditions, "pattern p" for periodic boundary conditions.
