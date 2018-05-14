# mimetic_evolve

Simple mimetic gravity background solver

Note on installation and compilation (on Mac):
- install gcc-8 with Homebrew
- compile using gcc-8 -O2 mimetic_evolve.c -o mimetic_evolve.exe
- run with ./mimetic_evolve.exe

Use of gcc from Homebrew is necessary to avoid some incompatibilities with inline functions on the default Mac compiler (clang), and for OpenMP to work properly.
