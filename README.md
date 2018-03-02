# mimetic_evolve

Simple mimetic gravity background solver

Note on installation and compilation (on Mac):
- install gcc-7 with Homebrew
- compile using gcc-7 -O2 mimetic_evolve.c -o mimetic_evolve.exe
This is necessary due to some incompatibilities with inline functions on the default Mac compiler (clang).
