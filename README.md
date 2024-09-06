# README

Here is a stripped down version of the interior strucuture solver STEREO-D. I just include an ideal equation of state and one test. 

## INSTALLATION

The program uses the **BOOST** (https://www.boost.org/) and **EIGEN** (https://eigen.tuxfamily.org) header libraries. Both are easy to install as they are simply header files. The user can either place them in a standard 'include' directory or specify their locations with the variable PDIR in 'makefile'. Either make a 'packages' directory, or change it to your own.

## OPERATION
To compile all executables run 

`make all`

or to compile individual executables run `make` followed by the executable name.

`make test` will compile and run the test (a static polytrope solution).

To run the executable type `./`+ the executable name
