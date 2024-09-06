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

## Code structure

The organisation of the code works as follows. The solution is held in C++ stuct called P_P (‘StructStruct.hpp’), which is passed around different parts of the code. An initial guess is set up (in my stripped version, this is in ‘Struct_test_boost.cc’, but the full code has more examples). This is passed to the class ‘timestepper’ in ‘timestepper.hpp’. In the full code this steps the evolution forward, in the stripped version, I just left a static structure version, because your own problem will have its particular way you will want to advance the evolution. ‘timestepper’ passes the solution on to the class ‘converge’ in ‘converger.hpp/cpp’. This then passes the solution to JacobMaker (‘Jacob_Maker.hpp/cpp’) which sets up the discretised equations and the jacobian matrix. Then ‘converge’ does the iteration steps. 
