<!-- ![BSA logo](./resources/images/BSA_logo_extended.PNG "BSA logo") -->

# BSA

BSA library provides a framework for the Bispectral Stochastic Analysis (BSA) of linear systems, 
under non-Gaussian stationary random excitations.

It is conceived to be a **plug-in** library for any personal/commercial Finite Element software capable 
of providing some basic features.

Its main API is defined in the main library interface module file `.\src\BsaLib\BsaLib.f90`.

The repository also provides a sample program `.\src\bsa.f90` which wants to help showing how the 
library interface could be integrated in your own program/library. 
This file could be adapted to the specific needs, if an executable is desired instead.

Currently, `BSA` only supports out-of-the box Windows build system (VS solution file).
Any contribution aimed at providing (also) a cross-platform build system (cmake, make-files, etc.) 
is welcome.
