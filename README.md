# SciFortran

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

### Dependencies

* gfortran > 4.9x **OR** ifort  > 13.0
* cmake > 3.0.0    
* lapack  ( https://github.com/aamaricci/Lapack )   
* blas  ( https://github.com/aamaricci/Blas )   
* MPI ( https://github.com/open-mpi/ompi )  [optional, recommended]
* scalapack  ( https://github.com/aamaricci/scalapack )  [optional]

If libraries are not available in your system, please use the provided links to install them. All libraries listed can be installed using `CMake`. Should Lapack/Blas not be available in your system, `scifor` will compile internal copies of such libraries. This option, however, in most cases cause a slight degradation of the performances with respect to optimized versions of the same libraries. Intel MKL support is offered using a custom `CMake` macro, contained in `cmake/FindMKL.cmake`, which however should be considered in a beta development.       



### Installation

Installation is now available using CMake. Experimental support for Intel MKL is provided but this still not universal and may end up in wrong linking directives. 

Clone the Scifor repo:

`git clone https://github.com/aamaricci/SciFortran scifor`

And from the repository directory (`cd scifor`) make a standard out-of-source CMake compilation:

`mkdir build`  
`cd build`  
`cmake ..`      (*)  
`make`     
`make install`   

(*) *In some cases CMake fails to find the MPI Fortran compiler, even if it is effectively installed and loaded into the system. An easy fix is to setup and export the `FC=mpif90` environment variable before invoking `cmake`.* 

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/scifor.pc`  
* environment module file `~/.modules.d/scifor/<PLAT>/<VERSION>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/scifor>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  

### UNINSTALL

`Cmake` does not officially provide uninstall procedure in the generate Makefile. Yet, it is possible to construct one such uninstall mode in a simple way. SCIFOR provides a way to uninstall the files generated inside any out-of-source compilation calling: 
`make uninstall`  



### PROBLEMS

`SciFortran` has been tested with success on several Unix/Linux platoforms. Support for Windows using Linux Bash Shell is experimental, although few people reported successful installation with minimal efforts. 

Some issues has reported concerning the wrong setup for the library `pkg-config` file, contained in  `$PREFIX/<PLAT>/<VERSION>/etc/scifor.pc`. The variable `Libs=-L${libdir} -lscifor <blas/lapack/scalapack>` produced by `cmake` during the configuration and installation process can be not properly defined for the part corresponding to third parties libraries such as Blas/Lapack/Scalapack. This breaks compilation against `scifor` using `pkg-config` produced linking options. 

FIX: edit the `scifor.pc` file manually, fixing the definition of the variable `Libs`. 

 

### CONTACT

For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright (C) Adriano Amaricci, Lorenzo Crippa, Giacomo Mazza, Massimo Capone

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LGPL for more details.

You should have received a copy of the GNU LGPL along with this program.  If not, see <http://www.gnu.org/licenses/>.

