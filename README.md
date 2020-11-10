# SciFortran

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

### Dependencies

* gfortran > 4.9x **OR** ifort  > 13.0
* cmake > 2.8.8    
* lapack  ( https://github.com/aamaricci/Lapack )   
* blas  ( https://github.com/aamaricci/Blas )   
* MPI ( https://github.com/open-mpi/ompi )  [optional, recommended]
* scalapack  ( https://github.com/aamaricci/scalapack )  [optional]

If libraries are not available in your system, please use the provided links to install them. All libraries listed can be installed using `CMake` 



### Installation

Installation is now available using CMake. Experimental support for Intel MKL is provided but this still not universal and may end up in wrong linking directives. 

Clone the Scifor repo:

`git clone https://github.com/aamaricci/SciFortran scifor`

And from the repository directory (`cd scifor`) make a standard out-of-source CMake compilation:

`mkdir build`
`cd build`
`cmake ..`     
`make`     
`make install`   
`make post-install`    

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/scifor.pc`  
* environment module file `~/.modules.d/scifor/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/scifor/VERSION/PLAT>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright (C) Adriano Amaricci, Lorenzo Crippa, Giacomo Mazza, Massimo Capone

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.