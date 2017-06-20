# SciFortran

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

###Installation 
Installation is now available using CMake.

Clone the repo:
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
* environment module file `~/.modules.d/scifor/gnu`  
* homebrew `bash` script `SF_PREFIX/bin/configvars.sh`
 

The `CMake` compilation can be tuned using the following additional variables:   

* `-DPREFIX=prefix directory <~/opt/scifor/SVERSION/PLAT>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright (C) 2012-2017  Adriano Amaricci

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