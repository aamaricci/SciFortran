# SciFortran

[![Ubuntu](https://img.shields.io/github/actions/workflow/status/QcmPlab/SciFortran/Ubuntu_Scheduled.yml?label=Ubuntu&logo=ubuntu&style=flat-square)](https://github.com/QcmPlab/SciFortran/actions/workflows/Scheduled.yml) 
[![MacOS](https://img.shields.io/github/actions/workflow/status/QcmPlab/SciFortran/MacOS_Scheduled.yml?label=macOS&logo=apple&style=flat-square)](https://github.com/QcmPlab/SciFortran/actions/workflows/Scheduled.yml) 
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/SciFortran)

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

### Dependencies

* [GNU Fortran (`gfortran`)](https://gcc.gnu.org/fortran/) > 5.0 **OR** [Intel Fortran Compiler Classic (`ifort`)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)  > 13.0
* [CMake](https://cmake.org/) ≥ 3.0 [> 3.16 for ninja support] 
* [Make](https://www.gnu.org/software/make/) **OR** [Ninja](https://ninja-build.org/) ≥ 1.10 
* [LAPACK](https://github.com/Reference-LAPACK/lapack) / [BLAS](https://netlib.org/blas/) [if not present in your system unoptimized[^1] versions will be compiled]  
* [MPI](https://github.com/open-mpi/ompi)  [optional, recommended]
* [ScaLAPACK](https://github.com/Reference-ScaLAPACK/scalapack)  [optional, recommended]

If any of the required libraries is not available in your system, or the version requirement is not satisfied, please install/upgrade them. We generally advice for pre-packaged versions, as provided by either [`apt`](https://en.wikipedia.org/wiki/APT_(software)), [`pip`](https://pypi.org/project/pip/), [`brew`](https://formulae.brew.sh/), [`conda`](https://docs.conda.io/en/latest/) or [`spack`](https://spack.io/). The latter may provide the best options for HPC environments (trivial installation without sudo, easy management of interdependent multiple versions, automatic generation of environment modules, etc.), but choose freely according to your needs.

[^1]: Should Lapack/Blas not be available in your system, SciFortran will compile internal copies of such libraries. This option, however, in most cases leads to noticeable performance degradation, with respect to optimized versions of the same libraries, such as [Intel's MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library), [Apple's vecLib](https://developer.apple.com/documentation/accelerate/veclib), [OpenBLAS](https://www.openblas.net/), etc. Support for linking MKL is offered via a [custom CMake macro](./cmake/FindMKL.cmake), which should however be considered in beta development, as it is not universal and may end up in wrong linking directives. Apple's vecLib is instead automatically recognized by CMake, in any standard macOS system.



## BUILD

Our build system relies on CMake. Experimental support for linking Intel's MKL is provided, although it may fail on some systems.

Clone the SciFortran repo:

```
git clone https://github.com/aamaricci/SciFortran scifor
```

Optionally[^2] define the fortran compiler:

```
export FC=mpif90 # or gfortran or ifort if you don't need MPI
```

From the repository directory (`cd scifor`) make a standard out-of-source CMake compilation:

<details>
<summary> Using <tt>make</tt> (click to expand) </summary>
Default CMake workflow, with widest version support (CMake > 3.0).

```
mkdir build 
cd build  
cmake .. 
make
```      

</details>

<details>
<summary> Using <tt>ninja</tt> (click to expand)</summary>

If a fortran-capable[^3] version of `ninja` ( https://ninja-build.org ) is available in your system (and CMake can[^4] take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

```
mkdir build    
cd build  
cmake -GNinja ..  
ninja
```       

</details>

In both cases you can pass additional options to the `cmake` command, to customize the build (default values between `< >`):

* `-DPREFIX=prefix/directory <~/opt/scifor>`, to specify a custom installation directory 

* `-DUSE_MPI=<yes>/no`, to (de)activate MPI support  

* `-DVERBOSE=yes/<no>`, to generate verbose makefiles

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`, to select predefined compiler options

* `-DWITH_BLAS_LAPACK=yes/<no>`, to skip search of preinstalled linear algebra libraries and enforce compilation from source

* `-DWITH_SCALAPACK=<yes>/no`, to link with preinstalled ScaLAPACK library, for parallel algebra support.

[^2]: In some cases CMake fails to find the MPI fortran compiler, even if it is effectively installed and loaded into the system. An easy fix is to setup and export the `FC=mpif90` environment variable before invoking the `cmake <options> ..` command. 

[^3]: Ninja did not support fortran before version 1.10, although Kitware has long mantained a fortran-capable fork, which might be obtained easily as a [Spack package](https://packages.spack.io/package.html?name=ninja-fortran). Nevertheless we note that as of fall 2022 `pip install ninja --user` [ships Ninja v1.10.2](https://pypi.org/project/ninja/), hence obtaining a suitable official Ninja release should be trivial.

[^4]: This depends on your CMake version. Comparing [this](https://cmake.org/cmake/help/v3.16/generator/Ninja.html#fortran-support) to [this](https://cmake.org/cmake/help/v3.17/generator/Ninja.html#fortran-support) would suggest that CMake started supporting Ninja's fortran features only after v3.17 but we have verified that at least v3.16.3 (current version shipped by `apt` on Ubuntu 20.04 LTS) does indeed work. For more information you can take a look to the [related issue](https://github.com/QcmPlab/SciFortran/issues/16). 

## INSTALL

System-wide installation is completed after the build step using either: 

```
make install
```  

or   

```
ninja install
```  
 
To actually link the library we provide some alternatives:

* A generated [pkg-config](https://github.com/freedesktop/pkg-config) file to, installed to `~/.pkgconfig.d/scifor.pc`  
* A generated [environment module](https://github.com/cea-hpc/modules), installed to `~/.modules.d/scifor/<PLAT>/<VERSION>`  
* A generated `bash` script at `<PREFIX>/bin/configvars.sh`, to be sourced for permanent loading.

which you can choose among by following the instructions printed on screen.

## UNINSTALL

CMake does not officially provide uninstall procedures in the generated Make/Ninja files. Hence SciFortran supplies a homebrew method to remove the generated files by calling (from the relevant build folder): `make uninstall` / `ninja uninstall`.

## DOCUMENTATION
The repository is configured for source-based, automated, API docs generation via [FORD](https://github.com/Fortran-FOSS-Programmers/ford), so that you can build the html files for your specific local version by running
```
ford docs.config
```
at the root directory of the project. Alternatively you can examine the [official docs](https://qcmplab.github.io/SciFortran) for the latest commit on master, as generated by our continuous deployment workflow.

## KNOWN PROBLEMS

`SciFortran` has been tested with success on several Unix/Linux platforms. Support for Windows, through [WSL](https://learn.microsoft.com/en-us/windows/wsl/install), is still experimental, although few people reported successful installation with minimal efforts. 

Some have reported issues concerning the wrong setup for the library `pkg-config` file, contained in  `$PREFIX/<PLAT>/<VERSION>/etc/scifor.pc`. The variable `Libs=-L${libdir} -lscifor <blas/lapack/scalapack>` produced by `cmake` during the configuration and installation process can be not properly defined for the part corresponding to third parties libraries such as Blas/Lapack/Scalapack. This breaks compilation against `scifor` whenever `pkg-config` is used to generate the linking options. 

> FIX: edit the `scifor.pc` file manually, overwriting the definition of the variable `Libs`, as appropriate for your system. 

 

### CONTACT

If you encounter bugs or difficulties, please [file an issue](https://github.com/QcmPlab/SciFortran/issues/new/choose). For any other communication, please reach out to:    
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright (C) Adriano Amaricci, Lorenzo Crippa, Gabriele Bellomia, Giacomo Mazza, Samuele Giuli, Massimo Capone

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LGPL for more details.

You should have received a copy of the GNU LGPL along with this program.  If not, see <http://www.gnu.org/licenses/>.

