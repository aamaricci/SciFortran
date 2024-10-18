Dependencies
============

**SciFortran** requires a number of third party software and, quite
obviously, a modern Fortran compiler. A restricted number of
this software is mandatory while another part can be considered
optional, though it is still strongly recommended to get maximum
performance and access all SciFortran features. 


* `GNU Fortran <https://gcc.gnu.org/fortran/>`_ > 5.0 **OR** `Intel Fortran Compiler <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html>`_  > 13.0
* `CMake <https://cmake.org/>`_ ≥ 3.5 
* `Make <https://www.gnu.org/software/make>`_ **OR** `Ninja
  <https://ninja-build.org/>`_ ≥ 1.10


Optional
-----------
* `MPI <https://github.com/open-mpi/ompi>`_  
* `LAPACK <https://github.com/Reference-LAPACK/lapack>`_  / `BLAS
  <https://netlib.org/blas>`_ [if not present in your system
  unoptimized [1]_ versions will be compiled]  
*  `ScaLAPACK <https://github.com/Reference-ScaLAPACK/scalapack>`_
   
If any of the required libraries is not available in your system, or
the version requirement is not satisfied, please install/upgrade
them. We generally advice for pre-packaged versions, as provided by
either
`apt <https://en.wikipedia.org/wiki/APT_(software)>`_,
`pip <https://pypi.org/project/pip/>`_,
`brew <https://formulae.brew.sh/>`_,
`conda <https://docs.conda.io/en/latest/>`_ or
`spack <https://spack.io/>`_.
The latter may provide the best options for HPC environments (trivial installation without sudo, easy management of interdependent multiple versions, automatic generation of environment modules, etc.), but choose freely according to your needs.


.. rubric:: Footnotes

.. [1] Should Lapack/Blas not be available in your system, SciFortran
       will compile internal copies of such libraries. This option,
       however, in most cases leads to noticeable performance
       degradation, with respect to optimized versions of the same
       libraries, such as
       `Intel's  MKL <https://en.wikipedia.org/wiki/Math_Kernel_Library>`_,
       `Apple's vecLib <https://developer.apple.com/documentation/accelerate/veclib>`_,
       `OpenBLAS <https://www.openblas.net/>`_, etc.
       Support for linking MKL is offered via a [custom CMake
       macro](./cmake/FindMKL.cmake), which should however be
       considered in beta development, as it is not universal and may
       end up in wrong linking directives. Apple's vecLib is instead
       automatically recognized by CMake, in any standard macOS
       system.
       



