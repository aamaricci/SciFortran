SHELL  = /bin/bash 

RANLIB =  ranlib 
CUDA   =  nvcc
CXX    =  gcc -c
NVCC   =  nvcc -c 
LINK   =  gcc
NVCCFLAGS += --compiler-bindir=/usr/bin/gcc-4.3

#CUDA_INSTALL_PATH = $(CUDA) #/usr/local/cuda
#MKL_INSTALL_PATH  = $(MKLDIR) #/opt/mkl_10
INCLUDES          = -I./ -I$(MKLDIR)/include -I$(CUDADIR)/include

LIB               = -L$(CUDADIR)/lib64 -L$(MKLDIR)/lib/em64t -lmkl -lguide -lpthread -lcublas -lcudart -lcuda

PREP              = -fpp -D__CULA   

CUDAFLAG   =  -I/opt/cuda_sdk_2.3/common/inc -I./ -arch sm_13 -O3 -m64 --use_fast_math --compiler-bindir=/usr/bin/gcc-4.3
CCFLAGS    =  -O3 -I./  -static -traceback 
CXXFLAGS   =  -O3 -I./ -I/opt/cuda_sdk_2.3/common/inc -static  
FCFLAGS    =  -O3 -I./ ${PREP}   
#  CUDAFLAG   =  -I/opt/cuda_sdk_2.3/common/inc -I./ -arch sm_13 -O0  
#  CCFLAGS    =  -O0 -I./  -static -traceback
#  CXXFLAGS   =  -O0 -I./ -I/opt/cuda_sdk_2.3/common/inc -static -traceback 
#  FCFLAGS    =  -O0 -I./ -traceback ${PREP}   

CC         =  gcc
C++        =  gcc
FC         =  ifort

###############################
##       CUDA INSTALLED      ##
##         DO NOTHING        ##
###############################
CUFILES    =  cuda_init.cu cuda_Real_invert.cu cuda_Cinvert.cu cuda_Matmul.cu cuda_Matmul_complex.cu binding_to_cuda.cu Lanczos_real.cu Lanczos_complex.cu Lanczos_real_fly.cu Lanczos_complex_fly.cu Lanczos_real_updn.cu Lanczos_complex_updn.cu cuda_matrix_vector.cu cuda_array_inverse.cu cuda_array_matmul.cu  
CCFILES    =
FCFILES    =  fortran_cuda_routines.f90 fortran_cuda_library.f90
CXXFILES   =  
PUT_OBJ    =  cp CUDA_FORTRAN_BINDING/*.o ./  

###############################
##     CUDA NOT INSTALLED    ##
##  UNCOMMENT THE FOLLOWING  ##
###############################
#CUFILES    =
#CCFILES    =
#FCFILES    =  fortran_cuda_routines_no_cuda_installed.f90 fortran_cuda_library.f90
#CXXFILES   =  no_cuda_installed.c
#PUT_OBJ    =  

all:    
	@echo '*********** Finished compilation **************'
	@echo '       use make deltemp to remove bad files    '
	@echo '***********************************************'
rout:  
	@rm -f ./include/*.mod  || echo "no such files..."
	@rm -f *.mod *.o *.a || echo "no such files..."
	${CUDA}   -c ${CUDAFLAG} ${CUFILES} || echo "no files.."
	${CC}     -c ${CCFLAGS}  ${CCFILES} || echo "no files.."
	${FC}     -c ${FCFLAGS}  ${FCFILES} || echo "no files.."
	${C++}    -c ${CXXFLAGS} ${CXXFILES} || echo "no files.."
	@mv *.mod include/
	@rm -f *.linkinfo || echo "no such files..."
	@echo "done...."

build: rout libdo

libdo: 
	@echo "build library now...."
	@ echo ${PUT_OBJ}
	@ ${PUT_OBJ}
	ar ru libgpu.a `ls *.o | sort | uniq `
	${RANLIB} libgpu.a
	@mv -v *.a lib/
	@rm -f *.o || echo "no such files...."
	@echo 'LIBRARY IS DONE.........'

all : benchmark

benchmark : benchmark.o residual_norm.o random_matrix_t.o gpu_lapack_internal.o gpu_sgetrf.o gpu_sgeqrf.o gpu_spotrf.o gpu_blas3.o gpu_batch_sswap.o gpu_transpose.o gpu_misc.o
	$(LINK) $(LIB) -o benchmark benchmark.o residual_norm.o random_matrix_t.o gpu_lapack_internal.o gpu_sgetrf.o gpu_sgeqrf.o gpu_spotrf.o gpu_blas3.o gpu_batch_sswap.o gpu_transpose.o gpu_misc.o

benchmark.o : benchmark.cpp
	$(CXX) $(INCLUDES) benchmark.cpp

residual_norm.o : residual_norm.cpp
	$(CXX) $(INCLUDES) residual_norm.cpp

random_matrix_t.o : random_matrix_t.cpp
	$(CXX) $(INCLUDES) random_matrix_t.cpp

gpu_misc.o : gpu_misc.cu
	$(NVCC) $(INCLUDES) -o gpu_misc.o gpu_misc.cu

gpu_blas3.o : gpu_blas3.cu
	$(NVCC) $(INCLUDES) -o gpu_blas3.o gpu_blas3.cu

gpu_batch_sswap.o : gpu_batch_sswap.cu
	$(NVCC) $(INCLUDES) -o gpu_batch_sswap.o gpu_batch_sswap.cu

gpu_transpose.o : gpu_transpose.cu
	$(NVCC) $(INCLUDES) -o gpu_transpose.o gpu_transpose.cu

gpu_lapack_internal.o : gpu_lapack_internal.cpp
	$(CXX) $(INCLUDES) gpu_lapack_internal.cpp

gpu_sgetrf.o : gpu_sgetrf.cpp
	$(CXX) $(INCLUDES) gpu_sgetrf.cpp

gpu_sgeqrf.o : gpu_sgeqrf.cpp
	$(CXX) $(INCLUDES) gpu_sgeqrf.cpp

gpu_spotrf.o : gpu_spotrf.cpp
	$(CXX) $(INCLUDES) gpu_spotrf.cpp

clean :
	rm -f *.o *.linkinfo 
