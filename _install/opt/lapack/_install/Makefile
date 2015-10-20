#
#  Top Level Makefile for LAPACK
#  Version 3.2.2
#  June 2010
#

include make.inc

all: lapack_install lib 
#lapack_install lapack_testing blas_testing

lapack_install:
	( cd install; $(MAKE)) 
#./testlsame; ./testslamch; ./testdlamch; ./testsecond; ./testdsecnd; ./testversion )

lib: lapacklib
#lib: blaslib variants lapacklib tmglib

lapacklib: #lapack_install
	( cd src; $(MAKE) )


clean: cleanlib #cleantesting cleanblas_testing 


blaslib:
	( cd blas/src; $(MAKE) )


variants:
	( cd src/VARIANTS ; $(MAKE))

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	( cd TESTING ; $(MAKE) )

variants_testing: lib variants
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/cholrl.a' ; \
	mv stest.out stest_cholrl.out ; mv dtest.out dtest_cholrl.out ; mv ctest.out ctest_cholrl.out ; mv ztest.out ztest_cholrl.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/choltop.a' ; \
	mv stest.out stest_choltop.out ; mv dtest.out dtest_choltop.out ; mv ctest.out ctest_choltop.out ; mv ztest.out ztest_choltop.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/lucr.a' ; \
	mv stest.out stest_lucr.out ; mv dtest.out dtest_lucr.out ; mv ctest.out ctest_lucr.out ; mv ztest.out ztest_lucr.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/lull.a' ; \
	mv stest.out stest_lull.out ; mv dtest.out dtest_lull.out ; mv ctest.out ctest_lull.out ; mv ztest.out ztest_lull.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/lurec.a' ; \
	mv stest.out stest_lurec.out ; mv dtest.out dtest_lurec.out ; mv ctest.out ctest_lurec.out ; mv ztest.out ztest_lurec.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='src/VARIANTS/LIB/qrll.a' ; \
	mv stest.out stest_qrll.out ; mv dtest.out dtest_qrll.out ; mv ctest.out ctest_qrll.out ; mv ztest.out ztest_qrll.out )

blas_testing:
	( cd blas/TESTING; $(MAKE) -f Makeblat1 )
	( cd blas; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd blas/TESTING; $(MAKE) -f Makeblat2 )
	( cd blas; ./xblat2s < sblat2.in     ; \
	           ./xblat2d < dblat2.in     ; \
	           ./xblat2c < cblat2.in     ; \
	           ./xblat2z < zblat2.in     )
	( cd blas/TESTING; $(MAKE) -f Makeblat3 )
	( cd blas; ./xblat3s < sblat3.in     ; \
	           ./xblat3d < dblat3.in     ; \
	           ./xblat3c < cblat3.in     ; \
	           ./xblat3z < zblat3.in     ) 

cleanlib:
	( cd install; $(MAKE) clean )
	( cd blas/src; $(MAKE) clean )
	( cd src; $(MAKE) clean )
	( cd src/VARIANTS; $(MAKE) clean )

cleanblas_testing:	
	( cd blas/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd blas/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd blas/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd blas; rm -f xblat* )

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm -f xlin* xeig* )

cleanall: cleanlib cleanblas_testing cleantesting
	rm -rfv src/obj_*
	rm -fv install/test*  blas/*.out

