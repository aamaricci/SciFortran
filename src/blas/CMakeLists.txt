set(SBLAS1 isamax.f sasum.f saxpy.f scopy.f sdot.f snrm2.f srot.f srotg.f sscal.f sswap.f sdsdot.f srotmg.f srotm.f)

set(CBLAS1 scabs1.f scasum.f scnrm2.f icamax.f caxpy.f ccopy.f cdotc.f cdotu.f csscal.f crotg.f cscal.f cswap.f csrot.f)
set(DBLAS1 idamax.f dasum.f daxpy.f dcopy.f ddot.f dnrm2.f drot.f drotg.f dscal.f dsdot.f dswap.f drotmg.f drotm.f)
set(ZBLAS1 dcabs1.f dzasum.f dznrm2.f izamax.f zaxpy.f zcopy.f zdotc.f zdotu.f zdscal.f zrotg.f zscal.f zswap.f zdrot.f)
set(CB1AUX isamax.f sasum.f saxpy.f scopy.f snrm2.f sscal.f)
set(ZB1AUX idamax.f dasum.f daxpy.f dcopy.f dnrm2.f dscal.f)
set(ALLBLAS  lsame.f xerbla.f)
set(SBLAS2 sgemv.f sgbmv.f ssymv.f ssbmv.f sspmv.f strmv.f stbmv.f stpmv.f strsv.f stbsv.f stpsv.f sger.f ssyr.f sspr.f ssyr2.f sspr2.f)
set(CBLAS2 cgemv.f cgbmv.f chemv.f chbmv.f chpmv.f ctrmv.f ctbmv.f ctpmv.f ctrsv.f ctbsv.f ctpsv.f cgerc.f cgeru.f cher.f chpr.f cher2.f chpr2.f)
set(DBLAS2 dgemv.f dgbmv.f dsymv.f dsbmv.f dspmv.f dtrmv.f dtbmv.f dtpmv.f dtrsv.f dtbsv.f dtpsv.f dger.f dsyr.f dspr.f dsyr2.f dspr2.f)
set(ZBLAS2 zgemv.f zgbmv.f zhemv.f zhbmv.f zhpmv.f ztrmv.f ztbmv.f ztpmv.f ztrsv.f ztbsv.f ztpsv.f zgerc.f zgeru.f zher.f zhpr.f zher2.f zhpr2.f)
set(SBLAS3 sgemm.f ssymm.f ssyrk.f ssyr2k.f strmm.f strsm.f)
set(CBLAS3 cgemm.f csymm.f csyrk.f csyr2k.f ctrmm.f ctrsm.f chemm.f cherk.f cher2k.f)
set(DBLAS3 dgemm.f dsymm.f dsyrk.f dsyr2k.f dtrmm.f dtrsm.f)
set(ZBLAS3 zgemm.f zsymm.f zsyrk.f zsyr2k.f ztrmm.f ztrsm.f zhemm.f zherk.f zher2k.f)

set(SINGLE  ${SBLAS1} ${ALLBLAS} ${SBLAS2} ${SBLAS3} )

set(DOUBLE ${DBLAS1} ${ALLBLAS} ${DBLAS2} ${DBLAS3} )

set(COMPLEX ${CBLAS1} ${CB1AUX} ${ALLBLAS} ${CBLAS2} ${CBLAS3} )

set(COMPLEX16 ${ZBLAS1} ${ZB1AUX} ${ALLBLAS} ${ZBLAS2} ${ZBLAS3} )

set(BLAS_LOCAL_SRC ${DOUBLE} ${COMPLEX16})


SET(BLAS_src ${BLAS_LOCAL_SRC})
ADD_LIBRARY(BLASLIB OBJECT ${BLAS_src})
