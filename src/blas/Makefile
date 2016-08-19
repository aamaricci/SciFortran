include ../../make.inc


SBLAS1 = isamax.o sasum.o saxpy.o scopy.o sdot.o snrm2.o srot.o srotg.o sscal.o sswap.o sdsdot.o srotmg.o srotm.o
CBLAS1 = scabs1.o scasum.o scnrm2.o icamax.o caxpy.o ccopy.o cdotc.o cdotu.o csscal.o crotg.o cscal.o cswap.o csrot.o
DBLAS1 = idamax.o dasum.o daxpy.o dcopy.o ddot.o dnrm2.o drot.o drotg.o dscal.o dsdot.o dswap.o drotmg.o drotm.o
ZBLAS1 = dcabs1.o dzasum.o dznrm2.o izamax.o zaxpy.o zcopy.o zdotc.o zdotu.o zdscal.o zrotg.o zscal.o zswap.o zdrot.o
CB1AUX = isamax.o sasum.o saxpy.o scopy.o snrm2.o sscal.o
ZB1AUX = idamax.o dasum.o daxpy.o dcopy.o dnrm2.o dscal.o
ALLBLAS  = lsame.o xerbla.o
SBLAS2 = sgemv.o sgbmv.o ssymv.o ssbmv.o sspmv.o strmv.o stbmv.o stpmv.o strsv.o stbsv.o stpsv.o sger.o ssyr.o sspr.o ssyr2.o sspr2.o
CBLAS2 = cgemv.o cgbmv.o chemv.o chbmv.o chpmv.o ctrmv.o ctbmv.o ctpmv.o ctrsv.o ctbsv.o ctpsv.o cgerc.o cgeru.o cher.o chpr.o cher2.o chpr2.o
DBLAS2 = dgemv.o dgbmv.o dsymv.o dsbmv.o dspmv.o dtrmv.o dtbmv.o dtpmv.o dtrsv.o dtbsv.o dtpsv.o dger.o dsyr.o dspr.o dsyr2.o dspr2.o
ZBLAS2 = zgemv.o zgbmv.o zhemv.o zhbmv.o zhpmv.o ztrmv.o ztbmv.o ztpmv.o ztrsv.o ztbsv.o ztpsv.o zgerc.o zgeru.o zher.o zhpr.o zher2.o zhpr2.o
SBLAS3 = sgemm.o ssymm.o ssyrk.o ssyr2k.o strmm.o strsm.o 
CBLAS3 = cgemm.o csymm.o csyrk.o csyr2k.o ctrmm.o ctrsm.o chemm.o cherk.o cher2k.o
DBLAS3 = dgemm.o dsymm.o dsyrk.o dsyr2k.o dtrmm.o dtrsm.o
ZBLAS3 = zgemm.o zsymm.o zsyrk.o zsyr2k.o ztrmm.o ztrsm.o zhemm.o zherk.o zher2k.o


SINGLE=$(SBLAS1) $(ALLBLAS) $(SBLAS2) $(SBLAS3)

DOUBLE= $(DBLAS1) $(ALLBLAS) $(DBLAS2) $(DBLAS3)

COMPLEX= $(CBLAS1) $(CB1AUX) $(ALLBLAS) $(CBLAS2) $(CBLAS3)

COMPLEX16= $(ZBLAS1) $(ZB1AUX) $(ALLBLAS) $(ZBLAS2) $(ZBLAS3)



.SUFFIXES:.f90 .f

#$(SINGLE) $(COMPLEX) 
lib: $(DOUBLE) $(COMPLEX16)

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
