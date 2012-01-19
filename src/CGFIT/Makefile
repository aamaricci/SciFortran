include $(HOME)/lib/opt.mk

ROUTINES=CGFIT.f90

all: libcgfit.a libcgfit_deb.a

libcgfit.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libcgfit_deb.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

# build: rout lib

# #-----------------------------------------------------------------------
# # Rules for compiling the .o files
# #-----------------------------------------------------------------------
# rout:  ${ROUTINES}
# 	@echo 'ROUTINES IS : '
# 	@echo ${ROUTINES}
# 	@echo 'MODULES ARE : '
# 	@echo ${MODULES}
# 	@rm -f *.mod *.o || echo "file first erased..."
# 	@rm -fv $(MYLIB)/libcgfit.a
# 	@echo "remove old library"
# 	ifort -c $(FFLAG) ${ROUTINES}
# 	@echo 'compiled' 

# #-----------------------------------------------------------------------
# # Target "lib" is used to built the subroutine library.
# #-----------------------------------------------------------------------

# lib: libdo $(SHARED_LIB)

# libdo : 
# 	ar ru libcgfit.a `ls *.o | sort | uniq `
# 	ranlib libcgfit.a
# 	@mv -v *.a $(MYLIB)/
# 	@mv -v *.mod $(MYINCLUDE)/
# 	@echo 'LIBRARY IS DONE.........'
# 	@rm *.o
# 	@echo ' '


# #-----------------------------------------------------------------------
# # Target "clean" is used to remove all the intermediate files.
# #-----------------------------------------------------------------------
# delfile: 
# 	@rm -f *.o






