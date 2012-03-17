include $(SFDIR)/etc/opt.mk

ROUTINES=d_ordered_list.f90 d_unordered_list.f90 z_unordered_list.f90 

all: liblist.a liblist_deb.a

liblist.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

liblist_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


