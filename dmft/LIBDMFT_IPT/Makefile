include $(HOME)/lib/opt.mk

#########################################################################
include $(HOME)/lib/lib.mk
#########################################################################

OBJS=IPT_VARS_GLOBAL.o \
IPT_MATS.o    \
IPT_KELDYSH.o \
IPT_SOPT.o    \
IPT_SC_SOPT.o \
IPT_SC_MATS.o \
DMFT_IPT.o



OBJS_DEB=IPT_VARS_GLOBAL_DEB.o \
IPT_MATS_DEB.o    \
IPT_KELDYSH_DEB.o \
IPT_SOPT_DEB.o    \
IPT_SC_SOPT_DEB.o \
IPT_SC_MATS_DEB.o \
DMFT_IPT_DEB.o


all: autoclean libdmft_ipt.a libdmft_ipt_deb.a

autoclean: 
	@rm -fv *.o *.mod $(MYLIB)/libdfmt_ipt*.a

libdmft_ipt.a: 	$(OBJS)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '	


libdmft_ipt_deb.a: $(OBJS_DEB)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@	  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


IPT_VARS_GLOBAL.o: IPT_VARS_GLOBAL.f90
	$(FC) $(FFLAG) -c IPT_VARS_GLOBAL.f90 $(MODS) 
IPT_KELDYSH.o: IPT_KELDYSH.f90
	$(FC) $(FFLAG) -c IPT_KELDYSH.f90 $(MODS) 
IPT_MATS.o: IPT_MATS.f90
	$(FC) $(FFLAG) -c IPT_MATS.f90 $(MODS)
IPT_SOPT.o: IPT_SOPT.f90
	$(FC) $(FFLAG) -c IPT_SOPT.f90 $(MODS)
IPT_ZEROT.o: IPT_ZEROT.f90
	$(FC) $(FFLAG) -c IPT_ZEROT.f90 $(MODS)
IPT_SC_SOPT.o: IPT_SC_SOPT.f90
	$(FC) $(FFLAG) -c IPT_SC_SOPT.f90 $(MODS)
IPT_SC_MATS.o: IPT_SC_MATS.f90
	$(FC) $(FFLAG) -c IPT_SC_MATS.f90 $(MODS)
#
MPT_SOPT.o: MPT_SOPT.f90
	$(FC) $(FFLAG) -c MPT_SOPT.f90 $(MODS)
MPT_MATS.o: MPT_MATS.f90
	$(FC) $(FFLAG) -c MPT_MATS.f90 $(MODS)
MPT_SC_SOPT.o: MPT_SC_SOPT.f90
	$(FC) $(FFLAG) -c MPT_SC_SOPT.f90 $(MODS)
MPT_SC_MATS.o: MPT_SC_MATS.f90
	$(FC) $(FFLAG) -c MPT_SC_MATS.f90 $(MODS)
#
DMFT_IPT.o: DMFT_IPT.f90 
	$(FC) $(FFLAG) -c DMFT_IPT.f90 $(MODS)


IPT_VARS_GLOBAL_DEB.o: IPT_VARS_GLOBAL.f90
	$(FC) $(DFLAG) -c IPT_VARS_GLOBAL.f90 $(MODS_DEB) -o IPT_VARS_GLOBAL_DEB.o
IPT_KELDYSH_DEB.o: IPT_KELDYSH.f90
	$(FC) $(DFLAG) -c IPT_KELDYSH.f90 $(MODS_DEB) -o IPT_KELDYSH_DEB.o
IPT_MATS_DEB.o: IPT_MATS.f90
	$(FC) $(DFLAG) -c IPT_MATS.f90 $(MODS_DEB) -o IPT_MATS_DEB.o
IPT_SOPT_DEB.o: IPT_SOPT.f90
	$(FC) $(DFLAG) -c IPT_SOPT.f90 $(MODS_DEB) -o IPT_SOPT_DEB.o
IPT_ZEROT_DEB.o: IPT_ZEROT.f90
	$(FC) $(DFLAG) -c IPT_ZEROT.f90 $(MODS_DEB) -o IPT_ZEROT_DEB.o
IPT_SC_SOPT_DEB.o: IPT_SC_SOPT.f90
	$(FC) $(DFLAG) -c IPT_SC_SOPT.f90 $(MODS_DEB) -o  IPT_SC_SOPT_DEB.o
IPT_SC_MATS_DEB.o: IPT_SC_MATS.f90
	$(FC) $(DFLAG) -c IPT_SC_MATS.f90 $(MODS_DEB) -o  IPT_SC_MATS_DEB.o
#
MPT_SOPT_DEB.o: MPT_SOPT.f90
	$(FC) $(DFLAG) -c MPT_SOPT.f90 $(MODS_DEB) -o MPT_SOPT_DEB.o
MPT_MATS_DEB.o: MPT_MATS.f90
	$(FC) $(DFLAG) -c MPT_MATS.f90 $(MODS_DEB) -o MPT_MATS_DEB.o
MPT_SC_SOPT_DEB.o: MPT_SC_SOPT.f90
	$(FC) $(DFLAG) -c MPT_SC_SOPT.f90 $(MODS_DEB) -o  MPT_SC_SOPT_DEB.o
MPT_SC_MATS_DEB.o: MPT_SC_MATS.f90
	$(FC) $(DFLAG) -c MPT_SC_MATS.f90 $(MODS_DEB) -o  MPT_SC_MATS_DEB.o
#
DMFT_IPT_DEB.o: DMFT_IPT.f90
	$(FC) $(DFLAG) -c DMFT_IPT.f90 $(MODS_DEB) -o DMFT_IPT_DEB.o



clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
