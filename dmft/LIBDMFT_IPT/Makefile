#########################################################################
include $(SFDIR)/etc/lib.mk
#########################################################################

OBJS=IPT_VARS_GLOBAL.o \
IPT_MATS.o    \
IPT_KELDYSH.o \
IPT_SOPT.o    \
IPT_SC_SOPT.o \
IPT_SC_MATS.o \
IPT_AF_MATS.o \
DMFT_IPT.o



OBJS_DEB=IPT_VARS_GLOBAL_DEB.o \
IPT_MATS_DEB.o    \
IPT_KELDYSH_DEB.o \
IPT_SOPT_DEB.o    \
IPT_SC_SOPT_DEB.o \
IPT_SC_MATS_DEB.o \
IPT_AF_MATS_DEB.o \
DMFT_IPT_DEB.o


all: autoclean libdmft_ipt.a libdmft_ipt_deb.a

autoclean: 
	@rm -fv *.o *.mod $(SFLIB)/libdfmt_ipt*.a

libdmft_ipt.a: 	$(OBJS)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '	


libdmft_ipt_deb.a: $(OBJS_DEB)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@	  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


IPT_VARS_GLOBAL.o: IPT_VARS_GLOBAL.f90
	$(FC) $(FFLAG) -c IPT_VARS_GLOBAL.f90 $(SFMODS) 
IPT_KELDYSH.o: IPT_KELDYSH.f90
	$(FC) $(FFLAG) -c IPT_KELDYSH.f90 $(SFMODS) 
IPT_MATS.o: IPT_MATS.f90
	$(FC) $(FFLAG) -c IPT_MATS.f90 $(SFMODS)
IPT_SOPT.o: IPT_SOPT.f90
	$(FC) $(FFLAG) -c IPT_SOPT.f90 $(SFMODS)
IPT_ZEROT.o: IPT_ZEROT.f90
	$(FC) $(FFLAG) -c IPT_ZEROT.f90 $(SFMODS)
IPT_SC_SOPT.o: IPT_SC_SOPT.f90
	$(FC) $(FFLAG) -c IPT_SC_SOPT.f90 $(SFMODS)
IPT_SC_MATS.o: IPT_SC_MATS.f90
	$(FC) $(FFLAG) -c IPT_SC_MATS.f90 $(SFMODS)
IPT_AF_MATS.o: IPT_AF_MATS.f90
	$(FC) $(FFLAG) -c IPT_AF_MATS.f90 $(SFMODS)


IPT_VARS_GLOBAL_DEB.o: IPT_VARS_GLOBAL.f90
	$(FC) $(DFLAG) -c IPT_VARS_GLOBAL.f90 $(SFMODS_DEB) -o IPT_VARS_GLOBAL_DEB.o
IPT_KELDYSH_DEB.o: IPT_KELDYSH.f90
	$(FC) $(DFLAG) -c IPT_KELDYSH.f90 $(SFMODS_DEB) -o IPT_KELDYSH_DEB.o
IPT_MATS_DEB.o: IPT_MATS.f90
	$(FC) $(DFLAG) -c IPT_MATS.f90 $(SFMODS_DEB) -o IPT_MATS_DEB.o
IPT_SOPT_DEB.o: IPT_SOPT.f90
	$(FC) $(DFLAG) -c IPT_SOPT.f90 $(SFMODS_DEB) -o IPT_SOPT_DEB.o
IPT_ZEROT_DEB.o: IPT_ZEROT.f90
	$(FC) $(DFLAG) -c IPT_ZEROT.f90 $(SFMODS_DEB) -o IPT_ZEROT_DEB.o
IPT_SC_SOPT_DEB.o: IPT_SC_SOPT.f90
	$(FC) $(DFLAG) -c IPT_SC_SOPT.f90 $(SFMODS_DEB) -o  IPT_SC_SOPT_DEB.o
IPT_SC_MATS_DEB.o: IPT_SC_MATS.f90
	$(FC) $(DFLAG) -c IPT_SC_MATS.f90 $(SFMODS_DEB) -o  IPT_SC_MATS_DEB.o
IPT_AF_MATS_DEB.o: IPT_AF_MATS.f90
	$(FC) $(DFLAG) -c IPT_AF_MATS.f90 $(SFMODS_DEB) -o  IPT_AF_MATS_DEB.o

#
DMFT_IPT.o: DMFT_IPT.f90 
	$(FC) $(FFLAG) -c DMFT_IPT.f90 $(SFMODS)
DMFT_IPT_DEB.o: DMFT_IPT.f90
	$(FC) $(DFLAG) -c DMFT_IPT.f90 $(SFMODS_DEB) -o DMFT_IPT_DEB.o



clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
