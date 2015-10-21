include ../make.inc

.SUFFIXES:.f90 .f

all: allmod lib


allmod: quadpack_qagi.o quadpack_qags.o quadpack_qawf.o quadpack_qaws.o quadpack_qag.o  quadpack_qagp.o quadpack_qawc.o quadpack_qawo.o quadpack_qng.o quadpack_aux.o MOD_QUADPACK.o

lib:
	ar cvq $(LIBQUADPACK) `ls *.o | sort | uniq`  
	${RANLIB} $(LIBQUADPACK)
	@echo "$(LIBQUADPACK) is done..."
	@echo ""
	rsync -avP *.o $(OBJ_INSTALL)/
	rsync -avP $(INC_INSTALL)/*.mod $(INC_TARGET)/

.f90.o:	
	$(FC) -c $(FFLAG) -I$(INC_INSTALL)/ $(MOPT)$(INC_INSTALL)/ $<

.f.o:	
	$(FC) -c $(FFLAG) $< 

clean:
	rm -fv *.o 
