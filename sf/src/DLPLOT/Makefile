include opt.mk

ROUTINES=DLPLOT.f90

# #Static compilation WANTS linking against X11-XORG libraries. 
# #PLEASE be sure to have them installed
# X11  = -L$(XLIB) -lX11 -lXt -lXext -lxcb -lX11-xcb -lXau -lXdmcp

# # #These variables are defined only to compile test program: 
# # #DISLIN=<path to DISLIN directory>
DSL      = -L${DISLIN} -ldislin
DLPLOT   = -L$(SFLIB) -ldlplot
TESTLIB  = ${DLPLOT} ${DSL} ${X11}
TESTMODS = -I$(HOME)/lib/src/MKL -I${DISLIN}/ifc -I$(SFINCLUDE)


#============COMPILE THE LIBRARY===============
# all: comp libdlplot.a
all: libdlplot.a libdlplot_deb.a

libdlplot.a:
	@rm -fv $(SFLIB)/$@
	$(FC) -c $(FFLAG) $(ROUTINES) 
	@mv -fv *.mod $(SFINCLUDE)/
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@rm -f *.o
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libdlplot_deb.a:
	@rm -fv $(SFLIB)/$@
	$(FC) -c $(DFLAG) $(ROUTINES)
	@mv -fv *.mod $(SFINCLUDE)/debug/
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

# #============COMPILE TEST PROGRAM===============
test:	
	@echo " ........... compile: TEST ........... "
	${FC} test${ROUTINES} -o ./test ${TESTMODS} ${TESTLIB}
	@mkdir TESTDIR
	@mv test TESTDIR
	@echo " ........... run TEST ........... "
	@cd TESTDIR;./test
	@echo "done............................."
