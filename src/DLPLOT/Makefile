include $(HOME)/lib/opt.mk

ROUTINES=DLPLOT.f90

# #Static compilation WANTS linking against X11-XORG libraries. 
# #PLEASE be sure to have them installed
# X11  = -L$(XLIB) -lX11 -lXt -lXext -lxcb -lX11-xcb -lXau -lXdmcp

# #These variables are defined only to compile test program: 
# #DISLIN=<path to DISLIN directory>
DSL      = -L${DISLIN} -ldislin
DLPLOT   = -L$(MYLIB) -ldlplot
TESTLIB  = ${DLPLOT} ${DSL} ${X11}
TESTMODS = -I$(HOME)/lib/src/MKL -I${DISLIN}/ifc -I$(MYINCLUDE)

# #Creation of this library needs to Include DISLIN object files and mods. 
# #If you use Intel Fortran Compiler 
# #to compile this is achieved with the following line, where:
# #DISLIN variable is a system wide variable containing the location of  DISLIN directory:
# #e.g. DISLIN=/usr/local/dislin or DISLIN=/opt/dislin-10.0
# #if you don't have this variable in your environment (.bashrc or .bash_profile) define it here
# LIBS=-I${DISLIN}/ifc ${X11}


#============COMPILE THE LIBRARY===============
# all: comp libdlplot.a
all: libdlplot.a libdlplot_deb.a

libdlplot.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libdlplot_deb.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

#============COMPILE TEST PROGRAM===============
test:	
	@echo " ........... compile: TEST ........... "
	${FC} test${ROUTINES} -o ./test ${TESTMODS} ${TESTLIB}
	@mkdir TESTDIR
	@mv test TESTDIR
	@echo " ........... run TEST ........... "
	@cd TESTDIR;./test
	@echo "done............................."
