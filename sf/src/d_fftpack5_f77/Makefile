# Several changes made by HCP so this would build without trouble 
# on a Linux/g77 system.
# (1) changed step to build library to use ar instead of update
#     (update must mean something different on someone elses Unix.)
# (2) Added make clean step 
# (3) In test step, changed a.out to ./a.out for cautious folk who don't 
#     have "." in their PATH.
# (4) Change FFLAGS from -O to -O2 -funroll-loops
# (5) Specify FC=gcc in case /usr/bin/f77 is not a link to g77 
#     (as it won't be if you have f77reorder installed)
# (6) Added targets shared and installshared to make and install a shared 
#     version of the library. You need /usr/local/lib in /etc/ld.so.conf 
#     for this to work   
# (7) Modified names for dble prec version
LIB=dfftpack

# Use these lines for Linux/g77
FC=ifort
FFLAGS=-O2

# Use these lines for Solaris
#FC=f77
#FFLAGS=-fast -O5

OBJ=\
zfftb.o\
cfftb1.o\
zfftf.o\
cfftf1.o\
zffti.o\
cffti1.o\
dcosqb.o\
cosqb1.o\
dcosqf.o\
cosqf1.o\
dcosqi.o\
dcost.o\
dcosti.o\
ezfft1.o\
dzfftb.o\
dzfftf.o\
dzffti.o\
passb.o\
passb2.o\
passb3.o\
passb4.o\
passb5.o\
passf.o\
passf2.o\
passf3.o\
passf4.o\
passf5.o\
radb2.o\
radb3.o\
radb4.o\
radb5.o\
radbg.o\
radf2.o\
radf3.o\
radf4.o\
radf5.o\
radfg.o\
dfftb.o\
rfftb1.o\
dfftf.o\
rfftf1.o\
dffti.o\
rffti1.o\
dsinqb.o\
dsinqf.o\
dsinqi.o\
dsint.o\
sint1.o\
dsinti.o

lib$(LIB).a:	$(OBJ)
	ar -rcs lib$(LIB).a $(OBJ)

shared:$(OBJ)
	$(FC)  -shared -o lib$(LIB).so $(OBJ)

install:	lib$(LIB).a
	mv lib$(LIB).a /usr/local/lib
	rm *.o

installshared:lib$(LIB).so
	mv lib$(LIB).so /usr/local/lib
	rm *.o
	ldconfig 

test: 
	$(FC) test/test.f -L./ -l$(LIB)
	time ./a.out

clean: 
	 rm -f -r *.o *.a *.so
