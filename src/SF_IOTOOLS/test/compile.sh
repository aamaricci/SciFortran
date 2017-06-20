#DFLAG="-O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none"
cd ../
echo "gfortran $DFLAG -c IOFILE.f90 IOPLOT.f90 IOREAD.f90 && mv -v *.o *.mod test/"
gfortran $DFLAG -c IOFILE.f90 IOPLOT.f90 IOREAD.f90 && mv -v *.o *.mod test/
cd test/
echo "gfortran $DFLAG *.o testSPLOT_SREAD.f90"
gfortran $DFLAG *.o testSPLOT_SREAD.f90
