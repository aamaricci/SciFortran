TSTDIR=src
TSTBIN=bin

all: bin
	make -C $(TSTDIR)

bin:
	mkdir $(TSTBIN)
	mkdir $(TSTBIN)/SF_ARRAYS
	mkdir $(TSTBIN)/SF_COLORS
	mkdir $(TSTBIN)/SF_CONSTANTS
	mkdir $(TSTBIN)/SF_FONTS
	mkdir $(TSTBIN)/SF_MISC
	mkdir $(TSTBIN)/SF_TIMER
	mkdir $(TSTBIN)/SF_INTEGRATE
	mkdir $(TSTBIN)/SF_PARSE_INPUT
	mkdir $(TSTBIN)/SF_IOTOOLS
	mkdir $(TSTBIN)/SF_DERIVATE
	mkdir $(TSTBIN)/SF_FFT
	mkdir $(TSTBIN)/SF_INTERPOLATE
# mkdir $(TSTBIN)/SF_LINALG
# mkdir $(TSTBIN)/SF_MPI
	mkdir $(TSTBIN)/SF_OPTIMIZE
# mkdir $(TSTBIN)/SF_RANDOM
# mkdir $(TSTBIN)/SF_SPECIAL
# mkdir $(TSTBIN)/SF_SP_LINALG
	echo SF_ARRAYS > $(TSTBIN)/list_dir
	echo SF_COLORS > $(TSTBIN)/list_dir
	echo SF_CONTANTS > $(TSTBIN)/list_dir
	echo SF_FONTS > $(TSTBIN)/list_dir
	echo SF_MISC > $(TSTBIN)/list_dir
	echo SF_TIMER > $(TSTBIN)/list_dir
	echo SF_INTEGRATE >> $(TSTBIN)/list_dir
	echo SF_PARSE_INPUT >> $(TSTBIN)/list_dir
	echo SF_IOTOOLS >> $(TSTBIN)/list_dir
	echo SF_DERIVATE >> $(TSTBIN)/list_dir
	echo SF_FFT >> $(TSTBIN)/list_dir
	echo SF_INTERPOLATE >> $(TSTBIN)/list_dir
# echo SF_LINALG >> $(TSTBIN)/list_dir
# echo SF_MPI >> $(TSTBIN)/list_dir
	echo SF_OPTIMIZE >> $(TSTBIN)/list_dir
# echo SF_RANDOM >> $(TSTBIN)/list_dir
# echo SF_SPECIAL >> $(TSTBIN)/list_dir
# echo SF_SP_LINALG >> $(TSTBIN)/list_dir

test:
	./test.sh

clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
