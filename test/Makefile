TSTDIR=src
TSTBIN=bin

all: bin
	make -C $(TSTDIR)

bin: 
	mkdir $(TSTBIN)

test:
	./test.sh

clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
