#REVISION SOFTWARE GIT:
REV = $(shell git rev-parse HEAD)
VER='character(len=41),parameter :: revision = "$(REV)"' > revision.inc
version:
	@echo $(VER)
