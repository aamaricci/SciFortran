include ../../make.inc
include ARmake.inc

serial: 
	@make -C $(UTILdir) $(PRECISIONS)
	@make -C $(SRCdir) $(PRECISIONS)

parallel: 
	@make -C $(UTILdir) $(PRECISIONS)
	@make -C $(SRCdir) $(PRECISIONS)
	@make -C $(PUTILdir) $(PRECISIONS)
	@make -C $(PSRCdir) $(PRECISIONS)

