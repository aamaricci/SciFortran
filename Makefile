.PHONY: serial parallel update

serial:
	@make -C ./src serial
parallel:
	@make -C ./src parallel
update:
	@make -C ./src update

portable:
	@make portable -C ./src/

clean: cleanlib
	@make clean -C ./src/

cleanlib:
	@rm -fv lib/*.a
