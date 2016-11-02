all: make.inc lib tests

make.inc:
	@echo 'Error cannot find make.inc file'
	@exit 1

lib:
	$(MAKE) -C lib/

tests: lib
	$(MAKE) -C tests/

clean:
	$(MAKE) -C lib/ clean
	$(MAKE) -C tests/ clean
	rm -rf ./config.log ./config.status ./autom4te.cache

.PHONY: all lib tests clean