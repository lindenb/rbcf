NAME=rbcf
.PHONY: all build install uninstall test doc

all:
	$(MAKE) uninstall build install test

build:
	R CMD build .

install: build
	R CMD INSTALL $(NAME)*.tar.gz

install2:
	(cd .. && R CMD INSTALL $(NAME))

uninstall:
	R CMD REMOVE $(NAME) || true

local:
	$(MAKE) uninstall || true
	$(MAKE) install2
	$(MAKE) test

test:
	cd tests && R --no-save < test.rbcf.R

doc:
	cd doc && bash generate.sh

clean:
	rm -f $(NAME)*.tar.gz
	$(MAKE) -C src clean
