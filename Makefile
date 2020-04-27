NAME=rbcf
.PHONY: all build install uninstall test doc

all:
	$(MAKE) uninstall build install test

build:
	rm -f $(NAME)*.tar.gz
	R CMD build .

install: build
	R CMD INSTALL $(NAME)*.tar.gz

install2:
	(cd .. && R CMD INSTALL $(NAME))

check:
	R CMD check --as-cran .

uninstall:
	R CMD REMOVE $(NAME) || true

local:
	$(MAKE) uninstall || true
	$(MAKE) install2
	$(MAKE) test

test: doc

doc:
	cd tests && bash generate.sh

clean:
	rm -f $(NAME)*.tar.gz
	$(MAKE) -C src clean
