NAME=rbcf
.PHONY: all build install uninstall test

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

test:
	cd tests && R --no-save < test.rbcf.R

clean:
	rm -f $(NAME)*.tar.gz
	$(MAKE) -C src clean
