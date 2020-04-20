CC:=gcc
HTSDIR=htslib
RCFLAGS:=-I/usr/share/R/include
CFLAGS= -fPIC -shared -I$(HTSDIR) -Wall ${RCFLAGS}

test: librbcf.so
	R --vanilla < test.R

librbcf.so : rbcf.c $(HTSDIR)/libhts.so
	$(CC) $(CFLAGS) -Wl,-soname,$(basename $@) -o $@ $< -L. -L $(HTSLIB) -lhts -lm -lz -lpthread

  
$(HTSDIR)/libhts.so:
	cd htslib && make libhts.so
	
clean:
	rm -f librbcf.so
	cd htslib && make clean
