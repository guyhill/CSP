CC = g++
CFLAGS = -g -Wall -DLPSC -DVSCIP
LDFLAGS = $(CFLAGS) -Wl,-Bsymbolic -Wl,-zdefs
LIBS = -L../build -lcsp -lscip

SRCDIR = src
BUILDDIR = ../build

.PHONY: all clean

all: $(BUILDDIR)/solve_csp

$(BUILDDIR):
	mkdir -p $@

$(BUILDDIR)/solve_csp: $(SRCDIR)/solve_csp.cpp | $(BUILDDIR)
	$(CC) $(SRCDIR)/solve_csp.cpp $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -rf $(BUILDDIR)
