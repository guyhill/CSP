CC = g++	# Some C++ is used, should remove that...
CFLAGS = -g -Wall -DVSCIP -Dsoplex -DMICROSOFT2 -fPIC
LDFLAGS = $(CFLAGS) -shared -Wl,-Bsymbolic -Wl,-zdefs
LIBS = -lscip
TARGET = libcsp.so

SRCDIR = src
BUILDDIR = build

HEADERS = $(wildcard $(SRCDIR)/*.h)
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SOURCES))

.PHONY: all clean

all: $(BUILDDIR)/$(TARGET)

$(BUILDDIR):
	mkdir -p $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HEADERS) | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/$(TARGET): $(OBJECTS) | $(BUILDDIR)
	$(CC) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -rf $(BUILDDIR)
