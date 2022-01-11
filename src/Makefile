# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = events map perturbation prof response variance

# Installation directory...
DESTDIR ?= ../bin

# Library directories...
LIBDIR += -L ../lib/build/lib -L ../lib/build/lib64

# Include directories...
INCDIR += -I ../lib/build/include

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Compiler...
CC ?= gcc

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp $(OPT) -g -fopenmp

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lnetcdf -lgsl -lgslcblas -lm

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c libgps.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libgps.o jurassic.o $(LDFLAGS)

libgps.o: libgps.c libgps.h Makefile
	$(CC) $(CFLAGS) -c -o libgps.o libgps.c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

clean:
	rm -f $(EXC) *.o *~

dist:
	zip gps_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC)
