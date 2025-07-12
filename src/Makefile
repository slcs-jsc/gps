# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = events map perturbation prof response variance

# Installation directory...
DESTDIR ?= ../bin

# Library directories...
LIBDIR += -L ../libs/build/lib -L ../libs/build/lib64

# Include directories...
INCDIR += -I ../libs/build/include

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
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)

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

.PHONY : all check clean coverage cppcheck dist doxygen indent install lizard mkdocs strip uninstall

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

cppcheck:
	cppcheck --enable=all ./

dist:
	cd ..  && tar cvfz gps_$(VERSION).tgz bin/* \
          docs/Doxyfile example/* libs/build.sh libs/*gz src/* tests/*

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
