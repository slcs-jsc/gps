# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = events map perturbation prof response variance

# Library directories...
LIBDIR = -L ../lib/build/lib -L ../lib/build/lib64

# Include directories...
INCDIR = -I ../lib/build/include

# Linking...
STATIC = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED $(DEFS) -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lnetcdf -lgsl -lgslcblas -lm

# Profiling...
ifdef PROF
  CFLAGS += -O2 -g -pg 
else
  CFLAGS += -O3
endif

# Linking...
ifdef STATIC
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
