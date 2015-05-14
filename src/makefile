OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) \
	-I/opt/intel/ics/2013/028/impi/4.1.0.024/intel64/include/ -DMPICH_IGNORE_CXX_SEEK

# Debug flags
#CFLAGS += -g

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
#LDFLAGS += -g
LDFLAGS += -L/opt/intel/ics/2013/028/impi/4.1.0.024/intel64/lib

CXX = mpicc

gem: gem.C
	$(CXX) $(CFLAGS) gem.C
	$(CXX) -o gem gem.o $(LDFLAGS)
	rm gem.o

test: random_number_test.C
	$(CXX) $(CFLAGS) random_number_test.C
	$(CXX) -o random_number_test random_number_test.o $(LDFLAGS)
	rm random_number_test.o
