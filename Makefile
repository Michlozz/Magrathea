DEBUG=-g

CC=g++ #-fopenmp

CFLAGS= -std=c++0x -O3 -Wall  # -I/opt/apps/spack/linux-centos7-x86_64/gcc/10.3.0/gsl-2.7-t7jnjnv7/include
LDFLAGS= -L/opt/apps/spack/linux-centos7-x86_64/gcc/10.3.0/gsl-2.7-t7jnjnv7/lib -lgsl -lgslcblas -lm

SRCDIR = src
LIBDIR = src

BIN = planet
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(patsubst $(SRCDIR)/%,%,$(SOURCES))
OBJ := $(patsubst %.cpp,%.o,$(OBJ))
OBJ := $(addprefix ./$(LIBDIR)/,$(OBJ))

planet: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp  
	$(CC) -o $@ -c $< $(CFLAGS) $(DEBUG)

clean:
	rm -f *.o ./src/*.o planet
