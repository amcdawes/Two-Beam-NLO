SRC	= defocus.cc fftw++.cc
OBJ	= $(SRC:.cc=.o)
LIBS	= -lfftw3 -lpng -lpngwriter -lfreetype -lz
EXE	= defocus

CC	= /usr/bin/g++
CFLAGS	= -funroll-loops -Wall -O4 -DNDEBUG -I/sw/include `freetype-config --cflags`
LIBPATH	= -L/sw/lib
LDFLAGS	= -o $(EXE) $(LIBPATH) $(LIBS)
CFDEBUG = -ansi -pedantic -Wall -g -DDEBUG $(LDFLAGS)
RM      = /bin/rm -f


# Compile and Assemble C Source Files into Object Files
%.o: %.cc
	$(CC) -c $(CFLAGS) $*.cc

# Link all Object Files with external Libraries into Binaries
$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS)

# Objects depend on these Libraries
$(OBJ): $(INCL)

# Create a gdb/dbx Capable Executable with DEBUG flags turned on
debug:
	$(CC) $(CFDEBUG) $(SRC)

# Clean Up Objects, Exectuables, Dumps out of source directory
clean:
	$(RM) $(OBJ) $(EXE) core a.out
