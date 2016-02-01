# C compiler
CC	= gcc

# Compiler options
CFLAGS	= -lproj -lnetcdf

all: .o
	$(CC) -o fixproj fixproj.o $(CFLAGS)
	$(CC) -o genseed genseed.o $(CFLAGS)

.o:
	$(CC) $(CFLAGS) -c fixproj.c
	$(CC) $(CFLAGS) -c genseed.c

clean:
	rm *.o

