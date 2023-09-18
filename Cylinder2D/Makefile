# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdCylinder2D.c  mdCylinder2D.h  mystring.c
	$(CC) $(CFLAGS) mdCylinder2D.c -o md  $(LIBS)

mdsqs: squareshoulder2D.c  squareshoulder2D.h  mystring.c
	$(CC) $(CFLAGS) squareshoulder2D.c -o mdsqs  $(LIBS)

mdhdd: harddisks.c  harddisks.h  mystring.c
	$(CC) $(CFLAGS) harddisks.c -o mdhdd  $(LIBS)
