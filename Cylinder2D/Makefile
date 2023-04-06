# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdCylinder2D.c  mdCylinder2D.h
	$(CC) $(CFLAGS) mdCylinder2D.c -o md  $(LIBS)
