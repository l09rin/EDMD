# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdCylinder.c  mdCylinder.h
	$(CC) $(CFLAGS) mdCylinder.c -o md  $(LIBS)
