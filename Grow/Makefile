# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdgrow.c  mdgrow.h
	$(CC) $(CFLAGS) mdgrow.c -o grow  $(LIBS)

md2D: mdgrow2D.c  mdgrow2D.h
	$(CC) $(CFLAGS) mdgrow2D.c -o grow2D  $(LIBS)
