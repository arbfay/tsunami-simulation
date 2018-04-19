#
#  Makefile for Linux
#
CC       = gcc
LD       = gcc
CFLAGS   = -g -O3 -Dgraphic -Wall
LFLAGS   = -Wall -g -O3
LIBS     = -lglfw -lm -lGL -lGLU

PROG     = launch
LISTEOBJ = \
  main.o   tsunami.o

all:
	$(CC) -c main.c -o main.o $(CFLAGS)
	$(CC) -c tsunami.c -o tsunami.o $(CFLAGS)
	$(LD) $(LISTEOBJ) -o $(PROG) $(LFLAGS) $(LIBS)

main:
	$(CC) -c main.c -o main.o $(CFLAGS)
	$(LD) $(LISTEOBJ) -o $(PROG) $(LFLAGS) $(LIBS)

clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out
