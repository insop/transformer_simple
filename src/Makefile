#*****************************************************************

CC = gcc
RM = rm

CFLAGS=-I../sml -g
LFLAGS=-L../sml

#*****************************************************************

all: transformer
.PHONY: all

transformer    : transformer.c
	$(CC) $(CFLAGS) -o $@  transformer.c  $(LFLAGS) ../sml/sml.lib  -lm

.PHONY: clean
clean:
	$(RM) -f *.o transformer
	$(RM) -rf *.dSYM
