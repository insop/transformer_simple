#*****************************************************************

CC = gcc
RM = rm

ROOT=../../
CFLAGS=-I${ROOT}/sml -g
LFLAGS=-L${ROOT}/sml ${ROOT}sml/sml.lib


#*****************************************************************

all: transformer
.PHONY: all

transformer    : transformer.c
	$(CC) $(CFLAGS) -o $@  transformer.c  $(LFLAGS) -lm

.PHONY: clean
clean:
	$(RM) -f *.o transformer
	$(RM) -rf *.dSYM
