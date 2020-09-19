#*****************************************************************

CC = gcc
RM = rm

ROOT=../../
LIB=${ROOT}/lib

CFLAGS=-I${LIB}/sml -g
LFLAGS=-L${LIB}/sml ${LIB}/sml/sml.lib

#*****************************************************************

all: transformer
.PHONY: all

transformer    : transformer.c
	$(CC) $(CFLAGS) -o $@  transformer.c  $(LFLAGS) -lm

.PHONY: clean
clean:
	$(RM) -f *.o transformer
	$(RM) -rf *.dSYM
