# makefile: make the SML library
#*****************************************************************

# macros
SMLOBJ = smlmem.o smlacf.o smlio.o smlerr.o smlcomp.o smlchol.o smlqr.o smlsvd.o smleig.o marmem.o mario.o
SMLBASE = smlbase.h sml.cfg  smlacf.h smlacm.h smlmem.h smlerr.h

CC = gcc
CC = wcl386

CC = CC
SMLCFLAGS = -fast -g

#*****************************************************************

# The Small Matrix Library (SML)

sml.lib : $(SMLOBJ)
	ar -r -u sml.lib $(SMLOBJ)
	rm -f *.o

smlacf.o  : smlacf.c  $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlacf.c

smlmem.o : smlmem.c   $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlmem.c

smlerr.o : smlerr.c   $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlerr.c

smlio.o : smlio.c   smlio.h       $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlio.c

smlcomp.o : smlcomp.c smlcompv.c smlcomp.h  $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlcomp.c

smlchol.o : smlchol.c smlcomp.h   $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlchol.c

smlqr.o     : smlqr.c smlcomp.h   $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlqr.c

smlsvd.o    : smlsvd.c smlcomp.h  $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smlsvd.c

smleig.o    : smleig.c smlcomp.h  $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) smleig.c

marmem.o : marmem.c mar.h         $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) marmem.c

mario.o : mario.c   mar.h smlio.h $(SMLBASE)
	$(CC) -c $(SMLCFLAGS) mario.c



# programs that use SML

demo    : demo.c
	$(CC)  $(SMLCFLAGS)  -o $@  demo.c   sml.lib  -lm

.PHONY: clean
clean	:
	rm -f sml.lib
