# 				M A K E  F I L E

#   main programs:	 class.c  cws.c  poly.c  nef.c  mori.c

SOURCES= Coord.c Rat.c Vertex.c Polynf.c
OBJECTS= $(SOURCES:.c=.o)

CLASS_SRC= Subpoly.c Subadd.c Subdb.c
CLASS_OBJ= $(CLASS_SRC:.c=.o)

NEF_SRC= E_Poly.c Nefpart.c LG.c
NEF_OBJ= $(NEF_SRC:.c=.o)

MORI_SRC= MoriCone.c SingularInput.c
MORI_OBJ= $(MORI_SRC:.c=.o)


CC=gcc

CFLAGS=-O3 -g -W -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
# CFLAGS=-O3 -g				      # add -g for GNU debugger gdb
# CFLAGS=-Ofast -O3 -mips4 -n32		      # SGI / 32 bit
# CFLAGS=-Ofast -O3 -mips4 -64                # SGI / 64 bit
# CFLAGS=-O3 -fast -non_shared -arch pca56    # statically linked for alpha_PC

#   targets : dependencies ; command
#             command
#             ...

all:	poly class cws nef mori cws0

clean:	;	rm -f *.o

cleanall: ;	rm -f *.o *.x core


poly:	poly.o $(OBJECTS) LG.o Global.h LG.h
	$(CC) $(CFLAGS) -o poly.x poly.o $(OBJECTS) LG.o

class:	class.o $(OBJECTS) $(CLASS_OBJ) Global.h Subpoly.h
	$(CC) $(CFLAGS) -o class.x class.o $(OBJECTS) $(CLASS_OBJ)

cws:    cws.o $(OBJECTS) LG.o Global.h LG.h
	$(CC) $(CFLAGS) -o cws.x cws.o $(OBJECTS) LG.o

cws0:   cws0.o $(OBJECTS) LG.o Global.h LG.h
	g++ $(CFLAGS) -o cws0.x weight_system_store_set.cpp cws0.o $(OBJECTS)

nef:    nef.o $(OBJECTS) $(NEF_OBJ) Global.h 
	$(CC) $(CFLAGS) -o nef.x nef.o $(OBJECTS) $(NEF_OBJ)

mori:   mori.o  $(OBJECTS) $(MORI_OBJ) LG.o Mori.h 
	$(CC) $(CFLAGS) -o mori.x mori.o $(OBJECTS) $(MORI_OBJ) LG.o



#			     D E P E N D E N C I E S

Coord.o:        Global.h Rat.h 
Polynf.o:		Global.h Rat.h
Rat.o:          Global.h Rat.h 
Subpoly.o:      Global.h Rat.h Subpoly.h  
Subadd.o:		Global.h Subpoly.h
Vertex.o:       Global.h Rat.h  
Subdb.o:		Global.h Subpoly.h
LG.o:           Global.h Rat.h LG.h

E_Poly.o:       Global.h Nef.h Rat.h
Nefpart.o:		Global.h Nef.h

MoriCone.o:      Global.h Rat.h Mori.h
SingularInput.o: Global.h Mori.h

poly.o:         Global.h LG.h
class.o:		Global.h Subpoly.h
cws.o:			Global.h LG.h Rat.h
cws0.o:		Global.h LG.h Rat.h
nef.o:          Global.h Nef.h LG.h
mori.o:     	Global.h LG.h Mori.h



#	experimental stuff ...
#
#bpoly: bpoly.o 	Global.h $(OBJECTS)
#	$(CC)  	$(CFLAGS) -o bpoly.x  bpoly.o $(OBJECTS)

#bpoly.o: bpoly.c Global.h

#vnl:	vnl.o  $(OBJECTS) Global.h
#	$(CC)   $(CFLAGS) -o  vnl.x  vnl.o  $(OBJECTS)

#vnl.o:  vnl.c	Global.h

#gen:    gen.o   $(OBJECTS) $(NEF_OBJ) Global.h
#	$(CC)   $(CFLAGS) -o  gen.x  gen.o  $(OBJECTS) $(NEF_OBJ)
#gen.o: 	gen.c
