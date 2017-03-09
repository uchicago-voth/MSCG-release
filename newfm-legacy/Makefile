HOME = /Users/lylu/newfm/newfm
MKLPATH = /usr/lib 
MKLINC = /usr/include
LSQRPATH = /Users/lylu/newfm
GSLPATH = /Users/lylu/gsl/lib
GSLINC = /Users/lylu/gsl/include
GMXPATH = /Users/lylu/xdrlib/lib
GMXINC = /Users/lylu/xdrlib/include/xdrfile
OPT = -O3 -std=c99
HEAD = common.h compute.h  control.h  matrix.h  output.h  top.h  trj.h  types.h

LIBS         = $(LSQRPATH)/pda_lsqr.o $(GSLPATH)/libgsl.a -framework Accelerate -lm -lxdrfile
LDFLAGS      = -L$(GMXPATH) -L$(MKLPATH) -L$(GSLPATH) 
CFLAGS	     = $(OPT) -I$(GSLINC) -I$(GMXINC) -I$(MKLINC)
CC           = gcc
LD           = gfortran-mp-4.4

newfm.o:	newfm.c  $(HEAD)
newfm:		newfm.o
		$(LD) $(LDFLAGS) -o $@ newfm.o $(LIBS)
range.o:	range.c range.h $(HEAD)
range:		range.o
		$(LD) $(LDFLAGS) -o $@ range.o $(LIBS)
avfm.o:		avfm.c avfm.h $(HEAD)
avfm:		avfm.o
		$(LD) $(LDFLAGS) -o $@ avfm.o $(LIBS)
disdis.o:	disdis.c  disdis.h $(HEAD)
disdis:		disdis.o
		$(LD) $(LDFLAGS) -o $@ disdis.o $(LIBS)
cg_gro.o:	cg_gro.c  cg_gro.h $(HEAD)
cg_gro:		cg_gro.o
		$(LD) $(LDFLAGS) -o $@ cg_gro.o $(LIBS)
force_comp.o:	force_comp.c  force_comp.h $(HEAD)
force_comp:	force_comp.o
		$(LD) $(LDFLAGS) -o $@ force_comp.o $(LIBS)
dis_3.o:	dis_3.c dis_3.h $(HEAD)
dis_3:		dis_3.o		
		$(LD) $(LDFLAGS) -o $@ dis_3.o $(LIBS)
ba_cor.o:	ba_cor.c ba_cor.h $(HEAD)
ba_cor:		ba_cor.o	
		$(LD) $(LDFLAGS) -o $@ ba_cor.o $(LIBS)
