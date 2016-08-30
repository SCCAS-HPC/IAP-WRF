#
.SUFFIXES: .f90 .f .o

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f
SFC             =       ifort 
DM_FC           =       mpif90 -f90=$(SFC) -mcmodel=medium 
FC              =       $(DM_FC)
#FCDEBUG         =       -g -traceback -check all -check noarg_temp_created 
FCDEBUG         =       -g -traceback -check all -check noarg_temp_created -fstack-protector -assume protect_parens
FORMAT_FREE     =       -FR
BYTESWAPIO      =      -convert big_endian 
FCBASEOPTS      =     -w -ftz -fp-model precise -IPF_fp_relaxed -ip -fno-alias -fno-fnalias $(BYTESWAPIO) $(FCDEBUG) 
AR              =      ar
ARFLAGS         =      ru

LIBTARGET = libgeatmlib.a

OBJS    = 

MODULES =	geatm_comp_mct.o \
	geatm_vartype.o \
	parall.o

$(LIBTARGET) :  $(MODULES)
	$(AR) $(ARFLAGS) $(LIBTARGET) $(MODULES) $(OBJS)

.f90.o :
	$(FC) $(FCBASEOPTS) -I/opt/netcdf/include -I$(OBJROOT)/mct/mct -I$(OBJROOT)/csm_share \
	-I$(LIBROOT)/include -I./ -L$(LIBROOT) -lmct -lnetcdf -lpio -lcsm_share \
	-L$(OBJROOT)/mct/mpeu -lmpeu  -c -module ./ $*.f90	
.f.o :
	$(FC)  $(FCBASEOPTS) -I ./ -c -module ./ $*.f

clean:
	rm -f *.o *.mod

.PHONY: geatm clean

geatm_comp_mct.o: ./geatm_vartype.o \
	./parall.o
geatm_vartype.o :
parall.o :
