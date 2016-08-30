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

MODULES =	adv_mark.o \
	advection.o \
	aerosol.o \
	alld.o \
	aqueous.o \
	cbmz.o \
	chemprod.o \
	convect43c.o \
	diffusion.o \
	diffusion_mark.o \
	drydep.o \
	eddyz.o \
	gaschem.o \
	geatm_comp_mct.o \
	geatm_vartype.o \
	getboundary.o \
	getclouddepth.o \
	getdu.o \
	getsmarkchange.o \
	hg.o \
	inout.o \
	inout2.o \
	inout3.o \
	model.o \
	modis.o \
	n_balance.o \
	parall.o \
	scavrat.o \
	setboundary.o \
	setboundtop1.o \
	setwindbound.o \
	termbal.o \
	tropause.o \
	wetdep_aer.o \
	wetdep_gas.o

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

geatm_comp_mct.o: ./adv_mark.o \
	./advection.o \
	./aerosol.o \
	./alld.o \
	./aqueous.o \
	./cbmz.o \
	./chemprod.o \
	./convect43c.o \
	./diffusion.o \
	./diffusion_mark.o \
	./drydep.o \
	./eddyz.o \
	./gaschem.o \
	./geatm_vartype.o \
	./getboundary.o \
	./getclouddepth.o \
	./getdu.o \
	./getsmarkchange.o \
	./hg.o \
	./inout.o \
	./inout2.o \
	./inout3.o \
	./model.o \
	./modis.o \
	./n_balance.o \
	./parall.o \
	./scavrat.o \
	./setboundary.o \
	./setboundtop1.o \
	./setwindbound.o \
	./termbal.o \
	./tropause.o \
	./wetdep_aer.o \
	./wetdep_gas.o
adv_mark.o : ./getsmarkchange.o
diffusion_mark.o : ./getsmarkchange.o
wetdep_aer.o : ./scavrat.o
wetdep_gas.o : ./scavrat.o \
	./drydep.o
advection.o :
aerosol.o :
alld.o :
aqueous.o :
cbmz.o :
chemprod.o :
convect43c.o :
diffusion.o : 
drydep.o :
eddyz.o :
gaschem.o :
geatm_vartype.o :
getboundary.o :
getclouddepth.o :
getdu.o :
getsmarkchange.o :
hg.o :
inout.o :
inout2.o :
inout3.o :
model.o :
modis.o :
n_balance.o :
parall.o :
scavrat.o :
setboundary.o :
setboundtop1.o :
setwindbound.o :
termbal.o :
tropause.o :
