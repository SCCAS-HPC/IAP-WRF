#!/bin/csh -f

#-------------------------------------------------------
# [1] Set necessary environment variables
#-------------------------------------------------------

#setenv FC xlf
#setenv FC pgi
setenv FC ifort

setenv EDGE_N    90.
setenv EDGE_E   180.
setenv EDGE_S   -90.
setenv EDGE_W  -180.

setenv NLON  360
setenv NLAT  180

setenv TASKS 1

setenv ROOTDIR $HOME/CoLM_RTM

# 1. set clm include directory root
setenv CLM_INCDIR $ROOTDIR/include

# 2. set clm raw land data directory root
setenv CLM_RAWDIR $ROOTDIR/rawdata

# 3. set clm surface data rectory root
setenv CLM_SRFDIR $ROOTDIR/mksrfdata

# 4. set clm input data directory root
setenv CLM_DATADIR $ROOTDIR/data

# 5. set clm initial directory root
setenv CLM_INIDIR $ROOTDIR/mkinidata

# 6. set clm source directory root
setenv CLM_SRCDIR $ROOTDIR/main

# 7. set executable directory
setenv CLM_EXEDIR $ROOTDIR/run

# 8. set output directory
setenv CLM_OUTDIR /p2/jidy/colm_rtm

mkdir -p $CLM_OUTDIR

#------------------------------------------------------
# [2] Build define.h in ./include directory
#------------------------------------------------------

\cat >! $CLM_INCDIR/define.h << EOF
#define LINUX
#undef  RDGRID
#undef  SOILINI
#undef  BATS
#undef  SIB2
#undef  IGBP
#define USGS
#define EcoDynamics
#define LANDONLY
#undef  LAND_SEA
#undef  SINGLE_POINT
#undef  GLACIER
#undef  MAPMASK
#undef  METMASK
#define PRINCETON
#undef  GSWP2
#define WR_DAILY
#define RTM
#define OFFLINE
#undef  COUP_CSM
EOF

if ($TASKS > 1) then
   \cat >> $CLM_INCDIR/define.h << EOF
#define SPMD
EOF
endif

#------------------------------------------------------------------------
# [3] Compling clm surface data making, clm initialization, clm timeloop
#------------------------------------------------------------------------

echo 'Compiling mksrfdata...'
cd $CLM_SRFDIR

#make -f Makefile.${FC} clean
make -f Makefile.${FC} >>& $CLM_EXEDIR/compile.log.clm || exit 5

cp -f $CLM_SRFDIR/srf.x $CLM_EXEDIR/srf.x

echo 'Compiling mkinidata...'
cd $CLM_INIDIR

#make -f Makefile.${FC} clean
make -f Makefile.${FC} >>& $CLM_EXEDIR/compile.log.clm || exit 5

cp -f $CLM_INIDIR/initial.x $CLM_EXEDIR/initial.x

echo 'Compiling main...'
cd $CLM_SRCDIR

#make -f Makefile.${FC} clean
make -f Makefile.${FC} >>& $CLM_EXEDIR/compile.log.clm || exit 5

cp -f $CLM_SRCDIR/clm.x $CLM_EXEDIR/clm.x

#------------------------------------------------------------------------
# [4] Executing clm surface data making, clm initialization, clm timeloop
#------------------------------------------------------------------------

cd $CLM_EXEDIR

# Create an input parameter namelist file for srf.x

\cat >! $CLM_EXEDIR/srfdat.stdin << EOF
&mksrfexp
fmetmask    = ''
fmapmask    = ''
fglacier    = ''
fgridname   = ''
fdemname    = '$CLM_RAWDIR/dem-usgs.30s'
fmaskname   = '$CLM_RAWDIR/lwmask-usgs.30s'
flandname   = '$CLM_RAWDIR/veg-usgs.30s'
fsolaname   = '$CLM_RAWDIR/soilcat.30s'
fsolbname   = '$CLM_RAWDIR/soilcatb.30s'
fsurdat     = '$CLM_DATADIR/srfdata.1deg'
lon_points  =  $NLON
lat_points  =  $NLAT
edgen       =  $EDGE_N
edgee       =  $EDGE_E
edges       =  $EDGE_S
edgew       =  $EDGE_W
/
EOF

echo 'Executing CLM Making Surface Data'

if($TASKS > 1)then
   mpirun -np $TASKS $CLM_EXEDIR/srf.x < $CLM_EXEDIR/srfdat.stdin >& $CLM_EXEDIR/clm.log.srf || exit 5
else
   $CLM_EXEDIR/srf.x < $CLM_EXEDIR/srfdat.stdin >& $CLM_EXEDIR/clm.log.srf || exit 5
endif

echo 'CLM Making Surface Data Completed'

# Create an input parameter namelist file for initial.x

\cat >! $CLM_EXEDIR/inidat.stdin << EOF
&clminiexp
site           = 'global'
greenwich      = .true.
start_yr       =  2006 
start_jday     =  1
start_sec      =  1800
fsrf           = '$CLM_DATADIR/srfdata.1deg'
fgrid          = '$CLM_DATADIR/gridata.1deg'
frivinp_rtm    = '$CLM_DATADIR/rdirc.05'
fsbcini        = '$CLM_DATADIR/sbcini.1deg'
flai           = ''
fsoilini       = ''
fmet           = '/p1/jidy/Data_save/princeton_30min'
fconst         = '$CLM_OUTDIR/global-rstTimeConst'
frestart       = '$CLM_OUTDIR/global-rstTimeVar'
fout           = '$CLM_OUTDIR/global'
finfolist      = '$CLM_EXEDIR/clmini.infolist'
lon_points     =  $NLON
lat_points     =  $NLAT
dtime          =  1800
mstep          =  1500
/
EOF

echo 'Executing CLM Initialization'

$CLM_EXEDIR/initial.x <$CLM_EXEDIR/inidat.stdin >& $CLM_EXEDIR/clm.log.initial || exit 5

echo 'CLM Initialization Completed'

# Create an input parameter namelist file for clm.x

mv -f $CLM_EXEDIR/clmini.infolist $CLM_EXEDIR/timeloop.stdin

# Create flux export namelist file for clm.x
# Don't change the sequence of the FLUX array elements !*!

set FLUX = ( +taux     +tauy     +fsena    +lfevpa    +fevpa    +fsenl    \
             +fevpl    +etr      +fseng    +fevpg     +fgrnd    +sabvsun  \
             +sabvsha  +sabg     +olrg     +rnet      +xerr     +zerr     \
             +rsur     +rnof     +assim    +respc     +tss_01   +tss_02   \
             +tss_03   +tss_04   +tss_05   +tss_06    +tss_07   +tss_08   \
             +tss_09   +tss_10   +wliq_01  +wliq_02   +wliq_03  +wliq_04  \
             +wliq_05  +wliq_06  +wliq_07  +wliq_08   +wliq_09  +wliq_10  \
             +wice_01  +wice_02  +wice_03  +wice_04   +wice_05  +wice_06  \
             +wice_07  +wice_08  +wice_09  +wice_10   +tg       +tlsun    \
             +tlsha    +ldew     +scv      +snowdp    +fsno     +sigf     \
             +green    +lai      +sai      +alb_11    +alb_12   +alb_21   \
             +alb_22   +emis     +z0ma     +trad      +ustar    +tstar    \
             +qstar    +zol      +rib      +fm        +fh       +fq       \
             +tref     +qref     +u10m     +v10m      +f10m     +us       \
             +vs       +tm       +qm       +prc       +prl      +pbot     \
             +frl      +solar    )

@ i = 0

set flux_exp = "flux_exp="

foreach str ($FLUX)
   @ i = $i + 1
   if("$str" =~ +*) then
      set flux_exp = "$flux_exp +$i"
   else
      set flux_exp = "$flux_exp -$i"
   endif
end

\cat >! $CLM_EXEDIR/flux.stdin << EOF
&flux_nml
$flux_exp
/
EOF

echo 'Executing CLM Time-looping'

if($TASKS > 1)then
   mpirun -np $TASKS $CLM_EXEDIR/clm.x < $CLM_EXEDIR/timeloop.stdin >& $CLM_EXEDIR/clm.log.timeloop || exit 5
else
   $CLM_EXEDIR/clm.x < $CLM_EXEDIR/timeloop.stdin >& $CLM_EXEDIR/clm.log.timeloop || exit 5
endif

echo 'CLM Running Completed'
