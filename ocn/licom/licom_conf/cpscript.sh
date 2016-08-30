#! /bin/csh -f
setenv CASENAME "testlicom3"
echo $CASENAME
cp ccsm_comp_mod.F90  ../{$CASENAME}/SourceMods/src.drv
cp Macros.cloud_intel ../{$CASENAME}/Macros.cloud_intel
cp pop2.buildexe.csh.bak ../{$CASENAME}/Buildconf/pop2.buildexe.csh
