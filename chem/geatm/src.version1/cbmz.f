! This computer software was prepared by Battelle Memorial Institute, 
c hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
c with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOT THE
c CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
c LIABILITY FOR THE USE OF THIS SOFTWARE.
c---------------------------------------------------------------------------------------------
c Box-Model     : CBMZ.1.2
c                 Regime-dependent CBM-Z + Condensed DMS Chemistry
c
c Integrator    : LSODES/LSODE
c
c Author        : Rahul A. Zaveri, PhD
c                 Pacific Northwest National Laboratory
c                 Atmospheric Sciences Technical Group
c                 P.O. Box 999, MSIN K9-30
c                 Richland, WA 99352
c                 Phone: (509) 372-6159  Fax: (509) 372-6168
c                 Email: Rahul.Zaveri@pnl.gov
c                 Website: www.pnl.gov/atmos_sciences/raz
c
c Update        : April 2006
c
c Support       : (1) U.S. Environmental Protection Agency (EPA) Aerosol Program
c                 (2) NASA Aerosol Program 
c                 (3) U.S. Department of Energy (DOE) Atmospheric Chemistry Program
c
c Terms of Use  : (1) CBM-Z may not be included in any commercial package, or used
c                     for any commercial applications without prior authorization
c                     from the author.
c                 (2) The code may be used for educational or non-profit purposes
c                     only. Any other usage must be first approved by the author.
c                 (3) The code may not be further distributed without the author's 
c                     prior consent.
c                 (4) No portion of the CBM-Z source code can be used in other codes
c                     without the author's prior permission.
c                 (5) The CBM-Z code is provided on an as-is basis, and the author
c                     bears no liability from its usage.
c                 (6) Publications resulting from the usage of CBM-Z should cite
c                     the reference below for proper acknowledgment.
c
c References    : Zaveri R.A. and L.K. Peters (1999) A new lumped structure photochemical
c                 mechanism for large-scale applications, J. Geophys. Res, 104, 30,387-30,415.
c
c---------------------------------------------------------------------------------------------
c
      subroutine  cbmz(cppb,factcld)
      implicit real(a-h,o-z), integer(i-n) 
      include 'chm.inc'
      include 'gas.inc'

      dimension cppb(ngas_max)


      
      call SetRunParameters

      call SetAirComposition
      call LoadPeroxyParameters			! Aperox and Bperox
      call DoMassBalance			! initial elemental mass balance

      call Print(cppb)
c------------------------------------------------------------------
c
c main time-loop begins...
      do 100 it = 1, ncstep

      call UpdateTime

      if(msolar.eq.1)then
       call SolarZenithAngle
      endif

      call UpdateMetFields

      call UpdateEmissions

      call IntegrateChemistry(factcld)

      call DoMassBalance

      call Print(cppb)
      
100   continue

c------------------------------------------------------------------


c
c      stop
      return !lijie
      end            


c---------------------------------------------------------------------


      subroutine ReadInputFile(lin)
      include 'chm.inc'
      include 'gas.inc'

      character*40 dword

      read(lin,*)dword

c begin time from 12:00 (noon) March 21 [min]
      read(lin,*)dword, tbeg_dd, tbeg_hh, tbeg_mm, tbeg_ss, dword
      read(lin,*)dword, trun_dd, trun_hh, trun_mm, trun_ss, dword
      read(lin,*)dword, dt_min,dword ! transport time-step [min]
      read(lin,*)dword, rlon,  dword ! longitude [deg]
      read(lin,*)dword, rlat,  dword ! latitude [deg]
      read(lin,*)dword, zalt_m,dword ! altitude  above mean sea level [m]
      read(lin,*)dword, RH,    dword ! relative humidity [%]
      read(lin,*)dword, te,    dword ! temperature [K]
      read(lin,*)dword, pr_atm,dword ! pressure [atm]
      read(lin,*)dword, msolar,dword ! msolar flag
      read(lin,*)dword, mphoto,dword ! mphoto flag
      read(lin,*)dword, iprint,dword ! freq of output
c
c
c read index, species name, initial conc, and emissions
      read(lin,*)dword
      read(lin,*)dword

      read(lin,*)dword ! GAS
      do i=1, ngas_max
      read(lin,*)k, species(k), cnn(k), emission(k)
      enddo

c      write(6,*)'Finished reading all inputs...'
      return
      end





c***********************************************************************
      subroutine IntegrateChemistry(factcld)
      include 'chm.inc'
      include 'gas.inc'

      t_in = told_sec
      t_out= tcur_sec

      call GasChemistry(t_in, t_out,factcld)

      return
      end




c***********************************************************************
      subroutine UpdateEmissions
      include 'chm.inc'

c update emission fluxes here

      return
      end




c***********************************************************************
      subroutine UpdateMetFields
      include 'chm.inc'

c update temperature, pressure, cloud flag, etc. here

      return
      end





c*****************************************************************
c     subroutine for printing output at iprint time steps
c 
      subroutine print(cppb)
      implicit real(a-h,o-z), integer(i-n)
      include 'chm.inc'
      include 'gas.inc'
      dimension cppb(ngas_max)

c
c=================================================================
c
c     Converting (molecules/cc) to (ppb)
c-----------------------------------------------------------------
      do l=1,ngas_max
        cppb(l) = (cnn(l)/cair_mlc)*ppb
      enddo

c=================================================================
c
c     Gas-Phase Species
c-----------------------------------------------------------------
c                                                            comm by lijie
c      if(it.eq.0)write(20,770)

c      write(20,7)tcur_hrs,cppb(ko3),cppb(kno),cppb(kno2),
c     &                    cppb(khno3),cppb(khono),cppb(kpan),
c     &                    cppb(konit),cppb(kh2o2),cppb(kolet),
c     &                    cppb(kolei),cppb(kpar),cppb(kald2),
c     &                    cppb(kisop)

c
c===================================================================
c
c     FORMATS
c-------------------------------------------------------------------

!7     format(f6.1, 13(2x,f7.3)) ! original
7       format(f6.1, 13(2x,f10.3)) !lijie modify
770   format('  Time     O3       NO       NO2      HNO3'
     &       '     HONO     PAN     ONIT      H2O2    OLET'
     &       '     OLEI     PAR      ALD2     ISOP')


      return
      end






c************************************************************************
c subroutine SelectGasRegime: selects an optimum combination of gas-phase
c                             mechanisms 
c
c input : cnn       = full species concentrations array (molec/cc)
c
c output: iregime   = 1     : com
c                   = 2     : com + urb
c                   = 3     : com + urb + bio
c                   = 4     : com + mar
c                   = 5     : com + urb + mar
c                   = 6     : com + urb + bio + mar    
c         ngas      = number of gas-phase species in the selected mechanism
c
c author: Rahul A. Zaveri
c date  : february 1996
c
c---------------------------------------------------------------------
      subroutine SelectGasRegime(ntot)
      include 'chm.inc'
      include 'gas.inc'

      cutoff_conc = 5.e+6     ! [molec/cc]

c
c---------------------------------------------
c initialize regime flags to zero...
      m_com = 1 ! 1 (always)
      m_urb = 0	! 0 or 1
      m_bio = 0	! 0 or 2
      m_mar = 0	! 0 or 3

c
c decide mechanism flags...
      do k = kpar, kxpar
      if( (cnn(k) .gt. cutoff_conc)  .or.
     &    (emission(k) .gt. 0.0)    ) m_urb=1
      enddo
c
      do k = kisop, kisopo2
      if( (cnn(k) .gt. cutoff_conc)  .or.
     &    (emission(k) .gt. 0.0)    ) m_bio=2
      enddo
c
      do k = kdms, ksulfhox
      if( (cnn(k) .gt. cutoff_conc)  .or.
     &    (emission(k) .gt. 0.0)     ) m_mar=3
      enddo
c
      iregime = m_com + m_urb*((2-m_bio)/2) + m_bio + m_mar

c-------------------------------------------------------------------------
      goto (1,2,3,4,5,6), iregime

1     ntot = nreg1
      return

2     ntot = nreg2
      return

3     ntot = nreg3
      return

4     ntot = nreg4
      return

5     ntot = nreg5
      return

6     ntot = nreg6
      return

      end







c***********************************************************************
      subroutine SetAirComposition
      include 'chm.inc'
      include 'gas.inc'

c-----------------------------------------------------------------------
c set bulk air composition in (molec/cc)
      cair_mlc = avogad*pr_atm/(82.056*te)	! air conc [molec/cc]
c      write(6,*)' air concentraton = ', cair_mlc, ' molec/cc'
      o2       = 0.21*cair_mlc
      h2       = 0.58*1.e-6*cair_mlc
      h2o      = WaterVapor(RH, cair_mlc, te, pr_atm)

c      write(6,*)'h2o = ', h2o

c-----------------------------------------------------------------------
c conversion factor for converting [ppb] to [molec/cc]
c
      ppb = 1.e+9
c
c-------------------------------------------------------------
c converting gas-phase conc from [ppb] to [molec/cc]
      do l=1, ngas_max
        cnn(l) = cnn(l)*cair_mlc/ppb
      enddo
c
c
c convert from ppb/hr to molec/cc/s
      do l=1, ngas_max
        emission(l) = emission(l)*cair_mlc/ppb/3600.
      enddo
c
      return
      end            






c***********************************************************************
c subroutine MapGas_Bio: maps cnn to and fro stot for the biogenic 
c                        gas-phase mechanism.
c
c nomenclature:
c cnn       = full species concentration array.
c stot      = subset of cnn. species concentration array to be supplied to
c             lsodes. length of stot depends on the selected mechanism
c iregime   = selected chemical regime (1-6)
c imap      = 0 : map cnn to stot
c           = 1 : map stot to cnn
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine MapGas_Bio(stot,imap)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)

      emit(iisop)	=emission(kisop)
      emit(iisoprd)	=emission(kisoprd)
      emit(iisopp)	=emission(kisopp)
      emit(iisopn)	=emission(kisopn)
      emit(iisopo2)	=emission(kisopo2)

      if(imap.eq.0)then    ! map cnn into stot
      stot(iisop)	=cnn(kisop)
      stot(iisoprd)	=cnn(kisoprd)
      stot(iisopp)	=cnn(kisopp)
      stot(iisopn)	=cnn(kisopn)
      stot(iisopo2)	=cnn(kisopo2)
      stot(iterp)       =cnn(kterp)
      stot(isv3)        =cnn(ksv3)
      stot(isv4)        =cnn(ksv4)
      stot(isv5)        =cnn(ksv5)
      stot(isv6)        =cnn(ksv6)
c
      else                 ! map stot back into cnn
      cnn(kisop)	=stot(iisop)
      cnn(kisoprd)	=stot(iisoprd)
      cnn(kisopp)	=stot(iisopp)
      cnn(kisopn)	=stot(iisopn)
      cnn(kisopo2)	=stot(iisopo2)
      cnn(kterp)        =stot(iterp)
      cnn(ksv3)         =stot(isv3)
      cnn(ksv4)         =stot(isv4)
      cnn(ksv5)         =stot(isv5)
      cnn(ksv6)         =stot(isv6) 
      endif

      return
      end





c***********************************************************************
c subroutine MapGas_Com: maps cnn to and fro stot for the common 
c                        gas-phase mechanism.
c
c nomenclature:
c cnn       = full species concentration array.
c stot      = subset of cnn. species concentration array to be supplied to
c             lsodes. length of stot depends on the selected mechanism
c iregime   = selected chemical regime (1-6)
c imap      = 0 : map cnn to stot
c           = 1 : map stot to cnn
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine MapGas_Com(stot,imap)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)

      emit(ih2so4)	=emission(kh2so4)
      emit(ihno3)	=emission(khno3)
      emit(ihcl)	=emission(khcl)
      emit(inh3)	=emission(knh3)
      emit(ino)		=emission(kno)
      emit(ino2)	=emission(kno2)
      emit(ino3)	=emission(kno3)
      emit(in2o5)	=emission(kn2o5)
      emit(ihono)	=emission(khono)
      emit(ihno4)	=emission(khno4)
      emit(io3)		=emission(ko3)
      emit(io1d)	=emission(ko1d)
      emit(io3p)	=emission(ko3p)
      emit(ioh)		=emission(koh)
      emit(iho2)	=emission(kho2)
      emit(ih2o2)	=emission(kh2o2)
      emit(ico)		=emission(kco)
      emit(iso2)	=emission(kso2)
      emit(ich4)	=emission(kch4)
      emit(ic2h6)	=emission(kc2h6)
      emit(ich3o2)	=emission(kch3o2)
      emit(iethp)	=emission(kethp)
      emit(ihcho)	=emission(khcho)
      emit(ich3oh)	=emission(kch3oh)
      emit(ianol)	=emission(kanol)
      emit(ich3ooh)	=emission(kch3ooh)
      emit(iethooh)	=emission(kethooh)
      emit(iald2)	=emission(kald2)
      emit(ihcooh)	=emission(khcooh)
      emit(ircooh)	=emission(krcooh)
      emit(ic2o3)	=emission(kc2o3)
      emit(ipan)	=emission(kpan)

      if(imap.eq.0)then    ! map cnn into stot
      stot(ih2so4)	=cnn(kh2so4)
      stot(ihno3)	=cnn(khno3)
      stot(ihcl)	=cnn(khcl)
      stot(inh3)	=cnn(knh3)
      stot(ino)		=cnn(kno)
      stot(ino2)	=cnn(kno2)
      stot(ino3)	=cnn(kno3)
      stot(in2o5)	=cnn(kn2o5)
      stot(ihono)	=cnn(khono)
      stot(ihno4)	=cnn(khno4)
      stot(io3)		=cnn(ko3)
      stot(io1d)	=cnn(ko1d)
      stot(io3p)	=cnn(ko3p)
      stot(ioh)		=cnn(koh)
      stot(iho2)	=cnn(kho2)
      stot(ih2o2)	=cnn(kh2o2)
      stot(ico)		=cnn(kco)
      stot(iso2)	=cnn(kso2)
      stot(ich4)	=cnn(kch4)
      stot(ic2h6)	=cnn(kc2h6)
      stot(ich3o2)	=cnn(kch3o2)
      stot(iethp)	=cnn(kethp)
      stot(ihcho)	=cnn(khcho)
      stot(ich3oh)	=cnn(kch3oh)
      stot(ianol)	=cnn(kanol)
      stot(ich3ooh)	=cnn(kch3ooh)
      stot(iethooh)	=cnn(kethooh)
      stot(iald2)	=cnn(kald2)
      stot(ihcooh)	=cnn(khcooh)
      stot(ircooh)	=cnn(krcooh)
      stot(ic2o3)	=cnn(kc2o3)
      stot(ipan)	=cnn(kpan)
c
      else                 ! map stot back into cnn
      cnn(kh2so4)	=stot(ih2so4)
      cnn(khno3)	=stot(ihno3)
      cnn(khcl)		=stot(ihcl)
      cnn(knh3)		=stot(inh3)
      cnn(kno)		=stot(ino)
      cnn(kno2)		=stot(ino2)
      cnn(kno3)		=stot(ino3)
      cnn(kn2o5)	=stot(in2o5)
      cnn(khono)	=stot(ihono)
      cnn(khno4)	=stot(ihno4)
      cnn(ko3)		=stot(io3)
      cnn(ko1d)		=stot(io1d)
      cnn(ko3p)		=stot(io3p)
      cnn(koh)		=stot(ioh)
      cnn(kho2)		=stot(iho2)
      cnn(kh2o2)	=stot(ih2o2)
      cnn(kco)		=stot(ico)
      cnn(kso2)		=stot(iso2)
      cnn(kch4)		=stot(ich4)
      cnn(kc2h6)	=stot(ic2h6)
      cnn(kch3o2)	=stot(ich3o2)
      cnn(kethp)	=stot(iethp)
      cnn(khcho)	=stot(ihcho)
      cnn(kch3oh)	=stot(ich3oh)
      cnn(kanol)	=stot(ianol)
      cnn(kch3ooh)	=stot(ich3ooh)
      cnn(kethooh)	=stot(iethooh)
      cnn(kald2)	=stot(iald2)
      cnn(khcooh)	=stot(ihcooh)
      cnn(krcooh)	=stot(ircooh)
      cnn(kc2o3)	=stot(ic2o3)
      cnn(kpan)		=stot(ipan)
      endif

      return
      end







c***********************************************************************
c subroutine MapGas_Mar: maps cnn to and fro stot for the marine
c                        gas-phase mechanism.
c
c nomenclature:
c cnn       = full species concentration array.
c stot      = subset of cnn. species concentration array to be supplied to
c             lsodes. length of stot depends on the selected mechanism
c iregime   = selected chemical regime (1-6)
c imap      = 0 : map cnn to stot
c           = 1 : map stot to cnn
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine MapGas_Mar(stot,imap)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)

      emit(idms)	=emission(kdms)
      emit(imsa)	=emission(kmsa)
      emit(idmso)	=emission(kdmso)
      emit(idmso2)	=emission(kdmso2)
      emit(ich3so2h)	=emission(kch3so2h)
      emit(ich3sch2oo)	=emission(kch3sch2oo)
      emit(ich3so2)	=emission(kch3so2)
      emit(ich3so3)	=emission(kch3so3)
      emit(ich3so2oo)	=emission(kch3so2oo)
      emit(ich3so2ch2oo)=emission(kch3so2ch2oo)
      emit(isulfhox)	=emission(ksulfhox)

      if(imap.eq.0)then    ! map cnn into stot
      stot(idms)	=cnn(kdms)
      stot(imsa)	=cnn(kmsa)
      stot(idmso)	=cnn(kdmso)
      stot(idmso2)	=cnn(kdmso2)
      stot(ich3so2h)	=cnn(kch3so2h)
      stot(ich3sch2oo)	=cnn(kch3sch2oo)
      stot(ich3so2)	=cnn(kch3so2)
      stot(ich3so3)	=cnn(kch3so3)
      stot(ich3so2oo)	=cnn(kch3so2oo)
      stot(ich3so2ch2oo)=cnn(kch3so2ch2oo)
      stot(isulfhox)	=cnn(ksulfhox)

      else                 ! map stot back into cnn
      cnn(kdms)		=stot(idms)
      cnn(kmsa)		=stot(imsa)
      cnn(kdmso)	=stot(idmso)
      cnn(kdmso2)	=stot(idmso2)
      cnn(kch3so2h)	=stot(ich3so2h)
      cnn(kch3sch2oo)	=stot(ich3sch2oo)
      cnn(kch3so2)	=stot(ich3so2)
      cnn(kch3so3)	=stot(ich3so3)
      cnn(kch3so2oo)	=stot(ich3so2oo)
      cnn(kch3so2ch2oo)	=stot(ich3so2ch2oo)
      cnn(ksulfhox)	=stot(isulfhox)
      endif

      return
      end







c***********************************************************************
c subroutine MapGasSpecies: maps cnn to and fro stot for the selected 
c                           gas-phase mechanism.
c
c nomenclature:
c cnn       = full species concentration array.
c stot      = subset of cnn. species concentration array to be supplied to
c             lsodes. length of stot depends on the selected mechanism
c iregime   = selected chemical regime (1-6)
c imap      = 0 : map cnn to stot
c           = 1 : map stot to cnn
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine MapGasSpecies(stot,imap)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)

      goto (1,2,3,4,5,6), iregime
c
c
1     call mapgas_com(stot,imap)
      return
c
c
2     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      return
c
c
3     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      return
c
c
4     call mapgas_com(stot,imap)
      call mapgas_mar(stot,imap)
      return
c
c 
5     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_mar(stot,imap)
      return
c
c
6     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      call mapgas_mar(stot,imap)
      return
c
      end







c***********************************************************************
c subroutine MapGas_Urb: maps cnn to and fro stot for the urban 
c                        gas-phase mechanism.
c
c nomenclature:
c cnn       = full species concentration array.
c stot      = subset of cnn. species concentration array to be supplied to
c             lsodes. length of stot depends on the selected mechanism
c iregime   = selected chemical regime (1-6)
c imap      = 0 : map cnn to stot
c           = 1 : map stot to cnn
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine MapGas_Urb(stot,imap)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)

      emit(ipar)	=emission(kpar)
      emit(iaone)	=emission(kaone)
      emit(imgly)	=emission(kmgly)
      emit(ieth)	=emission(keth)
      emit(iolet)	=emission(kolet)
      emit(iolei)	=emission(kolei)
      emit(itol)	=emission(ktol)
      emit(ixyl)	=emission(kxyl)
      emit(icres)	=emission(kcres)
      emit(ito2)	=emission(kto2)
      emit(icro)	=emission(kcro)
      emit(iopen)	=emission(kopen)
      emit(ionit)	=emission(konit)
      emit(irooh)	=emission(krooh)
      emit(iro2)	=emission(kro2)
      emit(iano2)	=emission(kano2)
      emit(inap)	=emission(knap)
      emit(ixo2)	=emission(kxo2)
      emit(ixpar)	=emission(kxpar)

      if(imap.eq.0)then    ! map cnn into stot
      stot(ipar)	=cnn(kpar)
      stot(iaone)	=cnn(kaone)
      stot(imgly)	=cnn(kmgly)
      stot(ieth)	=cnn(keth)
      stot(iolet)	=cnn(kolet)
      stot(iolei)	=cnn(kolei)
      stot(itol)	=cnn(ktol)
      stot(ixyl)	=cnn(kxyl)
      stot(icres)	=cnn(kcres)
      stot(ito2)	=cnn(kto2)
      stot(icro)	=cnn(kcro)
      stot(iopen)	=cnn(kopen)
      stot(ionit)	=cnn(konit)
      stot(irooh)	=cnn(krooh)
      stot(iro2)	=cnn(kro2)
      stot(iano2)	=cnn(kano2)
      stot(inap)	=cnn(knap)
      stot(ixo2)	=cnn(kxo2)
      stot(ixpar)	=cnn(kxpar)
      stot(isv1)        =cnn(ksv1)
      stot(isv2)        =cnn(ksv2)
c
      else                 ! map stot back into cnn
      cnn(kpar)		=stot(ipar)
      cnn(kaone)	=stot(iaone)
      cnn(kmgly)	=stot(imgly)
      cnn(keth)		=stot(ieth)
      cnn(kolet)	=stot(iolet)
      cnn(kolei)	=stot(iolei)
      cnn(ktol)		=stot(itol)
      cnn(kxyl)		=stot(ixyl)
      cnn(kcres)	=stot(icres)
      cnn(kto2)		=stot(ito2)
      cnn(kcro)		=stot(icro)
      cnn(kopen)	=stot(iopen)
      cnn(konit)	=stot(ionit)
      cnn(krooh)	=stot(irooh)
      cnn(kro2)		=stot(iro2)
      cnn(kano2)	=stot(iano2)
      cnn(knap)		=stot(inap)
      cnn(kxo2)		=stot(ixo2)
      cnn(kxpar)	=stot(ixpar)
      cnn(ksv1)         =stot(isv1)
      cnn(ksv2)         =stot(isv2)
      endif

      return
      end







c***********************************************************************
      subroutine ode_bio
      include 'chm.inc'
      include 'gas.inc'
c
c   
      p_bio(ino)= 0.0
      d_bio(ino)= r_bio(8)+r_bio(9)+r_bio(10)
c
      p_bio(ino2)= .91*r_bio(8)+1.2*r_bio(9)+r_bio(10)+0.47*r_bio(20)
      d_bio(ino2)= 0.0
c
      p_bio(ino3)= 0.0
      d_bio(ino3)= r_bio(3)+r_bio(7)+r_bio(20)
c
      p_bio(ihno3)= .07*r_bio(7)
      d_bio(ihno3)= 0.0
c
      p_bio(io3)= 0.0
      d_bio(io3)= r_bio(2)+r_bio(6)+r_bio(19)
c
      p_bio(ioh)= .27*r_bio(2)+.27*r_bio(6)+0.57*r_bio(19)
      d_bio(ioh)= r_bio(1)+r_bio(5)+r_bio(18)
c
      p_bio(iho2)= .07*r_bio(2)+.33*r_bio(4)+.1*r_bio(6)
     &        +.93*r_bio(7)+.91*r_bio(8)+.8*r_bio(9)+r_bio(10)
     &        +.75*r_bio(18)+.07*r_bio(19)+.28*r_bio(20) 
      d_bio(iho2)= r_bio(11)+r_bio(12)+r_bio(13)
c
      p_bio(ih2o2)= 0.0
      d_bio(ih2o2)= 0.0
c
      p_bio(ico)= .07*r_bio(2)+.33*r_bio(4)+.16*r_bio(6)
     &       +.64*r_bio(7)+.59*r_bio(10)+0.001*r_bio(19)
      d_bio(ico)= 0.0
c
      p_bio(ihcho)= .6*r_bio(2)+.2*r_bio(4)+.15*r_bio(6)
     &         +.28*r_bio(7)+.63*r_bio(8)+.25*r_bio(10)
     &         +.28*r_bio(18)+.24*r_bio(19)
      d_bio(ihcho)= 0.0
c
      p_bio(iald2)= .15*r_bio(2)+.07*r_bio(4)+.02*r_bio(6)
     &         +.28*r_bio(7)+.8*r_bio(9)+.55*r_bio(10)+r_bio(15)
     &         +.5*r_bio(16)+.15*r_bio(17)+.47*r_bio(18)
     &         +.21*r_bio(19)+0.47*r_bio(20)
      d_bio(iald2)= 0.0
c
      p_bio(ipar)= 1.86*r_bio(7)+.18*r_bio(8)+1.6*r_bio(9)+
     &             2.*r_bio(12)+2.*r_bio(15)+0.51*r_bio(17)+
     &             7.*r_bio(19) 
      d_bio(ipar)= 0.0
c
      p_bio(iaone)= .03*r_bio(4)+.09*r_bio(6)+.63*r_bio(10)
     &         +.5*r_bio(16)
      d_bio(iaone)= 0.0
c
      p_bio(imgly)= .85*r_bio(6)+.34*r_bio(10)
      d_bio(imgly)= 0.0
c
      p_bio(ionit)= .93*r_bio(7)+.09*r_bio(8)+.8*r_bio(9)+r_bio(12)
     &         +r_bio(15)+0.53*r_bio(20)
      d_bio(ionit)= 0.0
c
      p_bio(ircooh)= .39*r_bio(2)+.46*r_bio(6)
      d_bio(ircooh)= 0.0
c
      p_bio(irooh)= r_bio(11)+r_bio(13)
      d_bio(irooh)= 0.0
c
      p_bio(ich3o2)= .7*r_bio(4)+.05*r_bio(6)
      d_bio(ich3o2)= 0.0
c
      p_bio(ic2o3)= .2*r_bio(2)+.97*r_bio(4)+.5*r_bio(5)
     &         +.11*r_bio(6)+.07*r_bio(7)
      d_bio(ic2o3)= 0.0
c
      p_bio(ixo2)= .08*r_bio(1)+.2*r_bio(2)+.2*r_bio(5)+.07*r_bio(6)
     &        +.93*r_bio(7)+1.25*r_bio(18)+.76*r_bio(19)+.76*r_bio(20)
      d_bio(ixo2)= 0.0
c
      p_bio(iisop)= 0.0
      d_bio(iisop)= r_bio(1)+r_bio(2)+r_bio(3)
c
      p_bio(iisoprd)= .65*r_bio(2)+.91*r_bio(8)+.2*r_bio(9)+r_bio(14)
      d_bio(iisoprd)= r_bio(4)+r_bio(5)+r_bio(6)+r_bio(7)
c
      p_bio(iisopp)= r_bio(1)
      d_bio(iisopp)= r_bio(8)+r_bio(11)+r_bio(14)
c
      p_bio(iisopn)= r_bio(3)
      d_bio(iisopn)= r_bio(9)+r_bio(12)+r_bio(15)
c
      p_bio(iisopo2)= .5*r_bio(5)
      d_bio(iisopo2)= r_bio(10)+r_bio(13)+r_bio(16)
c
      p_bio(io1d) = 0.0
      d_bio(io1d) = r_bio(17)
c
      p_bio(isv3) = .110681*r_bio(17)+.063304*r_bio(18)
     &             +.10998*r_bio(19)+.092287*r_bio(20)
      d_bio(isv3) = 0.0
c
      p_bio(isv4) = .4172*r_bio(17)+.324763*r_bio(18)  
     *             +.389277*r_bio(19)+.395069*r_bio(20)
      d_bio(isv4) = 0.0
c
      p_bio(isv5) = .232*r_bio(1)
      d_bio(isv5) = 0.0
c
      p_bio(isv6) = .0288*r_bio(1)
      d_bio(isv6) = 0.0
c
      p_bio(iterp) = 0.0
      d_bio(iterp) = r_bio(17)+r_bio(18)+r_bio(19)+r_bio(20)
c      
      return
      end








c***********************************************************************
      subroutine ode_com
      include 'chm.inc'
      include 'gas.inc'

      p_com(ih2so4)= r_com(45)                                                          
      d_com(ih2so4)= 0.0                                                            
c
      p_com(ihno3)= r_com(24)+.3*r_com(41)+r_com(42)+r_com(42)                                
     &         +r_com(52)+r_com(68)                                                           
      d_com(ihno3)= r_com(4)+r_com(27)                                                      
c
      p_com(ihcl)= 0.0                                                              
      d_com(ihcl)= 0.0                                                              
c
      p_com(inh3)= 0.0                                                              
      d_com(inh3)= 0.0                                                              
c
      p_com(ino)= r_com(1)+0.11*r_com(2)+r_com(3)+r_com(15)+r_com(38)                                   
      d_com(ino)= r_com(17)+r_com(18)+r_com(23)+r_com(33)+r_com(37)                               
     &       +r_com(57)+r_com(58)+r_com(71)                                                       
c
      p_com(ino2)= 0.89*r_com(2)+r_com(4)+r_com(5)
     &        +r_com(6)+r_com(17)+r_com(18)                             
     &        +r_com(25)+r_com(26)+r_com(28)
     &        +r_com(33)+r_com(36)+r_com(37)                              
     &        +r_com(37)+r_com(38)+r_com(40)+r_com(40)+.7*r_com(41)                                 
     &        +r_com(43)+r_com(57)+r_com(58)+r_com(59)
     &        +r_com(60)+r_com(70)                              
     &        +r_com(71)+r_com(72)                                                      
      d_com(ino2)= r_com(1)+r_com(15)+r_com(16)+r_com(19)
     &        +r_com(24)+r_com(34)                               
     &        +r_com(35)+r_com(38)+r_com(39)+r_com(69)                                          
c
      p_com(ino3)= r_com(6)+r_com(16)+r_com(19)+r_com(27)+r_com(43)                                     
      d_com(ino3)= r_com(2)+r_com(25)+r_com(37)+r_com(38)
     &        +r_com(39)+r_com(40)                               
     &        +r_com(40)+r_com(41)+r_com(52)+r_com(59)
     &        +r_com(60)+r_com(68)                              
     &        +r_com(72)                                                            
c
      p_com(in2o5)= r_com(39)                                                           
      d_com(in2o5)= r_com(6)+r_com(42)+r_com(43)                                                
c
      p_com(ihono)= r_com(23)+r_com(35)                                                     
      d_com(ihono)= r_com(3)+r_com(26)                                                      
c
      p_com(ihno4)= r_com(34)                                                           
      d_com(ihno4)= r_com(5)+r_com(28)+r_com(36)                                                
c
      p_com(io3)= r_com(13)+.4*r_com(73)                                                    
      d_com(io3)= r_com(7)+r_com(8)+r_com(14)+r_com(18)
     &       +r_com(19)+r_com(20)+r_com(21)                                                             
c
      p_com(io1d)= r_com(8)                                                             
      d_com(io1d)= r_com(10)+r_com(11)+r_com(12)                                                
c
      p_com(io3p)= r_com(1)+0.89*r_com(2)+r_com(7)
     &        +r_com(10)+r_com(11)                                  
      d_com(io3p)= r_com(13)+r_com(14)+r_com(15)+r_com(16)+r_com(17)                                    
c
      p_com(ioh)= r_com(3)+r_com(4)+2*r_com(9)
     &       +2*r_com(12)+r_com(21)+r_com(33)                              
     &       +.7*r_com(41)+r_com(53)+r_com(54)
     &       +.3*r_com(55)+.5*r_com(56)                            
      d_com(ioh)= r_com(20)+r_com(22)+r_com(23)
     &       +r_com(24)+r_com(25)+r_com(26)                               
     &       +r_com(27)+r_com(28)+r_com(29)+r_com(30)
     &       +r_com(44)+r_com(45)                               
     &       +r_com(46)+r_com(47)+r_com(48)
     &       +r_com(51)+r_com(55)+r_com(56)                               
     &       +r_com(65)+r_com(67)                                                       
c
      p_com(iho2)= r_com(5)+r_com(20)+r_com(22)
     &        +r_com(25)+r_com(30)+r_com(36)                               
     &        +r_com(44)+r_com(45)+r_com(48)
     &        +2*r_com(49)+r_com(51)                                  
     &        +r_com(52)+r_com(53)+r_com(54)
     &        +r_com(57)+r_com(58)+r_com(59)                              
     &        +r_com(60)+.32*r_com(63)+.6*r_com(64)
     &        +r_com(65)+r_com(66)                             
      d_com(iho2)= r_com(21)+r_com(29)+r_com(31)
     &        +r_com(31)+r_com(32)+r_com(32)                              
     &        +r_com(33)+r_com(34)+r_com(35)
     &        +r_com(41)+r_com(61)+r_com(62)                              
     &        +r_com(73)                                                            
c
      p_com(ih2o2)= r_com(31)+r_com(32)                                                     
      d_com(ih2o2)= r_com(9)+r_com(30)                                                      
c
      p_com(ico)= r_com(49)+r_com(50)+r_com(51)
     &        +r_com(52)+r_com(66)                                     
      d_com(ico)= r_com(44)                                                             
c
      p_com(iso2)= 0.0                                                              
      d_com(iso2)= r_com(45)                                                            
c
      p_com(ich4)= 0.0                                                              
      d_com(ich4)= r_com(46)                                                            
c
      p_com(ic2h6)= .2*r_com(64)                                                        
      d_com(ic2h6)= r_com(47)                                                           
c
      p_com(ich3o2)= r_com(46)+.7*r_com(55)+r_com(66)
     &          +r_com(71)+r_com(72)                               
     &          +r_com(74)                                                          
      d_com(ich3o2)= r_com(57)+r_com(59)+r_com(61)+r_com(63)                                        
c
      p_com(iethp)= r_com(47)+.5*r_com(56)                                                  
      d_com(iethp)= r_com(58)+r_com(60)+r_com(62)+r_com(64)                                         
c
      p_com(ihcho)= r_com(48)+r_com(53)+.3*r_com(55)
     &         +r_com(57)+r_com(59)                                
     &         +.66*r_com(63)                                                       
      d_com(ihcho)= r_com(49)+r_com(50)+r_com(51)+r_com(52)                                         
c
      p_com(ich3oh)= .34*r_com(63)                                                      
      d_com(ich3oh)= r_com(48)                                                          
c
      p_com(ianol)= 0.0                                                             
      d_com(ianol)= r_com(65)                                                           
c
      p_com(ich3ooh)= r_com(61)                                                         
      d_com(ich3ooh)= r_com(53)+r_com(55)                                                   
c
      p_com(iethooh)= r_com(62)                                                         
      d_com(iethooh)= r_com(54)+r_com(56)                                                   
c
      p_com(iald2)= r_com(54)+.5*r_com(56)+r_com(58)
     &         +r_com(60)+.8*r_com(64)                             
     &         +r_com(65)                                                           
      d_com(iald2)= r_com(66)+r_com(67)+r_com(68)                                               
c
      p_com(ihcooh)= 0.0                                                            
      d_com(ihcooh)= 0.0                                                            
c
      p_com(ircooh)= .4*r_com(73)                                                       
      d_com(ircooh)= 0.0                                                            
c
      p_com(ic2o3)= r_com(67)+r_com(68)+r_com(70)                                               
      d_com(ic2o3)= r_com(69)+r_com(71)+r_com(72)
     &        +r_com(73)+r_com(74)                                   
c
      p_com(ipan)= r_com(69)                                                            
      d_com(ipan)= r_com(70)                                                            
c


      return
      end






c***************************************************************************
c subroutine ode_gas
c
c purpose: computes time derivatives of species concentrations ds/dt = sdot(*).
c          calls ode_com, ode_urb, ode_bio, ode_mar depending on the
c          chemical regime (iregime)
c
c author: Rahul A. Zaveri
c
c---------------------------------------------------------------------------
c
      subroutine ode_gas(ntot,tt,s,sdot)
      include 'chm.inc'
      include 'gas.inc'
c
      dimension s(nmax),sdot(nmax)
c
      do i=1,nrxn_com
        r_com(i) = 0.
      enddo

      do i=1,nrxn_urb
        r_urb(i) = 0.
      enddo

      do i=1,nrxn_bio
        r_bio(i) = 0.
      enddo

      do i=1,nrxn_mar
        r_mar(i) = 0.
      enddo

      do i=1,nrxn_het
        r_het(i) = 0.
      enddo


      call GasRates(s)
c
      do i=1,ngas_max
        p_com(i) = 0.
        p_urb(i) = 0.
        p_bio(i) = 0.
        p_mar(i) = 0.
        p_het(i) = 0.

        d_com(i) = 0.
        d_urb(i) = 0.
        d_bio(i) = 0.
        d_mar(i) = 0.
        d_het(i) = 0.
      enddo
c    
c
      goto (1,2,3,4,5,6), iregime

1     call ode_com
      call ode_het

      do i=1,nreg1
        sdot(i) = real( dble(p_com(i)+p_het(i)) - 
     &                  dble(d_com(i)+d_het(i)) ) 
!     &          + emit(i)
      enddo

      return
c
c----------------------------------------------------------
2     call ode_com
      call ode_urb
      call ode_het

      do i=1,nreg2
        sdot(i) = real( dble(p_com(i)+p_urb(i)+p_het(i)) - 
     &                  dble(d_com(i)+d_urb(i)+d_het(i)) )
!     &          + emit(i)
      enddo

      return
c
c----------------------------------------------------------
3     call ode_com
      call ode_urb
      call ode_bio
      call ode_het

      do i=1,nreg3
        sdot(i) = real( dble(p_com(i)+p_urb(i)+p_bio(i)+p_het(i)) -
     &                  dble(d_com(i)+d_urb(i)+d_bio(i)+d_het(i)) )
!     &          + emit(i)
      enddo

      return
c
c----------------------------------------------------------
4     call ode_com
      call ode_mar
      call ode_het

      do i=1,nreg4
        sdot(i) = real( dble(p_com(i)+p_mar(i)+p_het(i)) -
     &                  dble(d_com(i)+d_mar(i)+d_het(i)) )
!     &          + emit(i)
      enddo

      return
c
c----------------------------------------------------------
5     call ode_com
      call ode_urb
      call ode_mar
      call ode_het

      do i=1,nreg5
        sdot(i) = real( dble(p_com(i)+p_urb(i)+p_mar(i)+p_het(i)) -
     &                  dble(d_com(i)+d_urb(i)+d_mar(i)+d_het(i)) )
!     &          + emit(i)
      enddo

      return
c
c----------------------------------------------------------
6     call ode_com
      call ode_urb
      call ode_bio
      call ode_mar
      call ode_het

      do i=1,nreg6
        sdot(i) = real( dble(p_com(i)+p_urb(i)+p_bio(i)+p_mar(i)+
     &                       p_het(i)) -
     &                  dble(d_com(i)+d_urb(i)+d_bio(i)+d_mar(i)+
     &                       d_het(i)) )
!     &          + emit(i)
      enddo

      return

      end






c***********************************************************************
      subroutine ode_het
      include 'chm.inc'
      include 'gas.inc'

      p_het(io3)= 0.0
      d_het(io3)= r_het(1)

      p_het(in2o5)= 0.0
      d_het(in2o5)= r_het(4)
      
      p_het(ihno3)= 2.*r_het(4)
      d_het(ihno3)= 0.0
      
      p_het(ino3) = 0.0
      d_het(ino3) = r_het(7)
      
      p_het(ino2)  = r_het(7)
      d_het(ino2)  = 0.0
            
      return
      end






c***********************************************************************
      subroutine ode_mar
      include 'chm.inc'
      include 'gas.inc'
c
c
      a = 5.e+5/(5.e+5 + o2*3.e-12)
      b = 1.5e+7/(1.5e+7 + o2*1.2e-12)
c
c
      p_mar(ino)= r_mar(17)   
      d_mar(ino)= r_mar(5)+r_mar(9)+r_mar(24)+r_mar(28)
c
      p_mar(ino2)= r_mar(5)+r_mar(9)+r_mar(24)   
      d_mar(ino2)= r_mar(17)+r_mar(27)
c
      p_mar(ino3)=   0.0
      d_mar(ino3)= r_mar(2)+r_mar(12)   
c
      p_mar(ihono)= r_mar(28)   
      d_mar(ihono)=   0.0
c
      p_mar(ihno3)= r_mar(2)+r_mar(12)+r_mar(27)   
      d_mar(ihno3)=   0.0
c
      p_mar(io3)=   0.0
      d_mar(io3)= r_mar(18)   
c
      p_mar(io3p)=   0.0
      d_mar(io3p)= r_mar(3)    
c
      p_mar(ioh)= r_mar(19)   
      d_mar(ioh)= r_mar(1)+r_mar(4)+r_mar(7)+r_mar(8)+
     &            r_mar(14)+r_mar(21)   
c
      p_mar(iho2)= (1.-a)*r_mar(4)+r_mar(6)+(1.-b)*r_mar(7)+
     &             r_mar(10)+r_mar(20)+r_mar(25)+r_mar(30)   
      d_mar(iho2)= r_mar(11)+r_mar(19)+r_mar(29)   
c
      p_mar(ih2o2)= r_mar(11)   
      d_mar(ih2o2)=   0.0
c
      p_mar(iso2)= r_mar(16)   
      d_mar(iso2)=   0.0 + r_mar(33)
c
      p_mar(ih2so4)= r_mar(26)   
      d_mar(ih2so4)=   0.0
c
      p_mar(ich3o2)= r_mar(3)+a*r_mar(4)+b*r_mar(7)+r_mar(16)+
     &               r_mar(26)   
      d_mar(ich3o2)= r_mar(6)+r_mar(10)+r_mar(13)+r_mar(20)+
     &               r_mar(25)   
c
      p_mar(ich3ooh)= r_mar(13)   
      d_mar(ich3ooh)=   0.0
c
      p_mar(ihcho)= r_mar(5)+r_mar(6)+r_mar(6)+r_mar(9)+
     &              r_mar(10)+r_mar(10)+r_mar(20)+r_mar(25)   
      d_mar(ihcho)= r_mar(30)   
c
      p_mar(idms)= 0.0
      d_mar(idms)= r_mar(1)+r_mar(2)+r_mar(3)+r_mar(4)    
c
      p_mar(imsa)= r_mar(15)+r_mar(21)+r_mar(27)+r_mar(28)+
     &             r_mar(29)+r_mar(30)   
      d_mar(imsa)=   0.0
c
      p_mar(idmso)= (1.-a)*r_mar(4)    
      d_mar(idmso)= r_mar(7) + r_mar(34)
c
      p_mar(idmso2)= (1.-b)*r_mar(7)    
      d_mar(idmso2)= r_mar(8) + r_mar(35)
c
      p_mar(ich3so2h)= b*r_mar(7)    
      d_mar(ich3so2h)= r_mar(11)+r_mar(12)+r_mar(13)+r_mar(14)+
     &                 r_mar(15)   
c
      p_mar(ich3sch2oo)= r_mar(1)+r_mar(2)
      d_mar(ich3sch2oo)= r_mar(5)+r_mar(6)+r_mar(31)+2.*r_mar(32)
c
      p_mar(ich3so2)= r_mar(3)+a*r_mar(4)+r_mar(5)+r_mar(6)+
     &                r_mar(9)+r_mar(10)+r_mar(11)+r_mar(12)+
     &                r_mar(13)+r_mar(14)+r_mar(15)+r_mar(23)+
     &                r_mar(31)+1.85*r_mar(32)
      d_mar(ich3so2)= r_mar(16)+r_mar(17)+r_mar(18)+r_mar(19)+
     &                r_mar(20)+r_mar(21)+r_mar(22)+r_mar(31)
c
      p_mar(ich3so3)= r_mar(17)+r_mar(18)+r_mar(19)+r_mar(20)+
     &                r_mar(24)+r_mar(25)+r_mar(31)
      d_mar(ich3so3)= r_mar(15)+r_mar(26)+r_mar(27)+r_mar(28)+
     &                r_mar(29)+r_mar(30)   
c
      p_mar(ich3so2oo)= r_mar(22)   
      d_mar(ich3so2oo)= r_mar(23)+r_mar(24)+r_mar(25)   
c
      p_mar(ich3so2ch2oo)= r_mar(8)    
      d_mar(ich3so2ch2oo)= r_mar(9)+r_mar(10)
c
      p_mar(isulfhox)= 0.15*r_mar(32)
      d_mar(isulfhox)= 0.0
c

c
      return
      end







c***********************************************************************
      subroutine ODEsolver(ntot,stot,t_in,t_out)

      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)	! local species array

c lsodes parameters
      parameter(itoler = 2, itask = 1, iopt = 1, mf = 222)
      parameter(nnz = ngas_max*ngas_max)
      parameter(lwm = 3*ngas_max*ngas_max + 12*ngas_max)
      parameter(nrdim = 20 + 9*ngas_max + lwm)
      parameter(nidim = 31 + ngas_max + ngas_max*ngas_max)
      dimension rwork(nrdim),iwork(nidim), atol(ngas_max)
      external ode_gas,jac

c set LSODE parameters...
c relative tolerance
      rtol=1.e-6
c
c absolute tolerances for gas-phase species
      do i=1,ngas_max
        atol(i)=1.e+1	! [molec/cc]
      enddo     

      iwork(6) = 1000
      iwork(7) = 1
      istate   = 1
      rwork(6) = dt_sec
      iflag=0  ! juanxiong he
      if(iflag .eq. 0)then
      call xsetf(iflag)
      endif

      call LSODES(ode_gas,ntot,stot,t_in,t_out,itoler,
     & rtol,atol,itask,istate,iopt,rwork,nrdim,iwork,nidim,jac,mf)

c
c error message from LSODES...
      if(istate.lt.0)then
czifa        write(6,*)'WARNING from LSODES: istate=',istate
czifa        write(6,*)'Calling LSODE to solve this step'
c
c reset time
        t_in = told_sec
        t_out= tcur_sec
        istate = 1
c
        call MapGasSpecies(stot,0)	! map cnn into stot

        call LSODE(ode_gas,ntot,stot,t_in,t_out,itoler,
     &  rtol,atol,itask,istate,iopt,rwork,nrdim,iwork,nidim,jac,22)

c error message from LSODE...
        if(istate.lt.0)then
c        write(6,*)'WARNING zifa from LSODE: istate=',istate
c        stop
        else
c!zifa        write(6,*)'successful solution from LSODE' 
        endif
c
      else
c        write(6,*)'successful solution from LSODES'
      endif
c
c successful solution...
c set negative concentration values, if any, to zero...
      if(istate.eq.2)then
        do i=1,ntot
        stot(i)=max(0.,stot(i))
        enddo
      endif

      return
      end











c***********************************************************************
      subroutine ode_urb
      include 'chm.inc'
      include 'gas.inc'
c
      p_urb(ihno3)= r_urb(6)+r_urb(19)                                                      
      d_urb(ihno3)= 0.0                                                             
c
      p_urb(ino)= 0.0                                                               
      d_urb(ino)= r_urb(17)+r_urb(28)+r_urb(29)+r_urb(30)+r_urb(31)                                     
c
      p_urb(ino2)= .95*r_urb(17)+r_urb(27)+.84*r_urb(28)+r_urb(29)                                  
     &        +1.5*r_urb(30)+r_urb(31)+r_urb(32)+r_urb(33)
     &        +1.5*r_urb(34)+r_urb(35)+.5*r_urb(42)                                                   
      d_urb(ino2)= r_urb(20)                                                            
c
      p_urb(ino3)= 0.0                                                              
      d_urb(ino3)= r_urb(6)+r_urb(13)+r_urb(14)+r_urb(19)
     &        +r_urb(32)+r_urb(33)+r_urb(34)+r_urb(35)                                                      
c
      p_urb(io3)= 0.0                                                               
      d_urb(io3)= r_urb(7)+r_urb(9)+r_urb(10)+r_urb(23)                                             
c
      p_urb(ioh)= .12*r_urb(7)+.33*r_urb(9)+.6*r_urb(10)
     &       +.08*r_urb(23)+r_urb(24)+.23*r_urb(25)                                                   
      d_urb(ioh)= r_urb(1)+r_urb(3)+r_urb(5)
     &       +r_urb(8)+r_urb(11)+r_urb(12)+r_urb(15)                             
     &       +r_urb(16)+r_urb(18)+r_urb(21)+r_urb(25)+r_urb(26)                                     
c
      p_urb(iho2)= r_urb(4)+.22*r_urb(7)+r_urb(8)
     &        +.26*r_urb(9)+.22*r_urb(10)                            
     &        +r_urb(11)+r_urb(12)+.2*r_urb(15)
     &        +.55*r_urb(16)+.95*r_urb(17)                         
     &        +.6*r_urb(18)+2*r_urb(21)+r_urb(22)+.76*r_urb(23)                                 
     &        +.9*r_urb(24)+.9*r_urb(27)+.76*r_urb(28)                             
     &        +.5*r_urb(30)+.9*r_urb(32)+.5*r_urb(34)
     &        +.54*r_urb(40)                                      
      d_urb(iho2)= r_urb(36)+r_urb(37)+r_urb(38)+r_urb(39)                                          
c
      p_urb(ico)= r_urb(4)+r_urb(6)+.24*r_urb(7)
     &       +.31*r_urb(9)+.3*r_urb(10)                              
     &       +2*r_urb(21)+r_urb(22)+.69*r_urb(23)                                           
      d_urb(ico)= 0.0                                                               
c
      p_urb(ich3o2)= r_urb(2)+.07*r_urb(9)+.1*r_urb(10)                                         
      d_urb(ich3o2)= 0.0                                                            
c
      p_urb(iethp)= .06*r_urb(9)+.05*r_urb(10)+.1*r_urb(24)
     &         +.1*r_urb(27)                            
     &         +.08*r_urb(28)+.1*r_urb(32)+.06*r_urb(40)                                    
      d_urb(iethp)= 0.0                                                             
c
      p_urb(ihcho)= r_urb(7)+1.56*r_urb(8)+.57*r_urb(9)
     &         +r_urb(11)+r_urb(21)                             
     &         +.7*r_urb(23)+r_urb(29)+.5*r_urb(30)
     &         +r_urb(33)+.5*r_urb(34)                          
     &         +.7*r_urb(41)+.5*r_urb(42)                                               
      d_urb(ihcho)= 0.0                                                             
c
      p_urb(ich3oh)= .03*r_urb(9)+.04*r_urb(10)                                             
      d_urb(ich3oh)= 0.0                                                            
c
      p_urb(iald2)= .22*r_urb(8)+.47*r_urb(9)+1.03*r_urb(10)
     &         +r_urb(11)                              
     &         +1.77*r_urb(12)+.03*r_urb(23)+.3*r_urb(24)
     &         +.04*r_urb(25)                         
     &         +.3*r_urb(27)+.25*r_urb(28)+.5*r_urb(30)
     &         +.3*r_urb(32)                            
     &         +.5*r_urb(34)+.21*r_urb(40)+.5*r_urb(42)                                     
      d_urb(iald2)= 0.0                                                             
c
      p_urb(ihcooh)= .52*r_urb(7)+.22*r_urb(9)                                              
      d_urb(ihcooh)= 0.0                                                            
c
      p_urb(ircooh)= .09*r_urb(9)+.16*r_urb(10)                                             
      d_urb(ircooh)= 0.0                                                            
c
      p_urb(ic2o3)= r_urb(2)+r_urb(4)+r_urb(5)+r_urb(6)
     &         +.13*r_urb(9)+.19*r_urb(10)                          
     &         +r_urb(21)+r_urb(22)+.62*r_urb(23)
     &         +r_urb(29)+r_urb(33)                               
     &         +.7*r_urb(41)                                                        
      d_urb(ic2o3)= 0.0                                                             
c
      p_urb(ipar)= 1.1*r_urb(16)                                                        
      d_urb(ipar)= r_urb(1)+r_urb(44)                                                       
c
      p_urb(iaone)= .07*r_urb(10)+.23*r_urb(12)
     &         +.74*r_urb(24)+.74*r_urb(27)                         
     &         +.62*r_urb(28)+.74*r_urb(32)
     &         +.57*r_urb(40)+.15*r_urb(41)                         
      d_urb(iaone)= r_urb(2)+r_urb(3)                                                       
c
      p_urb(imgly)= .04*r_urb(9)+.07*r_urb(10)
     &         +.8*r_urb(16)+.2*r_urb(23)                            
     &         +.19*r_urb(25)+.15*r_urb(41)                                             
      d_urb(imgly)= r_urb(4)+r_urb(5)+r_urb(6)                                                  
c
      p_urb(ieth)= 0.0                                                              
      d_urb(ieth)= r_urb(7)+r_urb(8)                                                        
c
      p_urb(iolet)= 0.0                                                             
      d_urb(iolet)= r_urb(9)+r_urb(11)+r_urb(13)                                                
c
      p_urb(iolei)= 0.0                                                             
      d_urb(iolei)= r_urb(10)+r_urb(12)+r_urb(14)                                               
c
      p_urb(itol)= 0.0                                                              
      d_urb(itol)= r_urb(15)                                                            
c
      p_urb(ixyl)= 0.0                                                              
      d_urb(ixyl)= r_urb(16)                                                            
c
      p_urb(icres)= .12*r_urb(15)+.05*r_urb(16)                                             
      d_urb(icres)= r_urb(18)+r_urb(19)                                                     
c
      p_urb(ito2)= .8*r_urb(15)+.45*r_urb(16)                                               
      d_urb(ito2)= r_urb(17)                                                            
c
      p_urb(icro)= .4*r_urb(18)+r_urb(19)                                                   
      d_urb(icro)= r_urb(20)                                                            
c
      p_urb(iopen)= .95*r_urb(17)+.3*r_urb(18)                                              
      d_urb(iopen)= r_urb(21)+r_urb(22)+r_urb(23)                                               
c
      p_urb(ionit)= .05*r_urb(17)+r_urb(20)
     &         +.16*r_urb(28)+.5*r_urb(30)                              
     &         +.5*r_urb(34)+r_urb(38)+.5*r_urb(42)                                         
      d_urb(ionit)= r_urb(26)+r_urb(27)                                                     
c
      p_urb(irooh)= r_urb(36)+r_urb(37)                                                     
      d_urb(irooh)= r_urb(24)+r_urb(25)                                                     
c
      p_urb(iro2)= r_urb(1)+.03*r_urb(9)+.09*r_urb(10)
     &        +.77*r_urb(25)                                
      d_urb(iro2)= r_urb(28)+r_urb(32)+r_urb(36)+r_urb(40)                                          
c
      p_urb(iano2)= r_urb(3)+.11*r_urb(10)                                                  
      d_urb(iano2)= r_urb(29)+r_urb(33)+r_urb(37)+r_urb(41)                                         
c
      p_urb(inap)= r_urb(13)+r_urb(14)+r_urb(26)                                                
      d_urb(inap)= r_urb(30)+r_urb(34)+r_urb(38)+r_urb(42)                                          
c
      p_urb(ixo2)= r_urb(5)+r_urb(8)+r_urb(11)+r_urb(12)
     &        +.08*r_urb(15)                                  
     &        +.5*r_urb(16)+.6*r_urb(18)+r_urb(21)+.03*r_urb(23)                                
     &        +.4*r_urb(24)+.41*r_urb(27)+.34*r_urb(28)
     &        +.4*r_urb(32)+.24*r_urb(40)                                                        
      d_urb(ixo2)= r_urb(31)+r_urb(35)+r_urb(39)+r_urb(43)                                          
c
      p_urb(ixpar)= 1.06*r_urb(9)+2.26*r_urb(10)
     &         +r_urb(11)+2.23*r_urb(12)                           
     &         +1.98*r_urb(24)+.42*r_urb(25)+1.98*r_urb(27)                                 
     &         +1.68*r_urb(28)+r_urb(30)+1.98*r_urb(32)+r_urb(34)                               
     &         +1.25*r_urb(40)+r_urb(42)                                                
      d_urb(ixpar)= r_urb(44)          
c
      p_urb(isv1) = r_urb(15)*0.071+r_urb(16)*0.038
      d_urb(isv1) = 0.0
c
      p_urb(isv2) = r_urb(15)*0.138+r_urb(16)*0.167
      d_urb(isv2) = 0.0     
c      
      return
      end









c**************************************************************************
c subroutine PeroxyRateConstants: calculates parameterized thermal rate 
c                     constants for the alkylperoxy radical permutation 
c                     reactions for the entire mechanism.
c nomenclature:
c rk_param  = parameterized reaction rate constants (1/s)
c rk_perox  = individual permutation reaction rate constants (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c 
c author: Rahul A. Zaveri
c date  : June 1998
c
c-------------------------------------------------------------------------
      subroutine PeroxyRateConstants
      include 'chm.inc'
      include 'gas.inc'

      dimension sperox(nperox), rk_perox(nperox,nperox)
c
      sperox(jch3o2)  = cnn(kch3o2)
      sperox(jethp)   = cnn(kethp)
      sperox(jro2)    = cnn(kro2)
      sperox(jc2o3)   = cnn(kc2o3)
      sperox(jano2)   = cnn(kano2)
      sperox(jnap)    = cnn(knap)
      sperox(jisopp)  = cnn(kisopp)
      sperox(jisopn)  = cnn(kisopn)
      sperox(jisopo2) = cnn(kisopo2)
      sperox(jxo2)    = cnn(kxo2)

c
c initialize to zero
      do i = 1, nperox
      rk_param(i) = 0.0
      enddo

      do i = 1, nperox
      do j = 1, nperox
      rk_perox(i,j) = ARR(Aperox(i,j),Bperox(i,j))
      rk_param(i) = rk_param(i) + rk_perox(i,j)*sperox(j)
      enddo
      enddo
c
      return
      end








c**************************************************************************
      subroutine PhotoConstants_Fixed
      include 'chm.inc'
      include 'gas.inc'

      rk_photo(1)   = 0.0	! no2 + hv    --> no + o(3P)
      rk_photo(2)   = 0.0	! no3 + hv    --> .89no2 + .89o(3P) + .11no
      rk_photo(3)   = 0.0	! hono + hv   --> oh + no
      rk_photo(4)   = 0.0	! hno3 + hv   --> oh + no2
      rk_photo(5)   = 0.0	! hno4 + hv   --> ho2 + no2
      rk_photo(6)   = 0.0	! n2o5 + hv   --> no2 + no3
      rk_photo(7)   = 0.0	! o3 + hv     --> o(3P)
      rk_photo(8)   = 0.0	! o3 + hv     --> o(1D)
      rk_photo(9)   = 0.0	! h2o2 + hv   --> 2oh
      rk_photo(10)   = 0.0	! hcho + hv   --> 2ho2 + co
      rk_photo(11)  = 0.0	! hcho + hv   --> co
      rk_photo(12)  = 0.0	! ch3ooh + hv --> hcho + ho2 + oh
      rk_photo(13)  = 0.0	! ethooh + hv --> ald2 + ho2 + oh
      rk_photo(14)  = 0.0	! ald2 + hv   --> 
      rk_photo(15)  = 0.0	! aone + hv   --> 
      rk_photo(16)  = 0.0	! mgly + hv   --> 
      rk_photo(17)  = 0.0	! open + hv   --> 
      rk_photo(18)  = 0.0	! rooh + hv   -->
      rk_photo(19)  = 0.0	! onit + hv   --> 	
      rk_photo(20)  = 0.0	! isoprd + hv -->

      return
      end







c*************************************************************************
c subroutine PhotoConstants_Solar: calculates photochemical rate constants (1/s)
c
c input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
c        zalt_m (altitude above sea level in meters)
c
c-------------------------------------------------------------------------
 
      subroutine PhotoConstants_Solar(factcld)
      include 'chm.inc'
      include 'gas.inc'
      parameter (sza_cut = 89.00)		! cutoff solar zenith angle
      parameter (cos_sza_cut = 0.017452406)	! cos of sza_cut

!      factcld = 1.0
!        print*,factcld
      if(cos_sza .ge. cos_sza_cut)then	! daytime

         idaytime = 1	

         if(mphoto.eq.1)then
           call PhotoParam1
         elseif(mphoto.eq.2)then
           call PhotoParam2
         endif

c apply cloudiness correction factor
         do jphoto = 1, nphoto
           rk_photo(jphoto) = rk_photo(jphoto)*factcld
         enddo

      else				! nighttime

         idaytime = 0

         do jphoto = 1, nphoto
           rk_photo(jphoto) = 0.0
         enddo

      endif

      return
      end










c*************************************************************************
c subroutine PhotoParam1: calculates photochemical rate constants (1/s)
c
c input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
c        zalt_m (altitude above sea level in meters)
c
c-------------------------------------------------------------------------

      subroutine PhotoParam1
      include 'chm.inc'
      include 'gas.inc'

      z = min(1.e-4*zalt_m, 1.1)	! zalt in meters

c no2 + hv --> no + o3p
      alpha0=2.10223e-2-1.6e-3*(z-1.15)**2.
      alpha1=-12.2e-2+alpha0**.5
      beta0=1.258-.16*(z-2.15)**2.
      beta1=-1.2+beta0**.5
      rk_photo(jphoto_no2)    = alpha1*exp(beta1/cos_sza)

c no3 + hv --> .89 no2 + .89 o3p + .11 no
      rk_photo(jphoto_no3)    = 23.8*rk_photo(jphoto_no2)

c hono + hv --> oh + no
      rk_photo(jphoto_hono)   = 0.197*rk_photo(jphoto_no2)

c hno3 + hv --> oh + no2
      rk_photo(jphoto_hno3)   = 3.3e-5*rk_photo(jphoto_no2)

c hno4 + hv --> ho2 + no2
      rk_photo(jphoto_hno4)   = 5.9e-4*rk_photo(jphoto_no2)

c n2o5 + hv --> no2 + no3
      rk_photo(jphoto_n2o5)   = 0.0*rk_photo(jphoto_no2)

c o3 + hv --> o3p
      rk_photo(jphoto_o3a)    = 0.053*rk_photo(jphoto_no2)

c o3 + hv --> o1d
      alpha0=2.69924e-7-4.0e-8*(z-.375)**2.
      alpha7=-3.25e-4+sqrt(alpha0)
      beta0=4.173-.64*(z-2.0)**2.
      beta7=-3.2+sqrt(beta0)
      rk_photo(jphoto_o3b)    = alpha7*exp(beta7/cos_sza)

c h2o2 + hv --> 2 oh
      alpha0=2.540534e-9-4.e-11*(z-0.75)**2.
      alpha8=-3.e-5+sqrt(alpha0)
      beta0=.539284-.16*(z-1.75)**2.
      beta8=-1+sqrt(beta0)
      rk_photo(jphoto_h2o2)   = alpha8*exp(beta8/cos_sza)

c hcho + hv ---> 2 ho2 + co
      alpha0=7.1747e-9-1.6e-9*(z-.75)**2.
      alpha49=-4.e-5+sqrt(alpha0)
      beta0=.7631144-.16*(z-2.)**2.
      beta49=-1.2+sqrt(beta0)
      rk_photo(jphoto_hchoa)  = alpha49*exp(beta49/cos_sza)

c hcho + hv ---> co
      alpha0=1.693813e-7-4.e-8*(z-.875)**2.
      alpha50=-2.5e-4+sqrt(alpha0)
      beta0=.7631144-.16*(z-1.875)**2.
      beta50=-1.1+sqrt(beta0)
      rk_photo(jphoto_hchob)  = alpha50*exp(beta50/cos_sza)

      rk_photo(jphoto_ch3ooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ethooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ald2)   = 4.6e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_aone)   = 7.8e-5*rk_photo(jphoto_no2)
      rk_photo(jphoto_mgly)   = 9.64  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_open)   = 9.04  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_rooh)   = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_onit)   = 1.0e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_isoprd) = .025  *rk_photo(jphoto_hchob)

      return
      end









c*************************************************************************
c subroutine PhotoParam2: calculates photochemical rate constants (1/s)
c
c input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
c        zalt_m (altitude above sea level in meters)
c
c-------------------------------------------------------------------------

      subroutine PhotoParam2
      include 'chm.inc'
      include 'gas.inc'

      real kr(3,nphoto)	! kr(level, species #):

      cz = cos_sza


c surface photolysis rates:
c NO2
	KR(1,jphoto_no2)=
     +  -1.0184E-3+cz*1.8542E-2-cz**2*9.5368E-3+cz**3*1.8165E-3
c NO3       
	KR(1,jphoto_no3)=
     +  4.3945E-3+0.556446*cz-0.71996*cz**2+0.34253*cz**3
c O3      
        KR(1,jphoto_o3b)=(-3.2168E-8+cz*2.4588E-6-cz**2*2.6274E-5+
     + 	1.2005E-4*cz**3-7.6216E-5*cz**4+1.4628E-5*cz**5)
c HONO
	KR(1,jphoto_hono)=
     +  -1.7863E-4+3.2272E-3*cz-8.5989E-4*cz**2-1.8987E-4*cz**3
c HNO3       
	KR(1,jphoto_hno3)=
     +  1.9592E-8-2.8147E-7*cz+1.3533E-6*cz**2-4.2010E-7*cz**3
c H2O2     
	KR(1,jphoto_h2o2)=
     +  -2.1672E-7+4.0070E-6*cz+1.1189E-5*cz**2-6.4306E-6*cz**3
c HNO4      
	KR(1,jphoto_hno4)=
     +  2.1392E-8-2.0854E-7*cz+1.6131E-5*cz**2-7.2297E-6*cz**3
c HCHOB      
	KR(1,jphoto_hchoa)=
     +  -5.4160E-8+9.4694E-7*cz+6.4697E-5*cz**2-3.2594E-5*cz**3
c HCHOA      
	KR(1,jphoto_hchob)=
     +  -2.6855E-6+4.9102E-5*cz+5.7086E-5*cz**2-4.3040E-5*cz**3


c photolysis rates at 4 km:
c NO2	
  	KR(2,jphoto_no2)=
     +  -1.3136E-3+2.4948E-2*cz-1.9513E-2*cz**2+6.611E-3*cz**3
c NO3	       
	KR(2,jphoto_no3)=
     +  1.59E-2+0.54202*cz-0.72079*cz**2+0.34898*cz**3 
c O3
        KR(2,jphoto_o3b)=
     +    (1.6295E-7+cz*4.9940E-7-cz**2*2.9055E-5+1.8187E-4*cz**3
     +              -1.5627E-4*cz**4+4.5975E-5*cz**5)
c HONO       
	KR(2,jphoto_hono)=
     +  -2.6339E-4+4.6997E-3*cz-2.9408E-3*cz**2+7.4996E-4*cz**3
c HNO3       
	KR(2,jphoto_hno3)=
     +  2.2106E-8-3.4422E-7*cz+1.8449E-6*cz**2-6.7994E-7*cz**3
c H2O2       
	KR(2,jphoto_h2o2)=
     +  -4.73E-7+7.4881E-6*cz+9.7183E-6*cz**2-6.4955E-6*cz**3 
c HNO4       
	KR(2,jphoto_hno4)=
     +  -1.0672E-7+1.1698E-6*cz+1.9044E-5*cz**2-9.4072E-6*cz**3
c HCHOB       
	KR(2,jphoto_hchoa)=
     +  -7.4493E-7+8.7149E-6*cz+7.1885E-5*cz**2-3.9526E-5*cz**3 
c HCHOA       
	KR(2,jphoto_hchob)=
     +  -5.1681E-6+8.4398E-5*cz+2.6478E-5*cz**2-3.4452E-5*cz**3 


c photolysis rates at 8 km:
c NO2
       KR(3,jphoto_no2)=
     + -1.3748E-3+2.9757E-2*cz-2.8355E-2*cz**2+1.1168E-2*cz**3
c NO3
       KR(3,jphoto_no3)=
     + 2.80132E-2+0.51381*cz-0.68839*cz**2+0.33448*cz**3
c O3
       KR(3,jphoto_o3b)=
     + (1.6295E-7+cz*4.9940E-7-cz**2*2.9055E-5+1.8187E-4*cz**3
     +              -1.5627E-4*cz**4+4.5975E-5*cz**5)
c HONO
       KR(3,jphoto_hono)=
     + -3.1944E-4+6.0983E-3*cz-5.2694E-3*cz**2+1.9111E-3*cz**3
c HNO3
       KR(3,jphoto_hno3)=
     + 1.9176E-8-3.4083E-7*cz+2.1560E-6*cz**2-8.7941E-7*cz**3
c H2O2
       KR(3,jphoto_h2o2)=
     + -7.6642E-7+1.1717E-5*cz+5.3611E-6*cz**2-4.9358E-6*cz**3
c HNO4
       KR(3,jphoto_hno4)=
     + -3.2131E-7+3.6898E-6*cz+1.8481E-5*cz**2-9.826E-6*cz**3
c HCHOB
       KR(3,jphoto_hchoa)=
     + -1.7563E-6+2.0714E-5*cz+6.5668E-5*cz**2-3.9386E-5*cz**3
c HCHOA
       KR(3,jphoto_hchob)=
     + -7.9124E-6+1.258E-4*cz-2.8767E-5*cz**2-1.0505E-5*cz**3


c force all the kr to be non-negative
	do jdum = 1, nphoto
	do idum = 1, 3
	    kr(idum,jdum) = max( 0., kr(idum,jdum) )
	end do
	end do


c above 8km, use values at 8km
	if (zalt_m .ge. 8000) then
	 alpha = 1
	 k = 3
	 km1 = 3
c linear interpolation 
	else if (zalt_m .lt. 8000. .and.
     +           zalt_m .ge. 4000.) then
	 alpha = (8000. - zalt_m)/4000.
	 k = 3
	 km1 = 2
	else if (zalt_m .lt. 4000) then
	 alpha = (4000. - zalt_m)/4000
	 k = 2
	 km1 = 1
	end if

        a = alpha
        a1= 1. - alpha

      rk_photo(jphoto_no2)   = a1*kr(k,jphoto_no2)   +
     &                         a*kr(km1,jphoto_no2)
      rk_photo(jphoto_no3)   = a1*kr(k,jphoto_no3)   +
     &                         a*kr(km1,jphoto_no3)
      rk_photo(jphoto_o3a)   = 0.053*rk_photo(jphoto_no2)

      rk_photo(jphoto_o3b)   = a1*kr(k,jphoto_o3b)   +
     &                         a*kr(km1,jphoto_o3b)
      rk_photo(jphoto_hono)  = a1*kr(k,jphoto_hono)  +
     &                         a*kr(km1,jphoto_hono)
      rk_photo(jphoto_hno3)  = a1*kr(k,jphoto_hno3)  +
     &                         a*kr(km1,jphoto_hno3)
      rk_photo(jphoto_h2o2)  = a1*kr(k,jphoto_h2o2)  +
     &                         a*kr(km1,jphoto_h2o2)
      rk_photo(jphoto_hno4)  = a1*kr(k,jphoto_hno4)  +
     &                         a*kr(km1,jphoto_hno4)
      rk_photo(jphoto_hchoa) = a1*kr(k,jphoto_hchoa) +
     &                         a*kr(km1,jphoto_hchoa)
      rk_photo(jphoto_hchob) = a1*kr(k,jphoto_hchob) +
     &                         a*kr(km1,jphoto_hchob)

      rk_photo(jphoto_n2o5)   = 0.0   *rk_photo(jphoto_no2)
      rk_photo(jphoto_ch3ooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ethooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ald2)   = 4.6e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_aone)   = 7.8e-5*rk_photo(jphoto_no2)
      rk_photo(jphoto_mgly)   = 9.64  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_open)   = 9.04  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_rooh)   = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_onit)   = 1.0e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_isoprd) = .025  *rk_photo(jphoto_hchob)

      return
      end





             











c***********************************************************************
c subroutine SetGas_Bio: sets up gas-phase species indices for 
c the selected mechanism.
c
c author: Rahul A. Zaveri
c date  : february 1996
c-------------------------------------------------------------------------
      subroutine SetGas_Bio(ilast)
      include 'chm.inc'
      include 'gas.inc'

      iisop	= ilast + 1
      iisoprd	= ilast + 2
      iisopp	= ilast + 3
      iisopn	= ilast + 4
      iisopo2	= ilast + 5
      iterp     = ilast + 6
      isv3      = ilast + 7
      isv4      = ilast + 8
      isv5      = ilast + 9
      isv6      = ilast + 10
      
      ilast	= isv6

      return
      end









c***********************************************************************
c subroutine SetGas_Com: sets up gas-phase species indices for 
c the selected mechanism.
c
c author: Rahul A. Zaveri
c date  : february 1996
c-------------------------------------------------------------------------
      subroutine SetGas_Com(ilast)
      include 'chm.inc'
      include 'gas.inc'

      ih2so4	= 1
      ihno3	= 2
      ihcl	= 3
      inh3	= 4
      ino	= 5
      ino2	= 6
      ino3	= 7
      in2o5	= 8
      ihono	= 9
      ihno4	= 10
      io3	= 11
      io1d	= 12
      io3p	= 13
      ioh	= 14
      iho2	= 15
      ih2o2	= 16
      ico	= 17
      iso2	= 18
      ich4	= 19
      ic2h6	= 20
      ich3o2	= 21
      iethp	= 22
      ihcho	= 23
      ich3oh	= 24
      ianol	= 25
      ich3ooh	= 26
      iethooh	= 27
      iald2	= 28
      ihcooh	= 29
      ircooh	= 30
      ic2o3	= 31
      ipan	= 32
      ilast	= ipan
 
      return
      end













c***********************************************************************
c subroutine SetGasIndices: sets up gas-phase species indices for 
c the selected mechanism.
c
c input: iregime    = 1     : com
c                   = 2     : com + urb
c                   = 3     : com + urb + bio
c                   = 4     : com + mar
c                   = 5     : com + urb + mar
c                   = 6     : com + urb + bio + mar
c
c author: Rahul A. Zaveri
c date  : february 1996
c-------------------------------------------------------------------------
      subroutine SetGasIndices
      include 'chm.inc'
      include 'gas.inc'
c
      ilast = 0

      goto (1,2,3,4,5,6), iregime

1     call setgas_com(ilast)
      return
c
c
2     call setgas_com(ilast)
      call setgas_urb(ilast)
      return
c
c
3     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      return
c
c
4     call setgas_com(ilast)
      call setgas_mar(ilast)
      return
c
c 
5     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_mar(ilast)
      return
c
c
6     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      call setgas_mar(ilast)
      return
c
      end










c***********************************************************************
c subroutine SetGas_Mar: sets up gas-phase species indices for 
c the selected mechanism.
c
c author: Rahul A. Zaveri
c date  : february 1996
c-------------------------------------------------------------------------
      subroutine SetGas_Mar(ilast)
      include 'chm.inc'
      include 'gas.inc'

      idms         = ilast + 1
      imsa         = ilast + 2
      idmso        = ilast + 3
      idmso2       = ilast + 4
      ich3so2h     = ilast + 5
      ich3sch2oo   = ilast + 6
      ich3so2      = ilast + 7
      ich3so3      = ilast + 8
      ich3so2oo    = ilast + 9
      ich3so2ch2oo = ilast + 10
      isulfhox     = ilast + 11

      ilast	   = isulfhox

      return
      end












c***********************************************************************
c subroutine SetGas_Urb: sets up gas-phase species indices for 
c the selected mechanism.
c
c author: Rahul A. Zaveri
c date  : february 1996
c-------------------------------------------------------------------------
      subroutine SetGas_Urb(ilast)
      include 'chm.inc'
      include 'gas.inc'

      ipar	= ilast + 1
      iaone	= ilast + 2
      imgly	= ilast + 3
      ieth	= ilast + 4
      iolet	= ilast + 5
      iolei	= ilast + 6
      itol	= ilast + 7
      ixyl	= ilast + 8
      icres	= ilast + 9
      ito2	= ilast + 10
      icro	= ilast + 11
      iopen	= ilast + 12
      ionit	= ilast + 13
      irooh	= ilast + 14
      iro2	= ilast + 15
      iano2	= ilast + 16
      inap	= ilast + 17
      ixo2	= ilast + 18
      ixpar	= ilast + 19
      isv1      = ilast + 20
      isv2      = ilast + 21
      
      ilast	= isv2

      return
      end
      subroutine SetRunParameters(trunm, dtmin)
      include 'chm.inc'

      tbeg_sec = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss
      trun_sec = ((trun_dd*24.+trun_hh)*60.+trun_mm)*60.+trun_ss

      dt_sec   = 60.*dt_min	! time step [seconds]
      ncstep = INT(trun_sec/dt_sec + 0.01)	! number of time steps

      tcur_sec = tbeg_sec	! initialize current time to tbeg_sec
      tcur_min = tcur_sec/60.
      tcur_hrs = tcur_min/60.
      it       = 0		! initialize step counter to zero

c
c convert rlon and rlat to radians
      rlon = rlon*deg2rad
      rlat = rlat*deg2rad       
c
      if(msolar .eq. 2)then
        write(6,*)' You have selected fixed photolysis rates'
        write(6,*)' for the entire simulation (msolar = 2).'
        write(6,*)' Set the rate constants in file:'
        write(6,*)' photoconstants_fixed.f'
      endif
c
      return
      end















c*******************************************************************
c      SUBROUTINE SOLAR - calculates cosine of zenith angle
c                         for use in photochemical rate coefficient
c                         calculations.
c
c      Nomenclature:
c  
c      tmid_sec = time in seconds from Greenich Noon March 21
c
c      cos_za = cosine of zenith angle
c
c      rlon   = local longitude, W hemis. values are negative (radians)
c                   
c      rlat   = local latitude, S hemis. values are negative (radians)
c
c*******************************************************************
 
       subroutine SolarZenithAngle
       include 'chm.inc'
c 
       tlocal=tmid_sec                                                  
       tdec=0.4092797*sin(1.992385E-7*tlocal)                           
       sidec=sin(tdec)                                                  
       codec=cos(tdec)                                                  
       tloc=7.272205E-5*tlocal                                           
       thou=cos(rlon+tloc)                                      
       cos_sza=sin(rlat)*sidec+cos(rlat)*codec*thou         
       return
       end





c------------------------------------------------------------
      subroutine UpdateTime

      include 'chm.inc'

      tsav_sec = tcur_sec
      told_sec = tcur_sec
      tcur_sec = tcur_sec + dt_sec
      tcur_min = tcur_sec/60.
      tcur_hrs = tcur_min/60.
      tmid_sec = told_sec + 0.5*dt_sec

      return
      end




c*******************************************************************
        block data
c
c   place various constant or initial values into common
c
      include 'chm.inc'
      include 'gas.inc'

      common /eh0001/ mesflg, lunit
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu


c--------------------------------------------------
c     define fundamental constants...
      data pi		/3.141592654/
      data avogad	/6.02217e+23/
      data deg2rad	/0.017453293/

c---------------------------------
c define species indices
c
c species in inorganic chemistry
      data kh2so4       /  1/
      data khno3        /  2/
      data khcl         /  3/
      data knh3         /  4/
      data kno          /  5/
      data kno2         /  6/
      data kno3         /  7/
      data kn2o5        /  8/
      data khono        /  9/
      data khno4        / 10/
      data ko3          / 11/
      data ko1d         / 12/
      data ko3p         / 13/
      data koh          / 14/
      data kho2         / 15/
      data kh2o2        / 16/
      data kco          / 17/
      data kso2         / 18/
c
c species in methane, ethane, formaldehyde chemistry
      data kch4         / 19/
      data kc2h6        / 20/
      data kch3o2       / 21/
      data kethp        / 22/
      data khcho        / 23/
      data kch3oh       / 24/
      data kanol	/ 25/
      data kch3ooh      / 26/
      data kethooh	/ 27/
      data kald2        / 28/
      data khcooh	/ 29/
      data krcooh	/ 30/
      data kc2o3	/ 31/
      data kpan		/ 32/
c 
c species in hc1 mechanism. initialize indices to zero
      data kpar         / 33/
      data kaone        / 34/
      data kmgly        / 35/
      data keth         / 36/
      data kolet        / 37/
      data kolei        / 38/
      data ktol         / 39/
      data kxyl         / 40/
      data kcres        / 41/
      data kto2         / 42/
      data kcro         / 43/
      data kopen        / 44/
      data konit        / 45/
      data krooh	/ 46/
      data kro2         / 47/
      data kano2	/ 48/
      data knap		/ 49/
      data kxo2		/ 50/
      data kxpar	/ 51/
c
c species in hc2 mechanism. initialize indices to zero
      data kisop	/ 52/
      data kisoprd	/ 53/
      data kisopp	/ 54/
      data kisopn	/ 55/
      data kisopo2	/ 56/

c species in dms mechanism. initialize indices to zero
      data kdms         / 57/
      data kmsa         / 58/
      data kdmso        / 59/
      data kdmso2       / 60/
      data kch3so2h     / 61/
      data kch3sch2oo   / 62/
      data kch3so2      / 63/
      data kch3so3      / 64/
      data kch3so2oo    / 65/
      data kch3so2ch2oo / 66/
      data ksulfhox     / 67/
c species in soa formation 2 species in urban, 5 species in biogenic
      data kterp        / 68/ 
      data ksv1         / 69/
      data ksv2         / 70/
      data ksv3         / 71/
      data ksv4         / 72/
      data ksv5         / 73/
      data ksv6         / 74/
c
c regime-dependent chemistry definitions
c
      data iregime	/  1/
c
c     GAS
c
c species in common chemistry
      data ih2so4       /  1/
      data ihno3        /  2/
      data ihcl         /  3/
      data inh3         /  4/
      data ino          /  5/
      data ino2         /  6/
      data ino3         /  7/
      data in2o5        /  8/
      data ihono        /  9/
      data ihno4        / 10/
      data io3          / 11/
      data io1d         / 12/
      data io3p         / 13/
      data ioh          / 14/
      data iho2         / 15/
      data ih2o2        / 16/
      data ico          / 17/
      data iso2         / 18/
c
c species in methane, ethane, formaldehyde chemistry
      data ich4         / 19/
      data ic2h6        / 20/
      data ich3o2       / 21/
      data iethp        / 22/
      data ihcho        / 23/
      data ich3oh       / 24/
      data ianol	/ 25/
      data ich3ooh      / 26/
      data iethooh	/ 27/
      data iald2        / 28/
      data ihcooh	/ 29/
      data ircooh	/ 30/
      data ic2o3	/ 31/
      data ipan		/ 32/
c 
c species in hc1 mechanism. initialize indices to zero
      data ipar         / 33/
      data iaone        / 34/
      data imgly        / 35/
      data ieth         / 36/
      data iolet        / 37/
      data iolei        / 38/
      data itol         / 39/
      data ixyl         / 40/
      data icres        / 41/
      data ito2         / 42/
      data icro         / 43/
      data iopen        / 44/
      data ionit        / 45/
      data irooh	/ 46/
      data iro2         / 47/
      data iano2	/ 48/
      data inap		/ 49/
      data ixo2		/ 50/
      data ixpar	/ 51/
c
c species in hc2 mechanism. initialize indices to zero
      data iisop	/ 52/
      data iisoprd	/ 53/
      data iisopp	/ 54/
      data iisopn	/ 55/
      data iisopo2	/ 56/

c species in dms mechanism. initialize indices to zero
      data idms         / 57/
      data imsa         / 58/
      data idmso        / 59/
      data idmso2       / 60/
      data ich3so2h     / 61/
      data ich3sch2oo   / 62/
      data ich3so2      / 63/
      data ich3so3      / 64/
      data ich3so2oo    / 65/
      data ich3so2ch2oo / 66/
      data isulfhox     / 67/
c intial values for soa formation
      data iterp        / 68/
      data isv1         / 69/
      data isv2         / 70/
      data isv3         / 71/
      data isv4         / 72/
      data isv5         / 73/
      data isv6         / 74/
      
c alkylperoxy radical indices for parameterized permutation reactions
      data jch3o2	/  1/
      data jethp	/  2/
      data jro2		/  3/
      data jc2o3	/  4/
      data jano2	/  5/
      data jnap		/  6/
      data jisopp	/  7/
      data jisopn	/  8/
      data jisopo2	/  9/
      data jxo2		/ 10/

c photolyzing species indices
      data jphoto_no2	/  1/
      data jphoto_no3	/  2/
      data jphoto_hono	/  3/
      data jphoto_hno3	/  4/
      data jphoto_hno4	/  5/
      data jphoto_n2o5  /  6/
      data jphoto_o3a	/  7/
      data jphoto_o3b	/  8/
      data jphoto_h2o2	/  9/
      data jphoto_hchoa	/ 10/
      data jphoto_hchob	/ 11/
      data jphoto_ch3ooh/ 12/
      data jphoto_ethooh/ 13/
      data jphoto_ald2	/ 14/
      data jphoto_aone	/ 15/
      data jphoto_mgly	/ 16/
      data jphoto_open	/ 17/
      data jphoto_rooh	/ 18/
      data jphoto_onit	/ 19/
      data jphoto_isoprd/ 20/

c
c lsodes parameters
      data illin/0/, ntrep/0/
      data mesflg/1/, lunit/6/

      end
      subroutine DoMassBalance
      implicit real(a-h,o-z), integer(i-n) 
      include 'chm.inc'

c
c=================================================================
c Mass Balances
c=================================================================
c     initialize...
      totS = 0.0
      totN = 0.0
      totCl= 0.0
c
c__________________________________________________________________
c Sulfur Balance
c
      do i=kdms,ksulfhox
      totS = totS + cnn(i)
      enddo
      totS = totS + cnn(kso2) + cnn(kh2so4)
 
c__________________________________________________________________
c Nitrogen Balance
c
      do i=kno,khno4
      totN = totN + cnn(i)
      enddo
      totN = totN + cnn(kn2o5)+cnn(kpan)+cnn(konit)+cnn(knap)+
     &              cnn(kisopn)
c
c__________________________________________________________________
c Chlorine Balance
c
      totCl = cnn(khcl)
c
c__________________________________________________________________
c
c  Initial total elements (N,S,Cl)
      tNid = totN ! juanxiong he
      tSid = totS ! juanxiong he
      tClid= totCl ! juanxiong he
      if(it.eq.0) then

          if(totN.gt.1.E-30)then
          tNi = totN
          tNid= totN
          else
          tNi = totN
          tNid= 1.0
          endif

          if(totS.gt.1.E-30)then
          tSi = totS
          tSid= totS
          else
          tSi = totS
          tSid= 1.0
          endif

          if(totCl.gt.1.E-30)then
          tCli = totCl
          tClid= totCl
          else
          tCli = totCli
          tClid= 1.0
          endif

      end if

c
c  Calculate percent deviation in elemental mass balance
      DN = 100.*(totN-tNi)/tNid
      DS = 100.*(totS-tSi)/tSid
      DCl= 100.*(totCl-tCli)/tClid
 
      return
      end
      Function ARR(AA,BB)
      include 'chm.inc'	! for te
      ARR = AA*exp(BB/te)
      return
      end
	integer function nbllen( str )
c
c   returns the position of the last non-blank character in str
c
	character*(*) str

	j = len(str)

	if (j .gt. 0) then
1000	    if (str(j:j) .eq. ' ') then
		j = j - 1
		if (j .gt. 0) goto 1000
	    end if
	end if
	nbllen = max0( j, 0 )

	return
	end
      Function Troe(cair_mlc,te,rk0,rnn,rki,rmm)
      rk0 = rk0*cair_mlc*(te/300.)**(-rnn)
      rki = rki*(te/300.)**(-rmm)
      expo= 1./(1. + (ALOG10(rk0/rki))**2)
      troe  = (rk0*rki/(rk0+rki))*.6**expo
      return
      end





c**********************************************************************
c function WaterVapor
c purpose: calculates concentration of h2o using the method given in 
c          Seinfeld's book, pg. 181
c----------------------------------------------------------------------

      function WaterVapor(rh, cair_mlc, te, pr_atm)

      t_steam = 373.15			! steam temperature  [K]
      pr_std   = 1.0			! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)  			! [atm]
      WaterVapor = RH*(pr_h2o/pr_atm)*cair_mlc/100.	! [molec/cc]

      return
      end





c**********************************************************************
      subroutine GasChemistry(t_in, t_out,factcld)
      include 'chm.inc'
      include 'gas.inc'
      dimension stot(nmax)		! local species array

      call SelectGasRegime(ntot)	! selects iregime and calculates ntot

      call PeroxyRateConstants

      call GasRateConstants(factcld)

      call SetGasIndices		! set gas indices for selected iregime 

      call MapGasSpecies(stot,0)	! map cnn into stot for selected iregime

      call ODEsolver(ntot,stot,t_in,t_out)

      call MapGasSpecies(stot,1)	! map stot back into cnn

      return
      end





c**************************************************************************
c subroutine GasRateConstants_Bio: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_bio    = reaction rate constants for hc2 mechanism    (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants_Bio
      include 'chm.inc'
      include 'gas.inc'

      rk_bio(1)  = ARR(2.6e-11, 409.)
      rk_bio(2)  = ARR(1.2e-14, -2013.)
      rk_bio(3)  = ARR(3.0e-12, -446.)
      rk_bio(4)  = rk_photo(jphoto_isoprd)
      rk_bio(5)  = 3.3e-11
      rk_bio(6)  = 7.0e-18
      rk_bio(7)  = 1.0e-15
      rk_bio(8)  = 4.0e-12
      rk_bio(9)  = 4.0e-12
      rk_bio(10) = 4.0e-12
      rk_bio(11) = ARR(1.7e-13, 1300.)
      rk_bio(12) = ARR(1.7e-13, 1300.)
      rk_bio(13) = ARR(1.7e-13, 1300.)
      rk_bio(14) = rk_param(jisopp)
      rk_bio(15) = rk_param(jisopn)
      rk_bio(16) = rk_param(jisopo2)
      rk_bio(17) = 3.563e-11
      rk_bio(18) = ARR(6.712e-11,-449.)
      rk_bio(19) = ARR(7.44e-17,821.)
      rk_bio(20) = ARR(6.642e-12,-175.)
      
      return
      end





c**************************************************************************
c subroutine GasRateConstants_Com: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_com    = reaction rate constants for common mechanism (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c iregime = selected mechanism for the current chemical regime (1-6) 
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants_Com
      include 'chm.inc'
      include 'gas.inc'

      rk_com(1) = rk_photo(jphoto_no2)                                                         
      rk_com(2) = rk_photo(jphoto_no3)                                                         
      rk_com(3) = rk_photo(jphoto_hono)                                                        
      rk_com(4) = rk_photo(jphoto_hno3)                                                        
      rk_com(5) = rk_photo(jphoto_hno4)                                                       
      rk_com(6) = rk_photo(jphoto_n2o5)                                                        
      rk_com(7) = rk_photo(jphoto_o3a)                                                         
      rk_com(8) = rk_photo(jphoto_o3b)                                                         
      rk_com(9) = rk_photo(jphoto_h2o2)                                                        
      rk_com(10) = ARR(3.2e-11, 70.)                                                
      rk_com(11) = ARR(1.8e-11, 110.)                                               
      rk_com(12) = 2.2e-10                                                          
      rk_com(13) = cair_mlc*6.e-34*(te/300.)**(-2.3)                                                        
      rk_com(14) = ARR(8.0e-12, -2060.)                                             
      rk_com(15) = ARR(6.5e-12, -120.)                 

      rk0 = 9.0e-32
      rnn = 2.0
      rki = 2.2e-11
      rmm = 0.0                             
      rk_com(16) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)   

      rk0 = 9.0e-32
      rnn = 1.5
      rki = 3.0e-11
      rmm = 0.0                                 
      rk_com(17) = Troe(cair_mlc,te,rk0,rnn,rki,rmm) 
                                   
      rk_com(18) = ARR(2.0e-12, -1400.)                                             
      rk_com(19) = ARR(1.2e-13, -2450.)                                             
      rk_com(20) = ARR(1.6e-12, -940.)                                              
      rk_com(21) = ARR(1.1e-14, -500.)                                              
      rk_com(22) = ARR(5.5e-12, -2000.)  
                                           
      rk0 = 7.0e-31
      rnn = 2.6
      rki = 3.6e-11
      rmm = 0.1
      rk_com(23) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)   

      rk0 = 2.5e-30
      rnn = 4.4
      rki = 1.6e-11
      rmm = 1.7
      rk_com(24) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)                                    
      rk_com(25) = 2.2e-11                                                          
      rk_com(26) = ARR(1.8e-11, -390.)

             rko = 7.2e-15 * exp(785./te)
             rk2 = 4.1e-16 * exp(1440./te)
             rk3 = 1.9e-33 * exp(725./te)*cair_mlc
      RK_OH_HNO3 = rko + rk3/(1.+rk3/rk2)                                          
      rk_com(27) = RK_OH_HNO3                                                       
      rk_com(28) = ARR(1.3e-12, 380.)                                               
      rk_com(29) = ARR(4.8e-11, 250.)                                               
      rk_com(30) = ARR(2.9e-12, -160.)

      RK_2HO2    = 2.3e-13 * exp(600./te) + 	! ho2 + ho2 --> h2o2
     &             1.7e-33 * exp(1000./te)*cair_mlc                           
      rk_com(31) = RK_2HO2

      RK_2HO2_H2O= RK_2HO2*1.4e-21*exp(2200./te)! ho2 + ho2 + h2o --> h2o2                                                          
      rk_com(32) = RK_2HO2_H2O
                                      
      rk_com(33) = ARR(3.5e-12, 250.)     

      rk0 = 1.8e-31
      rnn = 3.2
      rki = 4.7e-12
      rmm = 1.4                                          
      rk_com(34) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)    
                                
      rk_com(35) = 5.0e-16                                                          
      rk_com(36) = rk_com(34)*ARR(4.8e26, -10900.)                                  
      rk_com(37) = ARR(1.5e-11, 170.)                                               
      rk_com(38) = ARR(4.5e-14, -1260.)     

      rk0 = 2.2e-30
      rnn = 3.9
      rki = 1.5e-12
      rmm = 0.7                  
      rk_com(39) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)         
                           
      rk_com(40) = ARR(8.5e-13, -2450.)                                             
      rk_com(41) = 3.5e-12                                                          
      rk_com(42) = 2.0e-21                                                          
      rk_com(43) = rk_com(39)*ARR(3.7e26, -11000.)     

      RK_CO_OH   = 1.5e-13 * (1.+8.18e-23*te*cair_mlc) ! co + oh --> ho2                             
      rk_com(44) = RK_CO_OH                   
                                      
      rk0 = 3.0e-31
      rnn = 3.3
      rki = 1.5e-12
      rmm = 0.0
      rk_com(45) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)  
                                  
      rk_com(46) = te**.667*ARR(2.8e-14, -1575.)                                    
      rk_com(47) = te**2*ARR(1.5e-17, -492.)                                        
      rk_com(48) = ARR(6.7e-12, -600.)                                              
      rk_com(49) = rk_photo(jphoto_hchoa)                                                      
      rk_com(50) = rk_photo(jphoto_hchob)                                                      
      rk_com(51) = 1.0e-11                                                          
      rk_com(52) = ARR(3.4e-13, -1900.)                                             
      rk_com(53) = rk_photo(jphoto_ch3ooh)                                                    
      rk_com(54) = rk_photo(jphoto_ethooh)                                                    
      rk_com(55) = ARR(3.8e-12, 200.)                                               
      rk_com(56) = ARR(3.8e-12, 200.)                                               
      rk_com(57) = ARR(3.0e-12, 280.)                                               
      rk_com(58) = ARR(2.6e-12, 365.)                                               
      rk_com(59) = 1.1e-12                                                          
      rk_com(60) = 2.5e-12                                                          
      rk_com(61) = ARR(3.8e-13, 800.)                                               
      rk_com(62) = ARR(7.5e-13, 700.)                                               
      rk_com(63) = rk_param(jch3o2)                                                 
      rk_com(64) = rk_param(jethp)                                                  
      rk_com(65) = ARR(7.0e-12, -235.)                                              
      rk_com(66) = rk_photo(jphoto_ald2)                                               
      rk_com(67) = ARR(5.6e-12, 270.)                                               
      rk_com(68) = ARR(1.4e-12, -1900.)       

      rk0 = 9.7e-29
      rnn = 5.6
      rki = 9.3e-12
      rmm = 1.5                                      
      rk_com(69) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)          
                          
      rk_com(70) = rk_com(69)*ARR(1.1e28, -14000.)                                  
      rk_com(71) = ARR(5.3e-12, 360.)                                               
      rk_com(72) = 4.0e-12                                                          
      rk_com(73) = ARR(4.5e-13, 1000.)                                              
      rk_com(74) = rk_param(jc2o3)      

      return
      end





c**************************************************************************
c subroutine GasRateConstants: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_com    = reaction rate constants for common mechanism (molec-cc-s)
c rk_urb    = reaction rate constants for hc1 mechanism    (molec-cc-s)
c rk_bio    = reaction rate constants for hc2 mechanism    (molec-cc-s)
c rk_mar    = reaction rate constants for marine mechanism (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c iregime = selected mechanism for the current chemical regime (1-6) 
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants(factcld)
      include 'chm.inc'
      include 'gas.inc'

      if(msolar.eq.1)then
       call SolarZenithAngle		! calculates cos_sza
       call PhotoConstants_Solar(factcld)	! natural diurnal variation
      elseif(msolar.eq.2)then
       call PhotoConstants_Fixed	! artificial as in a smog chamber
      endif

      call gasrateconstants_het

      goto (1,2,3,4,5,6), iregime

1     call gasrateconstants_com
      return

2     call gasrateconstants_com
      call gasrateconstants_urb
      return

3     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_bio
      return

4     call gasrateconstants_com
      call gasrateconstants_mar
      return
 
5     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_mar
      return

6     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_bio
      call gasrateconstants_mar
      return
c
      end




c**************************************************************************
c subroutine GasRateConstants_Het: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_het    = reaction rate constants for heterogeneous rxns (1/s)
c 
c author: Rahul A. Zaveri
c date  : August 2000
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants_Het
      include 'chm.inc'
      include 'gas.inc'

      real Kn_o3, kn_n2o5, lambda_o3, lambda_n2o5, Dp(15)	! Dp in cm

      data dp /.12e-4, .14e-4, .17e-4, .20e-4, .25e-4, .31e-4,
     &	       .40e-4, .54e-4, .76e-4, .95e-4,
     &         1.20e-4, 1.32e-4, 1.75e-4, 2.25e-4, 2.75e-4/



      alpha_o3  = 1.1e-3		! accommodation coefficient [-]
      alpha_n2o5= 0.02*RH/100
      
      Dg_o3     = 0.1				! diffusivity [cm^2/s]
      Dg_n2o5   = 0.1
      
      c_o3      = 1.455e4 * sqrt(te/48.)	! avg. molec speed [cm/s]
      c_n2o5    = 1.455e4 * sqrt(te/108.)
      c_no3     = 1.455e4 * sqrt(te/62.)
      
      lambda_o3  = 3.*Dg_o3/c_o3			! mean free path [cm]
      lambda_n2o5= 3.*Dg_n2o5/c_n2o5

      sum_o3 = 0.0
      sum_n2o5= 0.0
      do k = 1, 15
        Kn_o3  = 2.*lambda_o3/Dp(k)		! knudsen number [-]
        Kn_n2o5= 2.*lambda_n2o5/Dp(k)
        fkn_o3 = fuchs(kn_o3,alpha_o3)		! = J/Jc
        fkn_n2o5 = fuchs(kn_n2o5,alpha_n2o5)
        sum_o3 = sum_o3 + Dp(k)*fkn_o3*Npcasp(k)
        sum_n2o5= sum_n2o5 + Dp(k)*fkn_n2o5*Npcasp(k)
      enddo
      
      rk_fuchs_o3   = 2.*3.14159*Dg_o3*sum_o3
      rk_fuchs_n2o5 = 2.*3.14159*Dg_n2o5*sum_n2o5
c
c
c
c uptake coefficient theory

      gamma_n2o5 = 0.00
      gamma_o3   = 1.0e-6
      gamma_no3  = 0.0
      
      sum_area = 0.0
      do k = 1, 15
        sum_area = sum_area + Npcasp(k)*0.25*3.14159*Dp(k)**2	! Dp in cm
      enddo
      
      rk_uptake_o3   = c_o3  *gamma_o3  *sum_area
      rk_uptake_n2o5 = c_n2o5*gamma_n2o5*sum_area
      rk_uptake_no3  = c_no3 *gamma_no3 *sum_area
      
c      write(6,*)'rk_uptake_o3 = ', rk_uptake_o3





c Heterogeneous chemistry
      rk_het(1) = rk_uptake_o3			! O3 -->    [1/s]
      rk_het(2)	= 0.0				! HO2 --> 0.5H2O2
      rk_het(3)	= 0.0				! NO2 --> 0.5HONO + 0.5HNO3
      rk_het(4)	= rk_uptake_n2o5		! N2O5 --> 2HNO3
      rk_het(5)	= 0.0				! HNO3 --> NO2
      rk_het(6)	= 0.0				! HNO3 --> NO
      rk_het(7)	= rk_uptake_no3			! NO3 --> NO + O2

      return
      end




      function fuchs(Kn, a)
      real Kn

      fuchs = 0.75*a*(1. + Kn)/(Kn**2 + Kn + 0.283*Kn*a + 0.75*a)

      return
      end



c**************************************************************************
c subroutine GasRateConstants_Mar: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_mar    = reaction rate constants for marine mechanism (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants_Mar
      include 'chm.inc'
      include 'gas.inc'

c abstraction reaction
c Hynes et al. (1986)
      rk_tot_num =       te * exp(-234./te) + 
     &             8.46e-10 * exp(7230./te) +
     &             2.68e-10 * exp(7810./te)
      rk_tot_den = 1.04e+11 * te + 88.1 * exp(7460./te)
      rk_tot	 = rk_tot_num/rk_tot_den
c
      rk_mar(1)   = 9.60e-12 * exp(-234./te) ! ch3sch3 + oh --> ch3sch2
      Babs       = rk_mar(1)/rk_tot
      Badd	 = 1. - Babs
      rk_mar(2)   = 1.40e-13 * exp(500./te)  ! ch3sch3 + no3 --> 
      rk_mar(3)   = 1.26e-11 * exp(409./te)  ! ch3sch3 + o3p --> 
c
c addition reaction
      rk_mar(4)   = Badd*rk_tot		     ! ch3sch3 + oh --> ch3s(oh)ch3
      rk_mar(5)   = 8.0e-12
      rk_mar(6)   = 1.8e-13
      rk_mar(7)   = 5.8e-11
      rk_mar(8)   = 1.0e-14
      rk_mar(9)   = 5.0e-12
      rk_mar(10)  = 1.8e-13
      rk_mar(11)  = 1.0e-15
      rk_mar(12)  = 1.0e-13
      rk_mar(13)  = 1.0e-15
      rk_mar(14)  = 1.6e-11
      rk_mar(15)  = 1.0e-13
      rk_mar(16)  = 2.5e+13 * exp(-8686./te) ! ch3so2 --> so2 + ch3o2
      rk_mar(17)  = 1.0e-14
      rk_mar(18)  = 5.0e-15
      rk_mar(19)  = 2.5e-13
      rk_mar(20)  = 2.5e-13
      rk_mar(21)  = 5.0e-11
      rk_mar(22)  = 2.6e-18
      rk_mar(23)  = 3.3
      rk_mar(24)  = 1.0e-11
      rk_mar(25)  = 5.5e-12
      rk_mar(26)  = 2.0e+17 * exp(-12626./te) ! ch3so3 --> h2so4 + ch3o2
      rk_mar(27)  = 3.0e-15
      rk_mar(28)  = 3.0e-15
      rk_mar(29)  = 5.0e-11
      rk_mar(30)  = 1.6e-15
c
      rk_mar(31)  = 2.5e-13	! ch3sch2oo + ch3so2 --> ch3so3 + ch3so2
      rk_mar(32)  = 8.6e-14	! 2ch3sch2oo --> .15mtf + 1.85ch3so2
c
c dry deposition
      rk_mar(33)  = 0.0 ! 2.0e-5	! 1/s
      rk_mar(34)  = 0.0 ! 2.0e-5
      rk_mar(35)  = 0.0 ! 2.0e-5

      return
      end
c**************************************************************************
c subroutine GasRateConstants_Urb: generates thermal rate coefficients 
c                   for the selected mechanism
c nomenclature:
c rk_urb    = reaction rate constants for hc1 mechanism    (molec-cc-s)
c te        = ambient atmospheric temperature (K)
c 
c author: Rahul A. Zaveri
c date  : february 1996
c
c-------------------------------------------------------------------------
      subroutine GasRateConstants_Urb
      include 'chm.inc'
      include 'gas.inc'

      rk_urb(1) = 8.1e-13                                                           
      rk_urb(2) = rk_photo(jphoto_aone)                                                
      rk_urb(3) = te**2*ARR(5.3e-18, -230.)                                         
      rk_urb(4) = rk_photo(jphoto_mgly)                                                  
      rk_urb(5) = 1.7e-11                                                           
      rk_urb(6) = ARR(1.4e-12, -1900.)                                              
      rk_urb(7) = ARR(1.2e-14, -2630.)      

      rk0 = 1.0e-28
      rnn = 0.8
      rki = 8.8e-12
      rmm = 0.0                                        
      rk_urb(8) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)    
                                 
      rk_urb(9) = ARR(4.2e-15, -1800.)                                              
      rk_urb(10) = ARR(8.9e-16, -392.)                                              
      rk_urb(11) = ARR(5.8e-12, 478.)                                               
      rk_urb(12) = ARR(2.9e-11, 255.)                                               
      rk_urb(13) = ARR(3.1e-13, -1010.)                                             
      rk_urb(14) = 2.5e-12                                                          
      rk_urb(15) = ARR(2.1e-12, 322.)                                               
      rk_urb(16) = ARR(1.7e-11, 116.)                                               
      rk_urb(17) = 8.1e-12                                                          
      rk_urb(18) = 4.1e-11                                                          
      rk_urb(19) = 2.2e-11                                                          
      rk_urb(20) = 1.4e-11                                                          
      rk_urb(21) = 3.0e-11                                                          
      rk_urb(22) = rk_photo(jphoto_open)                                                 
      rk_urb(23) = ARR(5.4e-17, -500.)                                              
      rk_urb(24) = rk_photo(jphoto_rooh)                                                    
      rk_urb(25) = ARR(3.8e-12, 200.)                                               
      rk_urb(26) = ARR(1.6e-11, -540.)                                              
      rk_urb(27) = rk_photo(jphoto_onit)                                               
      rk_urb(28) = 4.0e-12                                                          
      rk_urb(29) = 4.0e-12                                                          
      rk_urb(30) = 4.0e-12                                                          
      rk_urb(31) = 4.0e-12                                                          
      rk_urb(32) = 2.5e-12                                                          
      rk_urb(33) = 1.2e-12                                                          
      rk_urb(34) = 4.0e-12                                                          
      rk_urb(35) = 2.5e-12                                                          
      rk_urb(36) = ARR(1.7e-13, 1300.)                                              
      rk_urb(37) = ARR(1.2e-13, 1300.)                                              
      rk_urb(38) = ARR(1.7e-13, 1300.)                                              
      rk_urb(39) = ARR(1.7e-13, 1300.)                                              
      rk_urb(40) = rk_param(jro2)                                                   
      rk_urb(41) = rk_param(jano2)                                                  
      rk_urb(42) = rk_param(jnap)                                                   
      rk_urb(43) = rk_param(jxo2)                                                   
      rk_urb(44) = 1.0e-11    

      return
      end
      subroutine GasRates_Bio(s)
      include 'chm.inc'
      include 'gas.inc'
      dimension s(nmax)

      r_bio(1)  = rk_bio(1) *s(iisop)*s(ioh)
      r_bio(2)  = rk_bio(2) *s(iisop)*s(io3)
      r_bio(3)  = rk_bio(3) *s(iisop)*s(ino3)
      r_bio(4)  = rk_bio(4) *s(iisoprd)
      r_bio(5)  = rk_bio(5) *s(iisoprd)*s(ioh)
      r_bio(6)  = rk_bio(6) *s(iisoprd)*s(io3)
      r_bio(7)  = rk_bio(7) *s(iisoprd)*s(ino3)
      r_bio(8)  = rk_bio(8) *s(iisopp)*s(ino)
      r_bio(9)  = rk_bio(9) *s(iisopn)*s(ino)
      r_bio(10) = rk_bio(10)*s(iisopo2)*s(ino)
      r_bio(11) = rk_bio(11)*s(iisopp)*s(iho2)
      r_bio(12) = rk_bio(12)*s(iisopn)*s(iho2)
      r_bio(13) = rk_bio(13)*s(iisopo2)*s(iho2)
      r_bio(14) = rk_bio(14)*s(iisopp)
      r_bio(15) = rk_bio(15)*s(iisopn)
      r_bio(16) = rk_bio(16)*s(iisopo2)
      r_bio(17) = rk_bio(17)*s(iterp)*s(io1d)
      r_bio(18) = rk_bio(18)*s(iterp)*s(ioh)
      r_bio(19) = rk_bio(19)*s(iterp)*s(io3)
      r_bio(20) = rk_bio(20)*s(iterp)*s(ino3)
      
      return
      end

      subroutine GasRates_Com(s)
      include 'chm.inc'
      include 'gas.inc'
      dimension s(nmax)

      r_com(1) = rk_com(1)*s(ino2)                                                      
      r_com(2) = rk_com(2)*s(ino3)                                                      
      r_com(3) = rk_com(3)*s(ihono)                                                     
      r_com(4) = rk_com(4)*s(ihno3)                                                     
      r_com(5) = rk_com(5)*s(ihno4)                                                     
      r_com(6) = rk_com(6)*s(in2o5)                                                     
      r_com(7) = rk_com(7)*s(io3)                                                       
      r_com(8) = rk_com(8)*s(io3)                                                       
      r_com(9) = rk_com(9)*s(ih2o2)                                                     
      r_com(10) = rk_com(10)*s(io1d)*.21*cair_mlc                                                 
      r_com(11) = rk_com(11)*s(io1d)*.79*cair_mlc                                                 
      r_com(12) = rk_com(12)*s(io1d)*H2O                                                
      r_com(13) = rk_com(13)*s(io3p)*.21*cair_mlc                                                 
      r_com(14) = rk_com(14)*s(io3p)*s(io3)                                             
      r_com(15) = rk_com(15)*s(io3p)*s(ino2)                                            
      r_com(16) = rk_com(16)*s(io3p)*s(ino2)                                            
      r_com(17) = rk_com(17)*s(io3p)*s(ino)                                             
      r_com(18) = rk_com(18)*s(io3)*s(ino)                                              
      r_com(19) = rk_com(19)*s(io3)*s(ino2)                                             
      r_com(20) = rk_com(20)*s(io3)*s(ioh)                                              
      r_com(21) = rk_com(21)*s(io3)*s(iho2)                                             
      r_com(22) = rk_com(22)*s(ioh)*H2                                                  
      r_com(23) = rk_com(23)*s(ioh)*s(ino)                                              
      r_com(24) = rk_com(24)*s(ioh)*s(ino2)                                             
      r_com(25) = rk_com(25)*s(ioh)*s(ino3)                                             
      r_com(26) = rk_com(26)*s(ioh)*s(ihono)                                            
      r_com(27) = rk_com(27)*s(ioh)*s(ihno3)                                            
      r_com(28) = rk_com(28)*s(ioh)*s(ihno4)                                            
      r_com(29) = rk_com(29)*s(ioh)*s(iho2)                                             
      r_com(30) = rk_com(30)*s(ioh)*s(ih2o2)                                            
      r_com(31) = rk_com(31)*s(iho2)*s(iho2)                                            
      r_com(32) = rk_com(32)*s(iho2)*s(iho2)*H2O                                        
      r_com(33) = rk_com(33)*s(iho2)*s(ino)                                             
      r_com(34) = rk_com(34)*s(iho2)*s(ino2)                                            
      r_com(35) = rk_com(35)*s(iho2)*s(ino2)                                            
      r_com(36) = rk_com(36)*s(ihno4)                                                   
      r_com(37) = rk_com(37)*s(ino3)*s(ino)                                             
      r_com(38) = rk_com(38)*s(ino3)*s(ino2)                                            
      r_com(39) = rk_com(39)*s(ino3)*s(ino2)                                            
      r_com(40) = rk_com(40)*s(ino3)*s(ino3)                                            
      r_com(41) = rk_com(41)*s(ino3)*s(iho2)                                            
      r_com(42) = rk_com(42)*s(in2o5)*H2O                                               
      r_com(43) = rk_com(43)*s(in2o5)                                                   
      r_com(44) = rk_com(44)*s(ico)*s(ioh)                                              
      r_com(45) = rk_com(45)*s(iso2)*s(ioh)                                             
      r_com(46) = rk_com(46)*s(ich4)*s(ioh)                                             
      r_com(47) = rk_com(47)*s(ic2h6)*s(ioh)                                            
      r_com(48) = rk_com(48)*s(ich3oh)*s(ioh)                                           
      r_com(49) = rk_com(49)*s(ihcho)                                                   
      r_com(50) = rk_com(50)*s(ihcho)                                                   
      r_com(51) = rk_com(51)*s(ihcho)*s(ioh)                                            
      r_com(52) = rk_com(52)*s(ihcho)*s(ino3)                                           
      r_com(53) = rk_com(53)*s(ich3ooh)                                                 
      r_com(54) = rk_com(54)*s(iethooh)                                                 
      r_com(55) = rk_com(55)*s(ich3ooh)*s(ioh)                                          
      r_com(56) = rk_com(56)*s(iethooh)*s(ioh)                                          
      r_com(57) = rk_com(57)*s(ich3o2)*s(ino)                                           
      r_com(58) = rk_com(58)*s(iethp)*s(ino)                                            
      r_com(59) = rk_com(59)*s(ich3o2)*s(ino3)                                          
      r_com(60) = rk_com(60)*s(iethp)*s(ino3)                                           
      r_com(61) = rk_com(61)*s(ich3o2)*s(iho2)                                          
      r_com(62) = rk_com(62)*s(iethp)*s(iho2)                                           
      r_com(63) = rk_com(63)*s(ich3o2)                                                  
      r_com(64) = rk_com(64)*s(iethp)                                                   
      r_com(65) = rk_com(65)*s(ianol)*s(ioh)                                            
      r_com(66) = rk_com(66)*s(iald2)                                                   
      r_com(67) = rk_com(67)*s(iald2)*s(ioh)                                            
      r_com(68) = rk_com(68)*s(iald2)*s(ino3)                                           
      r_com(69) = rk_com(69)*s(ic2o3)*s(ino2)                                           
      r_com(70) = rk_com(70)*s(ipan)                                                    
      r_com(71) = rk_com(71)*s(ic2o3)*s(ino)                                            
      r_com(72) = rk_com(72)*s(ic2o3)*s(ino3)                                           
      r_com(73) = rk_com(73)*s(ic2o3)*s(iho2)                                           
      r_com(74) = rk_com(74)*s(ic2o3)                                                   

      return
      end
c**************************************************************************
c subroutine gasrate: calculates reaction rates for the selected mechanism
c
c nomenclature:
c r_com, r_urb, r_bio, r_mar = reaction rates (molec/cc/sec)
c rk_com,rk_urb,rk_bio,rk_mar= rate constants in appropriate units
c s                          = species concentrations (molec/cc)
c o2                         = oxygen concentration   (molec/cc)
c cair_mlc (used for M)      = air concentration      (molec/cc)
c h2o                        = water vapor            (molec/cc)
c
c author: Rahul A. Zaveri
c date  : february 1996
c
c--------------------------------------------------------------------------
      subroutine GasRates(s)
      include 'chm.inc'
      include 'gas.inc'

      dimension s(nmax)
c

      call gasrates_het(s)

      goto (1,2,3,4,5,6), iregime

1     call gasrates_com(s)
      return

2     call gasrates_com(s)
      call gasrates_urb(s)
      return

3     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_bio(s)
      return

4     call gasrates_com(s)
      call gasrates_mar(s)
      return
 
5     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_mar(s)
      return

6     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_bio(s)
      call gasrates_mar(s)
      return
c
      end



      subroutine GasRates_Het(s)
      include 'chm.inc'
      include 'gas.inc'
      dimension s(nmax)


c heterogeneous chemistry
      r_het(1) = rk_het(1)*s(io3)	! O3 --> 
      r_het(2) = rk_het(2)*s(iho2)	! HO2 --> 0.5H2O2
      r_het(3) = rk_het(3)*s(ino2)	! NO2 --> 0.5HONO + 0.5HNO3
      r_het(4) = rk_het(4)*s(in2o5)	! N2O5 --> 2HNO3
      r_het(5) = rk_het(5)*s(ihno3)	! HNO3 --> NO2
      r_het(6) = rk_het(6)*s(ihno3)	! HNO3 --> NO
      r_het(7) = rk_het(7)*s(ino3)	! NO3 --> NO

      return
      end



      subroutine GasRates_Mar(s)
      include 'chm.inc'
      include 'gas.inc'
      dimension s(nmax)

      r_mar(1)    = rk_mar(1)   *s(idms)*s(ioh)
      r_mar(2)    = rk_mar(2)   *s(idms)*s(ino3)
      r_mar(3)    = rk_mar(3)   *s(idms)*s(io3p)
      r_mar(4)    = rk_mar(4)   *s(idms)*s(ioh)
      r_mar(5)    = rk_mar(5)   *s(ich3sch2oo)*s(ino)
      r_mar(6)    = rk_mar(6)   *s(ich3sch2oo)*s(ich3o2)
      r_mar(7)    = rk_mar(7)   *s(idmso)*s(ioh)
      r_mar(8)    = rk_mar(8)   *s(idmso2)*s(ioh)
      r_mar(9)    = rk_mar(9)   *s(ich3so2ch2oo)*s(ino)
      r_mar(10)   = rk_mar(10)  *s(ich3so2ch2oo)*s(ich3o2)
      r_mar(11)   = rk_mar(11)  *s(ich3so2h)*s(iho2)
      r_mar(12)   = rk_mar(12)  *s(ich3so2h)*s(ino3)
      r_mar(13)   = rk_mar(13)  *s(ich3so2h)*s(ich3o2)
      r_mar(14)   = rk_mar(14)  *s(ich3so2h)*s(ioh)
      r_mar(15)   = rk_mar(15)  *s(ich3so2h)*s(ich3so3)
      r_mar(16)   = rk_mar(16)  *s(ich3so2)
      r_mar(17)   = rk_mar(17)  *s(ich3so2)*s(ino2)
      r_mar(18)   = rk_mar(18)  *s(ich3so2)*s(io3)
      r_mar(19)   = rk_mar(19)  *s(ich3so2)*s(iho2)
      r_mar(20)   = rk_mar(20)  *s(ich3so2)*s(ich3o2)
      r_mar(21)   = rk_mar(21)  *s(ich3so2)*s(ioh)
      r_mar(22)   = rk_mar(22)  *s(ich3so2)*o2
      r_mar(23)   = rk_mar(23)  *s(ich3so2oo)
      r_mar(24)   = rk_mar(24)  *s(ich3so2oo)*s(ino)
      r_mar(25)   = rk_mar(25)  *s(ich3so2oo)*s(ich3o2)
      r_mar(26)   = rk_mar(26)  *s(ich3so3)
      r_mar(27)   = rk_mar(27)  *s(ich3so3)*s(ino2)
      r_mar(28)   = rk_mar(28)  *s(ich3so3)*s(ino)
      r_mar(29)   = rk_mar(29)  *s(ich3so3)*s(iho2)
      r_mar(30)   = rk_mar(30)  *s(ich3so3)*s(ihcho)
c
      r_mar(31)   = rk_mar(31)	*s(ich3sch2oo)*s(ich3so2)
      r_mar(32)   = rk_mar(32)  *s(ich3sch2oo)*s(ich3sch2oo)
c
      r_mar(33)   = rk_mar(33)  *s(iso2)
      r_mar(34)   = rk_mar(34)  *s(idmso)
      r_mar(35)   = rk_mar(35)  *s(idmso2)
c
      return
      end

      subroutine GasRates_Urb(s)
      include 'chm.inc'
      include 'gas.inc'
      dimension s(nmax)

      r_urb(1) = rk_urb(1)*s(ipar)*s(ioh)                                               
      r_urb(2) = rk_urb(2)*s(iaone)                                                     
      r_urb(3) = rk_urb(3)*s(iaone)*s(ioh)                                              
      r_urb(4) = rk_urb(4)*s(imgly)                                                     
      r_urb(5) = rk_urb(5)*s(imgly)*s(ioh)                                              
      r_urb(6) = rk_urb(6)*s(imgly)*s(ino3)                                             
      r_urb(7) = rk_urb(7)*s(ieth)*s(io3)                                               
      r_urb(8) = rk_urb(8)*s(ieth)*s(ioh)                                               
      r_urb(9) = rk_urb(9)*s(iolet)*s(io3)                                              
      r_urb(10) = rk_urb(10)*s(iolei)*s(io3)                                            
      r_urb(11) = rk_urb(11)*s(iolet)*s(ioh)                                            
      r_urb(12) = rk_urb(12)*s(iolei)*s(ioh)                                            
      r_urb(13) = rk_urb(13)*s(iolet)*s(ino3)                                           
      r_urb(14) = rk_urb(14)*s(iolei)*s(ino3)                                           
      r_urb(15) = rk_urb(15)*s(itol)*s(ioh)                                             
      r_urb(16) = rk_urb(16)*s(ixyl)*s(ioh)                                             
      r_urb(17) = rk_urb(17)*s(ito2)*s(ino)                                             
      r_urb(18) = rk_urb(18)*s(icres)*s(ioh)                                            
      r_urb(19) = rk_urb(19)*s(icres)*s(ino3)                                           
      r_urb(20) = rk_urb(20)*s(icro)*s(ino2)                                            
      r_urb(21) = rk_urb(21)*s(iopen)*s(ioh)                                            
      r_urb(22) = rk_urb(22)*s(iopen)                                                   
      r_urb(23) = rk_urb(23)*s(iopen)*s(io3)                                            
      r_urb(24) = rk_urb(24)*s(irooh)                                                   
      r_urb(25) = rk_urb(25)*s(irooh)*s(ioh)                                            
      r_urb(26) = rk_urb(26)*s(ionit)*s(ioh)                                            
      r_urb(27) = rk_urb(27)*s(ionit)                                                   
      r_urb(28) = rk_urb(28)*s(iro2)*s(ino)                                             
      r_urb(29) = rk_urb(29)*s(iano2)*s(ino)                                            
      r_urb(30) = rk_urb(30)*s(inap)*s(ino)                                             
      r_urb(31) = rk_urb(31)*s(ixo2)*s(ino)                                             
      r_urb(32) = rk_urb(32)*s(iro2)*s(ino3)                                            
      r_urb(33) = rk_urb(33)*s(iano2)*s(ino3)                                           
      r_urb(34) = rk_urb(34)*s(inap)*s(ino3)                                            
      r_urb(35) = rk_urb(35)*s(ixo2)*s(ino3)                                            
      r_urb(36) = rk_urb(36)*s(iro2)*s(iho2)                                            
      r_urb(37) = rk_urb(37)*s(iano2)*s(iho2)                                           
      r_urb(38) = rk_urb(38)*s(inap)*s(iho2)                                            
      r_urb(39) = rk_urb(39)*s(ixo2)*s(iho2)                                            
      r_urb(40) = rk_urb(40)*s(iro2)                                                    
      r_urb(41) = rk_urb(41)*s(iano2)                                                   
      r_urb(42) = rk_urb(42)*s(inap)                                                    
      r_urb(43) = rk_urb(43)*s(ixo2)                                                    
      r_urb(44) = rk_urb(44)*s(ipar)*s(ixpar)                                           

      return
      end



c external dummy jacobian evaluation for LSODES (when mf=222)
      subroutine jac(ngas,tt,s,j,ian,jan,pdj)
      dimension s(1),ian(1),jan(1),pdj(1)
      return
      end
c**************************************************************************
c subroutine LoadPeroxyParameters: loads thermal rate coefficients 
c                                  for peroxy-peroxy permutation reactions 
c
c nomenclature:
c Aperox  = Pre-exponential factor (molec-cc-s)
c Bperox  = activation energy (-E/R)  (K)
c 
c author: Rahul A. Zaveri
c date  : June 1998
c
c-------------------------------------------------------------------------
      subroutine LoadPeroxyParameters
      include 'chm.inc'
      include 'gas.inc'

      Aperox(jch3o2,jch3o2)   = 2.5e-13
      Aperox(jethp,jethp)     = 6.8e-14
      Aperox(jc2o3,jc2o3)     = 2.9e-12
      Aperox(jano2,jano2)     = 8.0e-12
      Aperox(jnap,jnap)       = 1.0e-12
      Aperox(jro2,jro2)       = 5.3e-16
      Aperox(jisopp,jisopp)   = 3.1e-14
      Aperox(jisopn,jisopn)   = 3.1e-14
      Aperox(jisopo2,jisopo2) = 3.1e-14
      Aperox(jxo2,jxo2)       = 3.1e-14

      Bperox(jch3o2,jch3o2)   = 190.
      Bperox(jethp,jethp)     = 0.0
      Bperox(jc2o3,jc2o3)     = 500.
      Bperox(jano2,jano2)     = 0.0
      Bperox(jnap,jnap)       = 0.0
      Bperox(jro2,jro2)       = 1980.
      Bperox(jisopp,jisopp)   = 1000.
      Bperox(jisopn,jisopn)   = 1000.
      Bperox(jisopo2,jisopo2) = 1000.
      Bperox(jxo2,jxo2)       = 1000.

      do i = 1, nperox
      do j = 1, nperox
        if(i.ne.j)then
          Aperox(i,j) = 2.0*sqrt(Aperox(i,i)*Aperox(j,j))
          Bperox(i,j) = 0.5*(Bperox(i,i) + Bperox(j,j))
        endif
      enddo
      enddo
c
c except for
      Aperox(jc2o3,jch3o2) = 1.3e-12
      Aperox(jch3o2,jc2o3) = 1.3e-12
      Bperox(jc2o3,jch3o2) = 640.
      Bperox(jch3o2,jc2o3) = 640.
c
      return
      end

























c***************************************************************************
c   lsodeoct91.f - from netlib on 14-oct-91
c	obtained using the following requests
c		send lsode from sodepack
c		send xsetun, xsatf, srcom from sodepack
c		send sgefa, sgesl, sgbfa, sgbsl from linpack
c		send saxpy, scopy, sdot, sscal, isamax, r1mech from core
c----------------------------------------------------------------------

      subroutine lsode (f, neq, y, t, tout, itol, rtol, atol, itask,
     1            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      real y, t, tout, rtol, atol, rwork
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsode.. livermore solver for ordinary differential equations.
c this version is in single precision.
c
c lsode solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c lsode is a package based on the gear and gearb packages, and on the
c october 23, 1978 version of the tentative odepack user interface
c standard, with minor modifications.
c-----------------------------------------------------------------------
c reference..
c     alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c-----------------------------------------------------------------------
c author and contact.. alan c. hindmarsh,
c                      computing and mathematics research div., l-316
c                      lawrence livermore national laboratory
c                      livermore, ca 94550.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsode package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including optional communication, nonstandard options,
c and instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. next determine (or guess) whether or not the problem is stiff.
c stiffness occurs when the jacobian matrix df/dy has an eigenvalue
c whose real part is negative and large in magnitude, compared to the
c reciprocal of the t span of interest.  if the problem is nonstiff,
c use a method flag mf = 10.  if it is stiff, there are four standard
c choices for mf, and lsode requires the jacobian matrix in some form.
c this matrix is regarded either as full (mf = 21 or 22),
c or banded (mf = 24 or 25).  in the banded case, lsode requires two
c half-bandwidth parameters ml and mu.  these are, respectively, the
c widths of the lower and upper parts of the band, excluding the main
c diagonal.  thus the band consists of the locations (i,j) with
c i-ml .le. j .le. i+mu, and the full bandwidth is ml+mu+1.
c
c c. if the problem is stiff, you are encouraged to supply the jacobian
c directly (mf = 21 or 24), but if this is not feasible, lsode will
c compute it internally by difference quotients (mf = 22 or 25).
c if you are supplying the jacobian, provide a subroutine of the form..
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
c               dimension y(neq), pd(nrowpd,neq)
c which supplies df/dy by loading pd as follows..
c     for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),
c the partial derivative of f(i) with respect to y(j).  (ignore the
c ml and mu arguments in this case.)
c     for a banded jacobian (mf = 24), load pd(i-j+mu+1,j) with
c df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
c pd from the top down.
c     in either case, only nonzero elements need be loaded.
c
c d. write a main program which calls subroutine lsode once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsode.  on the first call to lsode, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be roughly less (in magnitude) than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             20 + 16*neq                    for mf = 10,
c             22 +  9*neq + neq**2           for mf = 21 or 22,
c             22 + 10*neq + (2*ml + mu)*neq  for mf = 24 or 25.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least..
c             20        for mf = 10,
c             20 + neq  for mf = 21, 22, 24, or 25.
c          if mf = 24 or 25, input in iwork(1),iwork(2) the lower
c          and upper half-bandwidths ml,mu.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix (mf = 21 or 24).
c          if used, this name must be declared external in calling
c          program.  if not used, pass a dummy name.
c mf     = method flag.  standard values are..
c          10 for nonstiff (adams) method, no jacobian used.
c          21 for stiff (bdf) method, user-supplied full jacobian.
c          22 for stiff method, internally generated full jacobian.
c          24 for stiff method, user-supplied banded jacobian.
c          25 for stiff method, internally generated banded jacobian.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c e. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsode was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong mf).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of mf or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c
c f. to continue the integration after a successful return, simply
c reset tout and call lsode again.  no other parameters need be reset.
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsode.  the problem is from chemical
c kinetics, and consists of the following three rate equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c     dy3/dt = 3.e7*y2**2
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsode, using mf = 21
c and printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 2 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
c     external fex, jex
c     dimension y(3), atol(3), rwork(58), iwork(23)
c     neq = 3
c     y(1) = 1.
c     y(2) = 0.
c     y(3) = 0.
c     t = 0.
c     tout = .4
c     itol = 2
c     rtol = 1.e-4
c     atol(1) = 1.e-6
c     atol(2) = 1.e-10
c     atol(3) = 1.e-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     lrw = 58
c     liw = 23
c     mf = 21
c     do 40 iout = 1,12
c       call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
c    1     iopt,rwork,lrw,iwork,liw,jex,mf)
c       write(6,20)t,y(1),y(2),y(3)
c 20    format(7h at t =,e12.4,6h   y =,3e14.6)
c       if (istate .lt. 0) go to 80
c 40    tout = tout*10.
c     write(6,60)iwork(11),iwork(12),iwork(13)
c 60  format(/12h no. steps =,i4,11h  no. f-s =,i4,11h  no. j-s =,i4)
c     stop
c 80  write(6,90)istate
c 90  format(///22h error halt.. istate =,i3)
c     stop
c     end
c
c     subroutine fex (neq, t, y, ydot)
c     dimension y(3), ydot(3)
c     ydot(1) = -.04*y(1) + 1.e4*y(2)*y(3)
c     ydot(3) = 3.e7*y(2)*y(2)
c     ydot(2) = -ydot(1) - ydot(3)
c     return
c     end
c
c     subroutine jex (neq, t, y, ml, mu, pd, nrpd)
c     dimension y(3), pd(nrpd,3)
c     pd(1,1) = -.04
c     pd(1,2) = 1.e4*y(3)
c     pd(1,3) = 1.e4*y(2)
c     pd(2,1) = .04
c     pd(2,3) = -pd(1,3)
c     pd(3,2) = 6.e7*y(2)
c     pd(2,2) = -pd(1,2) - pd(3,2)
c     return
c     end
c
c the output of this program (on a cdc-7600 in single precision)
c is as follows..
c
c   at t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
c   at t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
c   at t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
c   at t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
c   at t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
c   at t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
c   at t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
c   at t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
c   at t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
c   at t =  4.0000e+08   y =  5.494529e-06  2.197824e-11  9.999945e-01
c   at t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
c   at t =  4.0000e+10   y = -7.170586e-08 -2.868234e-13  1.000000e+00
c
c   no. steps = 330  no. f-s = 405  no. j-s =  69
c-----------------------------------------------------------------------
c full description of user interface to lsode.
c
c the user interface to lsode consists of the following parts.
c
c i.   the call sequence to subroutine lsode, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsode package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the lsode package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
c and those used for both input and output are
c     y, t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsode to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c
c the descriptions of the call arguments are as follows.
c
c f      = the name of the user-supplied subroutine defining the
c          ode system.  the system must be put in the first-order
c          form dy/dt = f(t,y), where f is a vector-valued function
c          of the scalar t and the vector y.  subroutine f is to
c          compute the function f.  it is to have the form
c               subroutine f (neq, t, y, ydot)
c               dimension y(1), ydot(1)
c          where neq, t, and y are input, and the array ydot = f(t,y)
c          is output.  y and ydot are arrays of length neq.
c          (in the dimension statement above, 1 is a dummy
c          dimension.. it can be replaced by any value.)
c          subroutine f should not alter y(1),...,y(neq).
c          f must be declared external in the calling program.
c
c          subroutine f may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in f) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y below.
c
c          if quantities computed in the f routine are needed
c          externally to lsode, an extra call to f should be made
c          for this purpose, for consistent and accurate results.
c          if only the derivative dy/dt is needed, use intdv instead.
c
c neq    = the size of the ode system (number of first order
c          ordinary differential equations).  used only for input.
c          neq may be decreased, but not increased, during the problem.
c          if neq is decreased (with istate = 3 on input), the
c          remaining components of y should be left undisturbed, if
c          these are to be accessed in f and/or jac.
c
c          normally, neq is a scalar, and it is generally referred to
c          as a scalar in this user interface description.  however,
c          neq may be an array, with neq(1) set to the system size.
c          (the lsode package accesses only neq(1).)  in either case,
c          this parameter is passed as the neq argument in all calls
c          to f and jac.  hence, if it is an array, locations
c          neq(2),... may be used to store other integer data and pass
c          it to f and/or jac.  subroutines f and/or jac must include
c          neq in a dimension statement in that case.
c
c y      = a real array for the vector of dependent variables, of
c          length neq or more.  used for both input and output on the
c          first call (istate = 1), and only for output on other calls.
c          on the first call, y must contain the vector of initial
c          values.  on output, y contains the computed solution vector,
c          evaluated at t.  if desired, the y array may be used
c          for other purposes between calls to the solver.
c
c          this array is passed as the y argument in all calls to
c          f and jac.  hence its length may exceed neq, and locations
c          y(neq+1),... may be used to store other real data and
c          pass it to f and/or jac.  (the lsode package accesses only
c          y(1),...,y(neq).)
c
c t      = the independent variable.  on input, t is used only on the
c          first call, as the initial point of the integration.
c          on output, after each call, t is the value at which a
c          computed solution y is evaluated (usually the same as tout).
c          on an error return, t is the farthest point reached.
c
c tout   = the next value of t at which a computed solution is desired.
c          used only for input.
c
c          when starting the problem (istate = 1), tout may be equal
c          to t for one call, then should .ne. t for the next call.
c          for the initial t, an input value of tout .ne. t is used
c          in order to determine the direction of the integration
c          (i.e. the algebraic sign of the step sizes) and the rough
c          scale of the problem.  integration in either direction
c          (forward or backward in t) is permitted.
c
c          if itask = 2 or 5 (one-step modes), tout is ignored after
c          the first call (i.e. the first call with tout .ne. t).
c          otherwise, tout is required on every call.
c
c          if itask = 1, 3, or 4, the values of tout need not be
c          monotone, but a value of tout which backs up is limited
c          to the current internal t interval, whose endpoints are
c          tcur - hu and tcur (see optional outputs, below, for
c          tcur and hu).
c
c itol   = an indicator for the type of error control.  see
c          description below under atol.  used only for input.
c
c rtol   = a relative error tolerance parameter, either a scalar or
c          an array of length neq.  see description below under atol.
c          input only.
c
c atol   = an absolute error tolerance parameter, either a scalar or
c          an array of length neq.  input only.
c
c             the input parameters itol, rtol, and atol determine
c          the error control performed by the solver.  the solver will
c          control the vector e = (e(i)) of estimated local errors
c          in y, according to an inequality of the form
c                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
c          and the rms-norm (root-mean-square norm) here is
c          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
c          is a vector of weights which must always be positive, and
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
c
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting
c          user-supplied routines for the setting of ewt and/or for
c          the norm calculation.  see part iv below.
c
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly.
c
c itask  = an index specifying the task to be performed.
c          input only.  itask has the following values and meanings.
c          1  means normal computation of output values of y(t) at
c             t = tout (by overshooting and interpolating).
c          2  means take one step only and return.
c          3  means stop at the first internal mesh point at or
c             beyond t = tout and return.
c          4  means normal computation of output values of y(t) at
c             t = tout but without overshooting t = tcrit.
c             tcrit must be input as rwork(1).  tcrit may be equal to
c             or beyond tout, but not behind it in the direction of
c             integration.  this option is useful if the problem
c             has a singularity at or beyond t = tcrit.
c          5  means take one step, without passing tcrit, and return.
c             tcrit must be input as rwork(1).
c
c          note..  if itask = 4 or 5 and the solver reaches tcrit
c          (within roundoff), it will return t = tcrit (exactly) to
c          indicate this (unless itask = 4 and tout comes before tcrit,
c          in which case answers at t = tout are returned first).
c
c istate = an index used for input and output to specify the
c          the state of the calculation.
c
c          on input, the values of istate are as follows.
c          1  means this is the first call for the problem
c             (initializations will be done).  see note below.
c          2  means this is not the first call, and the calculation
c             is to continue normally, with no change in any input
c             parameters except possibly tout and itask.
c             (if itol, rtol, and/or atol are changed between calls
c             with istate = 2, the new values will be used but not
c             tested for legality.)
c          3  means this is not the first call, and the
c             calculation is to continue normally, but with
c             a change in input parameters other than
c             tout and itask.  changes are allowed in
c             neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,
c             and any of the optional inputs except h0.
c             (see iwork description for ml and mu.)
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 1 on input.
c
c          on output, istate has the following values and meanings.
c           1  means nothing was done, as tout was equal to t with
c              istate = 1 on input.  (however, an internal counter was
c              set to detect and prevent repeated calls of this type.)
c           2  means the integration was performed successfully.
c          -1  means an excessive amount of work (more than mxstep
c              steps) was done on this call, before completing the
c              requested task, but the integration was otherwise
c              successful as far as t.  (mxstep is an optional input
c              and is normally 500.)  to continue, the user may
c              simply reset istate to a value .gt. 1 and call again
c              (the excess work step counter will be reset to 0).
c              in addition, the user may increase mxstep to avoid
c              this error return (see below on optional inputs).
c          -2  means too much accuracy was requested for the precision
c              of the machine being used.  this was detected before
c              completing the requested task, but the integration
c              was successful as far as t.  to continue, the tolerance
c              parameters must be reset, and istate must be set
c              to 3.  the optional output tolsf may be used for this
c              purpose.  (note.. if this condition is detected before
c              taking any steps, then an illegal input return
c              (istate = -3) occurs instead.)
c          -3  means illegal input was detected, before taking any
c              integration steps.  see written message for details.
c              note..  if the solver detects an infinite loop of calls
c              to the solver with illegal input, it will cause
c              the run to stop.
c          -4  means there were repeated error test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              the problem may have a singularity, or the input
c              may be inappropriate.
c          -5  means there were repeated convergence test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              this may be caused by an inaccurate jacobian matrix,
c              if one is being used.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol(i)=0.0)
c              was requested on a variable which has now vanished.
c              the integration was successful as far as t.
c
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
c          also, since a negative input value of istate will be
c          regarded as illegal, a negative output value requires the
c          user to change it, and possibly other inputs, before
c          calling the solver again.
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used on this call.  input only.
c          the optional inputs are listed separately below.
c          iopt = 0 means no optional inputs are being used.
c                   default values will be used in all cases.
c          iopt = 1 means one or more optional inputs are being used.
c
c rwork  = a real working array (single precision).
c          the length of rwork must be at least
c             20 + nyh*(maxord + 1) + 3*neq + lwm    where
c          nyh    = the initial value of neq,
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lwm   = 0             if miter = 0,
c          lwm   = neq**2 + 2    if miter is 1 or 2,
c          lwm   = neq + 2       if miter = 3, and
c          lwm   = (2*ml+mu+1)*neq + 2 if miter is 4 or 5.
c          (see the mf description for meth and miter.)
c          thus if maxord has its default value and neq is constant,
c          this length is..
c             20 + 16*neq                  for mf = 10,
c             22 + 16*neq + neq**2         for mf = 11 or 12,
c             22 + 17*neq                  for mf = 13,
c             22 + 17*neq + (2*ml+mu)*neq  for mf = 14 or 15,
c             20 +  9*neq                  for mf = 20,
c             22 +  9*neq + neq**2         for mf = 21 or 22,
c             22 + 10*neq                  for mf = 23,
c             22 + 10*neq + (2*ml+mu)*neq  for mf = 24 or 25.
c          the first 20 words of rwork are reserved for conditional
c          and optional inputs and optional outputs.
c
c          the following word in rwork is a conditional input..
c            rwork(1) = tcrit = critical value of t which the solver
c                       is not to overshoot.  required if itask is
c                       4 or 5, and ignored otherwise.  (see itask.)
c
c lrw    = the length of the array rwork, as declared by the user.
c          (this will be checked by the solver.)
c
c iwork  = an integer work array.  the length of iwork must be at least
c             20        if miter = 0 or 3 (mf = 10, 13, 20, 23), or
c             20 + neq  otherwise (mf = 11, 12, 14, 15, 21, 22, 24, 25).
c          the first few words of iwork are used for conditional and
c          optional inputs and optional outputs.
c
c          the following 2 words in iwork are conditional inputs..
c            iwork(1) = ml     these are the lower and upper
c            iwork(2) = mu     half-bandwidths, respectively, of the
c                       banded jacobian, excluding the main diagonal.
c                       the band is defined by the matrix locations
c                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu
c                       must satisfy  0 .le.  ml,mu  .le. neq-1.
c                       these are required if miter is 4 or 5, and
c                       ignored otherwise.  ml and mu may in fact be
c                       the band parameters for a matrix to which
c                       df/dy is only approximately equal.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note..  the work arrays must not be altered between calls to lsode
c for the same problem, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsode between calls, if
c desired (but not for use by f or jac).
c
c jac    = the name of the user-supplied routine (miter = 1 or 4) to
c          compute the jacobian matrix, df/dy, as a function of
c          the scalar t and the vector y.  it is to have the form
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
c               dimension y(1), pd(nrowpd,1)
c          where neq, t, y, ml, mu, and nrowpd are input and the array
c          pd is to be loaded with partial derivatives (elements of
c          the jacobian matrix) on output.  pd must be given a first
c          dimension of nrowpd.  t and y have the same meaning as in
c          subroutine f.  (in the dimension statement above, 1 is a
c          dummy dimension.. it can be replaced by any value.)
c               in the full matrix case (miter = 1), ml and mu are
c          ignored, and the jacobian is to be loaded into pd in
c          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
c               in the band matrix case (miter = 4), the elements
c          within the band are to be loaded into pd in columnwise
c          manner, with diagonal lines of df/dy loaded into the rows
c          of pd.  thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).
c          ml and mu are the half-bandwidth parameters (see iwork).
c          the locations in pd in the two triangular areas which
c          correspond to nonexistent matrix elements can be ignored
c          or loaded arbitrarily, as they are overwritten by lsode.
c               jac need not provide df/dy exactly.  a crude
c          approximation (possibly with a smaller bandwidth) will do.
c               in either case, pd is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          each call to jac is preceded by a call to f with the same
c          arguments neq, t, and y.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by f and not recomputed by jac,
c          if desired.  also, jac may alter the y array, if desired.
c          jac must be declared external in the calling program.
c               subroutine jac may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in jac) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y above.
c
c mf     = the method flag.  used only for input.  the legal values of
c          mf are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, and 25.
c          mf has decimal digits meth and miter.. mf = 10*meth + miter.
c          meth indicates the basic linear multistep method..
c            meth = 1 means the implicit adams method.
c            meth = 2 means the method based on backward
c                     differentiation formulas (bdf-s).
c          miter indicates the corrector iteration method..
c            miter = 0 means functional iteration (no jacobian matrix
c                      is involved).
c            miter = 1 means chord iteration with a user-supplied
c                      full (neq by neq) jacobian.
c            miter = 2 means chord iteration with an internally
c                      generated (difference quotient) full jacobian
c                      (using neq extra calls to f per df/dy value).
c            miter = 3 means chord iteration with an internally
c                      generated diagonal jacobian approximation.
c                      (using 1 extra call to f per df/dy evaluation).
c            miter = 4 means chord iteration with a user-supplied
c                      banded jacobian.
c            miter = 5 means chord iteration with an internally
c                      generated banded jacobian (using ml+mu+1 extra
c                      calls to f per df/dy evaluation).
c          if miter = 1 or 4, the user must supply a subroutine jac
c          (the name is arbitrary) as described above under jac.
c          for other values of miter, a dummy argument can be used.
c-----------------------------------------------------------------------
c optional inputs.
c
c the following is a list of the optional inputs provided for in the
c call sequence.  (see also part ii.)  for each such input variable,
c this table lists its name as used in this documentation, its
c location in the call sequence, its meaning, and the default value.
c the use of any of these inputs requires iopt = 1, and in that
c case all of these inputs are examined.  a value of zero for any
c of these optional inputs will cause the default value to be used.
c thus to use a subset of the optional inputs, simply preload
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
c then set those of interest to nonzero values.
c
c name    location      meaning and default value
c
c h0      rwork(5)  the step size to be attempted on the first step.
c                   the default value is determined by the solver.
c
c hmax    rwork(6)  the maximum absolute step size allowed.
c                   the default value is infinite.
c
c hmin    rwork(7)  the minimum absolute step size allowed.
c                   the default value is 0.  (this lower bound is not
c                   enforced on the final step before reaching tcrit
c                   when itask = 4 or 5.)
c
c maxord  iwork(5)  the maximum order to be allowed.  the default
c                   value is 12 if meth = 1, and 5 if meth = 2.
c                   if maxord exceeds the default value, it will
c                   be reduced to the default value.
c                   if maxord is changed during the problem, it may
c                   cause the current order to be reduced.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from lsode, the variables listed
c below are quantities related to the performance of lsode
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsode, and on any return with
c istate = -1, -2, -4, -5, or -6.  on an illegal input return
c (istate = -3), they will be unchanged from their existing values
c (if any), except possibly for tolsf, lenrw, and leniw.
c on any error return, outputs relevant to the error will be defined,
c as noted below.
c
c name    location      meaning
c
c hu      rwork(11) the step size in t last used (successfully).
c
c hcur    rwork(12) the step size to be attempted on the next step.
c
c tcur    rwork(13) the current value of the independent variable
c                   which the solver has actually reached, i.e. the
c                   current internal mesh point in t.  on output, tcur
c                   will always be at least as far as the argument
c                   t, but may be farther (if interpolation was done).
c
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
c                   computed when a request for too much accuracy was
c                   detected (istate = -3 if detected at the start of
c                   the problem, istate = -2 otherwise).  if itol is
c                   left unaltered but rtol and atol are uniformly
c                   scaled up by a factor of tolsf for the next call,
c                   then the solver is deemed likely to succeed.
c                   (the user may also ignore tolsf and alter the
c                   tolerance parameters in any other way appropriate.)
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far.
c
c nje     iwork(13) the number of jacobian evaluations (and of matrix
c                   lu decompositions) for the problem so far.
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c the following two arrays are segments of the rwork array which
c may also be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address in rwork, and its description.
c
c name    base address      description
c
c yh      21             the nordsieck history array, of size nyh by
c                        (nqcur + 1), where nyh is the initial value
c                        of neq.  for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.
c
c acor     lenrw-neq+1   array of size neq used for the accumulated
c                        corrections on each step, scaled on output
c                        to represent the estimated local error in y
c                        on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from lsode.
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsode.
c (the routines xsetun and xsatf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsode, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsatf(mflag)         set a flag to control the printing of
c                             messages by lsode.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcom(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsode (see part iii below).
c                             rsav must be a real array of length 218
c                             or more, and isav must be an integer
c                             array of length 41 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcom is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsode.
c
c   call intdv(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsode.
c
c the detailed instructions for using intdv are as follows.
c the form of the call is..
c
c   call intdv (t, k, rwork(21), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsode).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsode directly.  since nqcur .ge. 1, the first
c             derivative dy/dt is always available with intdv.
c rwork(21) = the base address of the history array yh.
c nyh       = column length of yh, equal to the initial value of neq.
c
c the output parameters are..
c
c dky       = a real array of length neq containing the computed value
c             of the k-th derivative of y(t).
c iflag     = integer flag, returned as 0 if k and t were legal,
c             -1 if k was illegal, and -2 if t was illegal.
c             on an error return, a message is also written.
c-----------------------------------------------------------------------
c part iii.  common blocks.
c
c if lsode is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsode,
c   (2) the two internal common blocks
c         /ls0001/  of length  257  (218 single precision words
c                         followed by 39 integer words),
c         /eh0001/  of length  2 (integer words).
c
c if lsode is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above two common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsode is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsode call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsode call for that problem.  to save and restore the common
c blocks, use subroutine srcom (see part ii above).
c
c note.. in this version of lsode, there are two data statements,
c in subroutines lsode and xarrwv, which load variables into these
c labeled common blocks.  on some systems, it may be necessary to
c move these to a separate block data subprogram.
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the lsode package which
c relate to the measurement of errors.  either routine can be
c replaced by a user-supplied version, if desired.  however, since such
c a replacement may have a major impact on performance, it should be
c done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) ewsat.
c the following subroutine is called just before each internal
c integration step, and sets the array of error weights, ewt, as
c described under itol/rtol/atol above..
c     subroutine ewsat (neq, itol, rtol, atol, ycur, ewt)
c where neq, itol, rtol, and atol are as in the lsode call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewsat.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewsat is passed to the
c vnerm routine (see below), and also used by lsode in the computation
c of the optional output imxer, the diagonal jacobian approximation,
c and the increments for difference quotient jacobians.
c
c in the user-supplied version of ewsat, it may be desirable to use
c the current values of derivatives of y.  derivatives up to order nq
c are available from the history array yh, described above under
c optional outputs.  in ewsat, yh is identical to the ycur array,
c extended to nq + 1 columns with a column length of nyh and scale
c factors of h**j/factorial(j).  on the first call for the problem,
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.
c the quantities nq, nyh, h, and nst can be obtained by including
c in ewsat the statements..
c     common /ls0001/ rls(218),ils(39)
c     nq = ils(35)
c     nyh = ils(14)
c     nst = ils(36)
c     h = rls(212)
c thus, for example, the current value of dy/dt can be obtained as
c ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
c unnecessary when nst = 0).
c
c (b) vnerm.
c the following is a real function routine which computes the weighted
c root-mean-square norm of a vector v..
c     d = vnerm (n, v, w)
c where..
c   n = the length of the vector,
c   v = real array of length n containing the vector,
c   w = real array of length n containing weights,
c   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
c vnerm is called with n = neq and with w(i) = 1.0/ewt(i), where
c ewt is as set by subroutine ewsat.
c
c if the user supplies this function, it should return a non-negative
c value of vnerm suitable for use in the error control in lsode.
c none of the arguments should be altered by vnerm.
c for example, a user-supplied vnerm routine might..
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsode package.
c
c in addition to subroutine lsode, the lsode package includes the
c following subroutines and function routines..
c  intdv    computes an interpolated value of the y vector at t = tout.
c  stade    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfade    sets all method coefficients and test constants.
c  prepj    computes and preprocesses the jacobian matrix j = df/dy
c           and the newton iteration matrix p = i - h*l0*j.
c  solsy    manages solution of linear system in chord iteration.
c  ewsat    sets the error weight vector ewt before each step.
c  vnerm    computes the weighted r.m.s. norm of a vector.
c  srcom    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  sgefa and sgesl   are routines from linpack for solving full
c           systems of linear algebraic equations.
c  sgbfa and sgbsl   are routines from linpack for solving banded
c           linear systems.
c  saxpy, sscal, isamax, and sdot   are basic linear algebra modules
c           (blas) used by the above linpack routines.
c  r1mech   computes the unit roundoff in a machine-independent manner.
c  xarrwv, xsetun, and xsatf   handle the printing of all error
c           messages and warnings.  xarrwv is machine-dependent.
c note..  vnerm, isamax, sdot, and r1mech are function routines.
c all the others are subroutines.
c
c the intrinsic and external routines used by lsode are..
c abs, amax1, amin1, float, max0, min0, mod, sign, sqrt, and write.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on llnl compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prepj, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   r1mech, vnerm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsode first,
c then those local to subroutine stade, and finally those used
c for communication.  the block is declared in subroutines
c lsode, intdv, stade, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c      data illin/0/, ntrep/0/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 5) go to 608
      if (miter .le. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0e0
      hmxi = 0.0e0
      hmin = 0.0e0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0e0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0e0) go to 615
      hmxi = 0.0e0
      if (hmax .gt. 0.0e0) hmxi = 1.0e0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0e0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2
      if (miter .eq. 3) lenwm = n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0e0) go to 619
        if (atoli .lt. 0.0e0) go to 620
 70     continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stade. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
      do 80 i = 1,n
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0e0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = r1mech(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0e0) go to 625
      if (h0 .ne. 0.0e0 .and. (t + h0 - tcrit)*h0 .gt. 0.0e0)
     1   h0 = tcrit - t
 110  jstart = 0
      if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0e0
      nqu = 0
      ccmax = 0.3e0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0e0
      call ewsat (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 120    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0e0) go to 180
      tdist = abs(tout - t)
      w0 = amax1(abs(t),abs(tout))
      if (tdist .lt. 2.0e0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = amax1(tol,rtol(i))
 140  if (tol .gt. 0.0e0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0e0) tol = amax1(tol,atoli/ayi)
 150    continue
 160  tol = amax1(tol,100.0e0*uround)
      tol = amin1(tol,0.001e0)
      sum = vnerm (n, rwork(lf0), rwork(lewt))
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = amin1(h0,tdist)
      h0 = sign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdv (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0e0 + 100.0e0*uround)
      if ((tp - tout)*h .gt. 0.0e0) go to 623
      if ((tn - tout)*h .lt. 0.0e0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
      if ((tcrit - tout)*h .lt. 0.0e0) go to 625
      if ((tn - tout)*h .lt. 0.0e0) go to 245
      call intdv (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
 245  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stade.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewsat (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 510
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 270  tolsf = uround*vnerm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0e0) go to 280
      tolsf = tolsf*2.0e0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xarrwv(50hlsode--  warning..internal t (=r1) and h (=r2) are,
     1   50, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(
     1  60h      such that in the machine, t + h = t on the next step  ,
     1   60, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      (h = step size). solver will continue anyway,
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xarrwv(50hlsode--  above warning has been issued i1 times.  ,
     1   50, 102, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      it will not be issued again for this problem,
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
 290  continue
c-----------------------------------------------------------------------
c     call stade(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prepj,solsy)
c-----------------------------------------------------------------------
      call stade (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prepj, solsy)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdv (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0e0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0e0) go to 345
      call intdv (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsode.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xarrwv(
     1  60hlsode--  repeated calls with istate = 1 and tout = t (=r1)  ,
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0e0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xarrwv(50hlsode--  at current t (=r1), mxstep (=i1) steps   ,
     1   50, 201, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      taken on this call before reaching tout     ,
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0e0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xarrwv(50hlsode--  at t (=r1), ewt(i1) has become r2 .le. 0.,
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xarrwv(50hlsode--  at t (=r1), too much accuracy requested  ,
     1   50, 203, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      for precision of machine..  see tolsf (=r2) ,
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xarrwv(50hlsode--  at t(=r1) and step size h(=r2), the error,
     1   50, 204, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      test failed repeatedly or with abs(h) = hmin,
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xarrwv(50hlsode--  at t (=r1) and step size h (=r2), the    ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(50h      corrector convergence failed repeatedly     ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(30h      or with abs(h) = hmin   ,
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0e0
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xarrwv(30hlsode--  istate (=i1) illegal ,
     1   30, 1, 0, 1, istate, 0, 0, 0.0e0, 0.0e0)
      go to 700
 602  call xarrwv(30hlsode--  itask (=i1) illegal  ,
     1   30, 2, 0, 1, itask, 0, 0, 0.0e0, 0.0e0)
      go to 700
 603  call xarrwv(50hlsode--  istate .gt. 1 but lsode not initialized  ,
     1   50, 3, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      go to 700
 604  call xarrwv(30hlsode--  neq (=i1) .lt. 1     ,
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0e0, 0.0e0)
      go to 700
 605  call xarrwv(50hlsode--  istate = 3 and neq increased (i1 to i2)  ,
     1   50, 5, 0, 2, n, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 606  call xarrwv(30hlsode--  itol (=i1) illegal   ,
     1   30, 6, 0, 1, itol, 0, 0, 0.0e0, 0.0e0)
      go to 700
 607  call xarrwv(30hlsode--  iopt (=i1) illegal   ,
     1   30, 7, 0, 1, iopt, 0, 0, 0.0e0, 0.0e0)
      go to 700
 608  call xarrwv(30hlsode--  mf (=i1) illegal     ,
     1   30, 8, 0, 1, mf, 0, 0, 0.0e0, 0.0e0)
      go to 700
 609  call xarrwv(50hlsode--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2),
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 610  call xarrwv(50hlsode--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2),
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 611  call xarrwv(30hlsode--  maxord (=i1) .lt. 0  ,
     1   30, 11, 0, 1, maxord, 0, 0, 0.0e0, 0.0e0)
      go to 700
 612  call xarrwv(30hlsode--  mxstep (=i1) .lt. 0  ,
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0e0, 0.0e0)
      go to 700
 613  call xarrwv(30hlsode--  mxhnil (=i1) .lt. 0  ,
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
      go to 700
 614  call xarrwv(40hlsode--  tout (=r1) behind t (=r2)      ,
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xarrwv(50h      integration direction is given by h0 (=r1)  ,
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0e0)
      go to 700
 615  call xarrwv(30hlsode--  hmax (=r1) .lt. 0.0  ,
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0e0)
      go to 700
 616  call xarrwv(30hlsode--  hmin (=r1) .lt. 0.0  ,
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0e0)
      go to 700
 617  call xarrwv(
     1  60hlsode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2),
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 618  call xarrwv(
     1  60hlsode--  iwork length needed, leniw (=i1), exceeds liw (=i2),
     1   60, 18, 0, 2, leniw, liw, 0, 0.0e0, 0.0e0)
      go to 700
 619  call xarrwv(40hlsode--  rtol(i1) is r1 .lt. 0.0        ,
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0e0)
      go to 700
 620  call xarrwv(40hlsode--  atol(i1) is r1 .lt. 0.0        ,
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0e0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xarrwv(40hlsode--  ewt(i1) is r1 .le. 0.0         ,
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0e0)
      go to 700
 622  call xarrwv(
     1  60hlsode--  tout (=r1) too close to t(=r2) to start integration,
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xarrwv(
     1  60hlsode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ,
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xarrwv(
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ,
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xarrwv(
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ,
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xarrwv(50hlsode--  at start of problem, too much accuracy   ,
     1   50, 26, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xarrwv(
     1  60h      requested for precision of machine..  see tolsf (=r1) ,
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0e0)
      rwork(14) = tolsf
      go to 700
 627  call xarrwv(50hlsode--  trouble from intdv. itask = i1, tout = r1,
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0e0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xarrwv(50hlsode--  repeated occurrences of illegal input    ,
     1   50, 302, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
c
 800  call xarrwv(50hlsode--  run aborted.. apparent infinite loop     ,
     1   50, 303, 2, 0, 0, 0, 0, 0.0e0, 0.0e0)
      return
c----------------------- end of subroutine lsode -----------------------
      end
      subroutine prepj (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm,
     1   f, jac)
clll. optimize
      external f, jac
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,
     1   mba, mband, meb1, meband, ml, ml3, mu, np1
      real y, yh, ewt, ftem, savf, wm
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con, di, fac, hl0, r, r0, srur, yi, yj, yjj,
     1   vnerm
      dimension neq(1), y(1), yh(nyh,1), ewt(1), ftem(1), savf(1),
     1   wm(1), iwm(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c prepj is called by stade to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c here j is computed by the user-supplied routine jac if
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c if miter = 3, a diagonal approximation to j is used.
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c subjected to lu decomposition in preparation for later solution
c of linear systems with p as coefficient matrix. this is done
c by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with prepj uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stade).
c savf  = array containing f evaluated at predicted y.
c wm    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3 and the lu decomposition
c         of p if miter is 1, 2 , 4, or 5.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround), used in numerical jacobian increments.
c         wm(2) = h*el0, saved for later use if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0   = el(1) (input).
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
c         p matrix found to be singular.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje.
c-----------------------------------------------------------------------
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0e0
      call jac (neq, tn, y, 0, 0, wm(3), n)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnerm (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = amax1(srur*abs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call f (neq, tn, y, ftem)
        do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      np1 = n + 1
      do 250 i = 1,n
        wm(j) = wm(j) + 1.0e0
 250    j = j + np1
c do lu decomposition on p. --------------------------------------------
      call sgefa (wm(3), n, n, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(2) = hl0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wm(3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1e0*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = 1.0e0
        if (abs(r0) .lt. uround/ewt(i)) go to 320
        if (abs(di) .eq. 0.0e0) go to 330
        wm(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ierpj = 1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0e0
      call jac (neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnerm (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = amax1(srur*abs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        call f (neq, tn, y, ftem)
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = amax1(srur*abs(yjj),r0/ewt(jj))
          fac = -hl0/r
          i1 = max0(jj-mu,1)
          i2 = min0(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
        wm(ii) = wm(ii) + 1.0e0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call sgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c----------------------- end of subroutine prepj -----------------------
      end
      subroutine solsy (wm, iwm, x, tem)
clll. optimize
      integer iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real wm, x, tem
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real di, hl0, phl0, r
      dimension wm(1), iwm(1), x(1), tem(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls sgesl to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c if miter is 4 or 5, it calls sgbsl.
c communication with solsy uses the following variables..
c wm    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround) (not used here),
c         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.
c         iersl = 1 if a singular matrix arose with miter = 3.
c this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
      iersl = 0
      go to (100, 100, 300, 400, 400), miter
 100  call sgesl (wm(3), n, n, iwm(21), x, 0)
      return
c
 300  phl0 = wm(2)
      hl0 = h*el0
      wm(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0e0 - r*(1.0e0 - 1.0e0/wm(i+2))
        if (abs(di) .eq. 0.0e0) go to 390
 320    wm(i+2) = 1.0e0/di
 330  do 340 i = 1,n
 340    x(i) = wm(i+2)*x(i)
      return
 390  iersl = 1
      return
c
 400  ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      call sgbsl (wm(3), meband, n, ml, mu, iwm(21), x, 0)
      return
c----------------------- end of subroutine solsy -----------------------
      end
      subroutine stade (neq, y, yh, nyh, yh1, ewt, savf, acor,
     1   wm, iwm, f, jac, pjac, slvs)
clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real y, yh, yh1, ewt, savf, acor, wm
      real conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnerm
      dimension neq(1), y(1), yh(nyh,1), yh1(1), ewt(1), savf(1),
     1   acor(1), wm(1), iwm(1)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c stade performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stade is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stade is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0e0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0e0
      rc = 0.0e0
      el0 = 1.0e0
      crate = 0.7e0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfade is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfade (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      ddn = vnerm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/float(l)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = amin1(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = amin1(rh,abs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfade is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfade (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = amax1(rh,hmin/abs(h))
 175  rh = amin1(rh,rmax)
      rh = rh/amax1(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnerm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnerm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = amax1(0.2e0*crate,del/delp)
      dcon = del*amin1(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnerm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnerm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnerm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = amin1(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = amax1(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0e0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stade -----------------------
      end
      real function vnerm (n, v, w)
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnerm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
c-----------------------------------------------------------------------
      integer n,   i
      real v, w,   sum
      dimension v(n), w(n)
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnerm = sqrt(sum/float(n))
      return
c----------------------- end of function vnerm -------------------------
      end
      subroutine xarrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      integer msg, nmes, nerr, level, ni, i1, i2, nr,
     1   i, lun, lunit, mesflg, ncpw, nch, nwds
      real r1, r2
      dimension msg(nmes)
c-----------------------------------------------------------------------
c subroutines xarrwv, xsatf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in single precision and are printed
c    in e21.13 format.
c 5. the common block /eh0001/ below is data-loaded (a machine-
c    dependent feature) with default values.
c    this block is needed for proper retention of parameters used by
c    this routine which the user can reset by calling xsatf or xsetun.
c    the variables in this block are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement below.
c
c for some systems, the data statement below must be replaced
c by a separate block data subprogram.
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh0001/ mesflg, lunit
c
c      data mesflg/1/, lunit/6/
c-----------------------------------------------------------------------
c the following data-loaded value of ncpw is valid for the cdc-6600
c and cdc-7600 computers.
c     data ncpw/10/
c the following is valid for the cray-1 computer.
c     data ncpw/8/
c the following is valid for the burroughs 6700 and 7800 computers.
c     data ncpw/6/
c the following is valid for the pdp-10 computer.
c     data ncpw/5/
c the following is valid for the vax computer with 4 bytes per integer,
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
      data ncpw/4/
c the following is valid for the pdp-11, or vax with 2-byte integers.
c     data ncpw/2/
c-----------------------------------------------------------------------
c
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, 10) (msg(i),i=1,nwds)
c-----------------------------------------------------------------------
c the following format statement is to have the form
c 10  format(1x,mmann)
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
c the following is valid for ncpw = 10.
c 10  format(1x,6a10)
c the following is valid for ncpw = 8.
c 10  format(1x,8a8)
c the following is valid for ncpw = 6.
c 10  format(1x,10a6)
c the following is valid for ncpw = 5.
c 10  format(1x,12a5)
c the following is valid for ncpw = 4.
  10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,e21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,e21.13,3x,4hr2 =,e21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c----------------------- end of subroutine xarrwv ----------------------
      end
      subroutine xsatf (mflag)
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsatf -----------------------
      end
      subroutine xsetun (lun)
c
c this routine resets the logical unit number for messages.
c
      integer lun, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (lun .gt. 0) lunit = lun
      return
c----------------------- end of subroutine xsetun ----------------------
      end

      subroutine srcom (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001 and eh0001, which are used internally
c by one or more odepack solvers.
c
c rsav = real array of length 218 or more.
c isav = integer array of length 41 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils
      integer i, lenils, lenrls
      real rsav,   rls
      dimension rsav(1), isav(1)
      common /ls0001/ rls(218), ils(39)
      common /eh0001/ ieh(2)
      data lenrls/218/, lenils/39/
c
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      return
c----------------------- end of subroutine srcom -----------------------
      end

      subroutine cfade (meth, elco, tesco)
clll. optimize
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
c-----------------------------------------------------------------------
c cfade is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for all orders and saved.
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfade is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array contains the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the implicit adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array contains test constants used for the
c local error test and the selection of step size and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c nq + 1 if k = 3.
c-----------------------------------------------------------------------
      dimension pc(12)
c
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0e0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0e0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0e0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = float(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0e0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
c----------------------- end of subroutine cfade -----------------------
      end
      subroutine ewsat (n, itol, rtol, atol, ycur, ewt)
clll. optimize
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
      integer n, itol
      integer i
      real rtol, atol, ycur, ewt
      dimension rtol(1), atol(1), ycur(n), ewt(n)
c
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine ewsat -----------------------
      end
      subroutine intdv (t, k, yh, nyh, dky, iflag)
clll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp
      dimension yh(nyh,1), dky(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c intdv computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0e0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xarrwv(30hintdv--  k (=i1) illegal      ,
     1   30, 51, 0, 1, k, 0, 0, 0.0e0, 0.0e0)
      iflag = -1
      return
 90   call xarrwv(30hintdv--  t (=r1) illegal      ,
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0e0)
      call xarrwv(
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine intdv -----------------------
      end

c----------------------------------------------------------------------------
c   the following routines are from the netlib linpack library
c----------------------------------------------------------------------------
      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      real abd(lda,1)
c
c     sgbfa factors a real band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c     fortran max0,min0
c
c     internal variables
c
      real t
      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0e0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0e0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0e0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0e0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call saxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0e0) info = n
      return
      end

c----------------------------------------------------------------------------
      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      real abd(lda,1),b(1)
c
c     sgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on entry
c
c        abd     real(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c     fortran min0
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call saxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call saxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end

c----------------------------------------------------------------------------
      subroutine sgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      real a(lda,1)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end

c----------------------------------------------------------------------------
      subroutine sgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      real a(lda,1),b(1)
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

c----------------------------------------------------------------------------
c   the following routines are from the netlib core library
c----------------------------------------------------------------------------

      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      real sx(1),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      smax = abs(sx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end

c---------------------------------------------------------------------------
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end

c---------------------------------------------------------------------------
      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end

c----------------------------------------------------------------------------
      real function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end

c---------------------------------------------------------------------------
      subroutine sscal(n,sa,sx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      real sa,sx(1)
      integer i,incx,ix,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      do 10 i = 1,n 
        sx(ix) = sa*sx(ix)
        ix = ix + incx 
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end

c----------------------------------------------------------------------------
      real function r1mech(i)
c
c  single-precision machine constants
c
c  r1mech(1) = b**(emin-1), the smallest positive magnitude.
c
c  r1mech(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  r1mech(3) = b**(-t), the smallest relative spacing.
c
c  r1mech(4) = b**(1-t), the largest relative spacing.
c
c  r1mech(5) = log10(b)
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.
c  on rare machines a static statement may need to be added.
c  (but probably more systems prohibit it than require it.)
c
c  for ieee-arithmetic machines (binary standard), the first
c  set of constants below should be appropriate.
c
c  where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  sometimes this requires using
c  equivalent integer arrays.  if your compiler uses half-word
c  integers by default (sometimes called integer*2), you may need to
c  change integer to integer*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc
c
      real rmach(5)
c
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
c
c     machine constants for ieee arithmetic machines, such as the at&t
c     3b series, motorola 68000 based machines (e.g. sun 3 and at&t
c     pc 7300), and 8087 based micros (e.g. ibm pc and at&t 6300).
c
       data small(1) /     8388608 /
       data large(1) /  2139095039 /
       data right(1) /   864026624 /
       data diver(1) /   872415232 /
       data log10(1) /  1050288283 /, sc/987/
c
c     machine constants for amdahl machines.
c
c      data small(1) /    1048576 /
c      data large(1) / 2147483647 /
c      data right(1) /  990904320 /
c      data diver(1) / 1007681536 /
c      data log10(1) / 1091781651 /, sc/987/
c
c     machine constants for the burroughs 1700 system.
c
c      data rmach(1) / z400800000 /
c      data rmach(2) / z5ffffffff /
c      data rmach(3) / z4e9800000 /
c      data rmach(4) / z4ea800000 /
c      data rmach(5) / z500e730e8 /, sc/987/
c
c     machine constants for the burroughs 5700/6700/7700 systems.
c
c      data rmach(1) / o1771000000000000 /
c      data rmach(2) / o0777777777777777 /
c      data rmach(3) / o1311000000000000 /
c      data rmach(4) / o1301000000000000 /
c      data rmach(5) / o1157163034761675 /, sc/987/
c
c     machine constants for ftn4 on the cdc 6000/7000 series.
c
c      data rmach(1) / 00564000000000000000b /
c      data rmach(2) / 37767777777777777776b /
c      data rmach(3) / 16414000000000000000b /
c      data rmach(4) / 16424000000000000000b /
c      data rmach(5) / 17164642023241175720b /, sc/987/
c
c     machine constants for ftn5 on the cdc 6000/7000 series.
c
c      data rmach(1) / o"00564000000000000000" /
c      data rmach(2) / o"37767777777777777776" /
c      data rmach(3) / o"16414000000000000000" /
c      data rmach(4) / o"16424000000000000000" /
c      data rmach(5) / o"17164642023241175720" /, sc/987/
c
c     machine constants for convex c-1.
c
c      data rmach(1) / '00800000'x /
c      data rmach(2) / '7fffffff'x /
c      data rmach(3) / '34800000'x /
c      data rmach(4) / '35000000'x /
c      data rmach(5) / '3f9a209b'x /, sc/987/
c
c     machine constants for the cray 1, xmp, 2, and 3.
c
c      data rmach(1) / 200034000000000000000b /
c      data rmach(2) / 577767777777777777776b /
c      data rmach(3) / 377224000000000000000b /
c      data rmach(4) / 377234000000000000000b /
c      data rmach(5) / 377774642023241175720b /, sc/987/
c
c     machine constants for the data general eclipse s/200.
c
c     note - it may be appropriate to include the following line -
c     static rmach(5)
c
c      data small/20k,0/,large/77777k,177777k/
c      data right/35420k,0/,diver/36020k,0/
c      data log10/40423k,42023k/, sc/987/
c
c     machine constants for the harris slash 6 and slash 7.
c
c      data small(1),small(2) / '20000000, '00000201 /
c      data large(1),large(2) / '37777777, '00000177 /
c      data right(1),right(2) / '20000000, '00000352 /
c      data diver(1),diver(2) / '20000000, '00000353 /
c      data log10(1),log10(2) / '23210115, '00000377 /, sc/987/
c
c     machine constants for the honeywell dps 8/70 series.
c
c      data rmach(1) / o402400000000 /
c      data rmach(2) / o376777777777 /
c      data rmach(3) / o714400000000 /
c      data rmach(4) / o716400000000 /
c      data rmach(5) / o776464202324 /, sc/987/
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c      data rmach(1) / z00100000 /
c      data rmach(2) / z7fffffff /
c      data rmach(3) / z3b100000 /
c      data rmach(4) / z3c100000 /
c      data rmach(5) / z41134413 /, sc/987/
c
c     machine constants for the interdata 8/32
c     with the unix system fortran 77 compiler.
c
c     for the interdata fortran vii compiler replace
c     the z's specifying hex constants with y's.
c
c      data rmach(1) / z'00100000' /
c      data rmach(2) / z'7effffff' /
c      data rmach(3) / z'3b100000' /
c      data rmach(4) / z'3c100000' /
c      data rmach(5) / z'41134413' /, sc/987/
c
c     machine constants for the pdp-10 (ka or ki processor).
c
c      data rmach(1) / "000400000000 /
c      data rmach(2) / "377777777777 /
c      data rmach(3) / "146400000000 /
c      data rmach(4) / "147400000000 /
c      data rmach(5) / "177464202324 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers (expressed in integer and octal).
c
c      data small(1) /    8388608 /
c      data large(1) / 2147483647 /
c      data right(1) /  880803840 /
c      data diver(1) /  889192448 /
c      data log10(1) / 1067065499 /, sc/987/
c
c      data rmach(1) / o00040000000 /
c      data rmach(2) / o17777777777 /
c      data rmach(3) / o06440000000 /
c      data rmach(4) / o06500000000 /
c      data rmach(5) / o07746420233 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     16-bit integers  (expressed in integer and octal).
c
c      data small(1),small(2) /   128,     0 /
c      data large(1),large(2) / 32767,    -1 /
c      data right(1),right(2) / 13440,     0 /
c      data diver(1),diver(2) / 13568,     0 /
c      data log10(1),log10(2) / 16282,  8347 /, sc/987/
c
c      data small(1),small(2) / o000200, o000000 /
c      data large(1),large(2) / o077777, o177777 /
c      data right(1),right(2) / o032200, o000000 /
c      data diver(1),diver(2) / o032400, o000000 /
c      data log10(1),log10(2) / o037632, o020233 /, sc/987/
c
c     machine constants for the sequent balance 8000.
c
c      data small(1) / $00800000 /
c      data large(1) / $7f7fffff /
c      data right(1) / $33800000 /
c      data diver(1) / $34000000 /
c      data log10(1) / $3e9a209b /, sc/987/
c
c     machine constants for the univac 1100 series.
c
c      data rmach(1) / o000400000000 /
c      data rmach(2) / o377777777777 /
c      data rmach(3) / o146400000000 /
c      data rmach(4) / o147400000000 /
c      data rmach(5) / o177464202324 /, sc/987/
c
c     machine constants for the vax unix f77 compiler.
c
c      data small(1) /       128 /
c      data large(1) /    -32769 /
c      data right(1) /     13440 /
c      data diver(1) /     13568 /
c      data log10(1) / 547045274 /, sc/987/
c
c     machine constants for the vax-11 with
c     fortran iv-plus compiler.
c
c      data rmach(1) / z00000080 /
c      data rmach(2) / zffff7fff /
c      data rmach(3) / z00003480 /
c      data rmach(4) / z00003500 /
c      data rmach(5) / z209b3f9a /, sc/987/
c
c     machine constants for vax/vms version 2.2.
c
c      data rmach(1) /       '80'x /
c      data rmach(2) / 'ffff7fff'x /
c      data rmach(3) /     '3480'x /
c      data rmach(4) /     '3500'x /
c      data rmach(5) / '209b3f9a'x /, sc/987/
c
c  ***  issue stop 778 if all data statements are commented...
      if (sc .ne. 987) stop 778
      if (i .lt. 1  .or.  i .gt. 5) goto 999
      r1mech = rmach(i)
      return
  999 write(*,1999) i
 1999 format(' r1mech - i out of bounds',i10)
      stop
      end

c Obtained Oct 11, 1994 from ODEPACK in NETLIB by RDS
      subroutine lsodes (f, neq, y, t, tout, itol, rtol, atol, itask,
     1            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      real y, t, tout, rtol, atol, rwork
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsodes.. livermore solver for ordinary differential equations
c          with general sparse jacobian matrices.
c this version is in single precision.
c
c lsodes solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c lsodes is a variant of the lsode package, and is intended for
c problems in which the jacobian matrix df/dy has an arbitrary
c sparse structure (when the problem is stiff).
c
c authors..      alan c. hindmarsh,
c                computing and mathematics research division, l-316
c                lawrence livermore national laboratory
c                livermore, ca 94550.
c
c and            andrew h. sherman
c                j. s. nolen and associates
c                houston, tx 77084
c-----------------------------------------------------------------------
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c
c 2.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
c     yale sparse matrix package.. i. the symmetric codes,
c     int. j. num. meth. eng., 18 (1982), pp. 1145-1151.
c
c 3.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
c     yale sparse matrix package.. ii. the nonsymmetric codes,
c     research report no. 114, dept. of computer sciences, yale
c     university, 1977.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsodes package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including optional communication, nonstandard options,
c and instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. next determine (or guess) whether or not the problem is stiff.
c stiffness occurs when the jacobian matrix df/dy has an eigenvalue
c whose real part is negative and large in magnitude, compared to the
c reciprocal of the t span of interest.  if the problem is nonstiff,
c use a method flag mf = 10.  if it is stiff, there are two standard
c for the method flag, mf = 121 and mf = 222.  in both cases, lsodes
c requires the jacobian matrix in some form, and it treats this matrix
c in general sparse form, with sparsity structure determined internally.
c (for options where the user supplies the sparsity structure, see
c the full description of mf below.)
c
c c. if the problem is stiff, you are encouraged to supply the jacobian
c directly (mf = 121), but if this is not feasible, lsodes will
c compute it internally by difference quotients (mf = 222).
c if you are supplying the jacobian, provide a subroutine of the form..
c               subroutine jac (neq, t, y, j, ian, jan, pdj)
c               dimension y(1), ian(1), jan(1), pdj(1)
c here neq, t, y, and j are input arguments, and the jac routine is to
c load the array pdj (of length neq) with the j-th column of df/dy.
c i.e., load pdj(i) with df(i)/dy(j) for all relevant values of i.
c the arguments ian and jan should be ignored for normal situations.
c lsodes will call the jac routine with j = 1,2,...,neq.
c only nonzero elements need be loaded.  usually, a crude approximation
c to df/dy, possibly with fewer nonzero elements, will suffice.
c
c d. write a main program which calls subroutine lsodes once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsodes.  on the first call to lsodes, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be roughly less (in magnitude) than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             20 + 16*neq            for mf = 10,
c             20 + (2 + 1./lenrat)*nnz + (11 + 9./lenrat)*neq
c                                    for mf = 121 or 222,
c          where..
c          nnz    = the number of nonzero elements in the sparse
c                   jacobian (if this is unknown, use an estimate), and
c          lenrat = the real to integer wordlength ratio (usually 1 in
c                   single precision and 2 in double precision).
c          in any case, the required size of rwork cannot generally
c          be predicted in advance if mf = 121 or 222, and the value
c          above is a rough estimate of a crude lower bound.  some
c          experimentation with this size may be necessary.
c          (when known, the correct required length is an optional
c          output, available in iwork(17).)
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least 30.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix (mf = 121).
c          if used, this name must be declared external in calling
c          program.  if not used, pass a dummy name.
c mf     = method flag.  standard values are..
c          10  for nonstiff (adams) method, no jacobian used.
c          121 for stiff (bdf) method, user-supplied sparse jacobian.
c          222 for stiff method, internally generated sparse jacobian.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c e. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsodes was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong mf).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of mf or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 means a fatal error return flag came from the sparse
c             solver cdrv by way of prjs or slss.  should never happen.
c          a return with istate = -1, -4, or -5 may result from using
c          an inappropriate sparsity structure, one that is quite
c          different from the initial structure.  consider calling
c          lsodes again with istate = 3 to force the structure to be
c          reevaluated.  see the full description of istate below.
c
c f. to continue the integration after a successful return, simply
c reset tout and call lsodes again.  no other parameters need be reset.
c
c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsodes.  the problem is from chemical
c kinetics, and consists of the following 12 rate equations..
c    dy1/dt  = -rk1*y1
c    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
c                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
c    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
c                + rk11*rk14*y4 + rk12*rk14*y6
c    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
c    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
c    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
c    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
c    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
c    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
c    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
c                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
c                - rk6*y10 - rk9*y10
c    dy11/dt = rk10*y8
c    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
c                - rk15*y2*y12 - rk17*y10*y12
c
c with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
c      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
c      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
c      rk15 = rk17 = 100.0.
c
c the t interval is from 0 to 1000, and the initial conditions
c are y1 = 1, y2 = y3 = ... = y12 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsodes, using mf = 121
c and printing results at t = .1, 1., 10., 100., 1000.  it uses
c itol = 1 and mixed relative/absolute tolerance controls.
c during the run and at the end, statistical quantities of interest
c are printed (see optional outputs in the full description below).
c
c     external fex, jex
c     dimension y(12), rwork(500), iwork(30)
c     data lrw/500/, liw/30/
c     neq = 12
c     do 10 i = 1,neq
c 10    y(i) = 0.0e0
c     y(1) = 1.0e0
c     t = 0.0e0
c     tout = 0.1e0
c     itol = 1
c     rtol = 1.0e-4
c     atol = 1.0e-6
c     itask = 1
c     istate = 1
c     iopt = 0
c     mf = 121
c     do 40 iout = 1,5
c       call lsodes (fex, neq, y, t, tout, itol, rtol, atol,
c    1     itask, istate, iopt, rwork, lrw, iwork, liw, jex, mf)
c       write(6,30)t,iwork(11),rwork(11),(y(i),i=1,neq)
c 30    format(//7h at t =,e11.3,4x,
c    1    12h no. steps =,i5,4x,12h last step =,e11.3/
c    2    13h  y array =  ,4e14.5/13x,4e14.5/13x,4e14.5)
c       if (istate .lt. 0) go to 80
c       tout = tout*10.0e0
c 40    continue
c     lenrw = iwork(17)
c     leniw = iwork(18)
c     nst = iwork(11)
c     nfe = iwork(12)
c     nje = iwork(13)
c     nlu = iwork(21)
c     nnz = iwork(19)
c     nnzlu = iwork(25) + iwork(26) + neq
c     write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
c 70  format(//22h required rwork size =,i4,15h   iwork size =,i4/
c    1   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
c    2   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
c    3   26h   no. of nonzeros in lu =,i5)
c     stop
c 80  write(6,90)istate
c 90  format(///22h error halt.. istate =,i3)
c     stop
c     end
c
c     subroutine fex (neq, t, y, ydot)
c     real t, y, ydot
c     real rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
c    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
c     dimension y(12), ydot(12)
c     data rk1/0.1e0/, rk2/10.0e0/, rk3/50.0e0/, rk4/2.5e0/, rk5/0.1e0/,
c    1   rk6/10.0e0/, rk7/50.0e0/, rk8/2.5e0/, rk9/50.0e0/, rk10/5.0e0/,
c    2   rk11/50.0e0/, rk12/50.0e0/, rk13/50.0e0/, rk14/30.0e0/,
c    3   rk15/100.0e0/, rk16/2.5e0/, rk17/100.0e0/, rk18/2.5e0/,
c    4   rk19/50.0e0/, rk20/50.0e0/
c     ydot(1)  = -rk1*y(1)
c     ydot(2)  = rk1*y(1) + rk11*rk14*y(4) + rk19*rk14*y(5)
c    1           - rk3*y(2)*y(3) - rk15*y(2)*y(12) - rk2*y(2)
c     ydot(3)  = rk2*y(2) - rk5*y(3) - rk3*y(2)*y(3) - rk7*y(10)*y(3)
c    1           + rk11*rk14*y(4) + rk12*rk14*y(6)
c     ydot(4)  = rk3*y(2)*y(3) - rk11*rk14*y(4) - rk4*y(4)
c     ydot(5)  = rk15*y(2)*y(12) - rk19*rk14*y(5) - rk16*y(5)
c     ydot(6)  = rk7*y(10)*y(3) - rk12*rk14*y(6) - rk8*y(6)
c     ydot(7)  = rk17*y(10)*y(12) - rk20*rk14*y(7) - rk18*y(7)
c     ydot(8)  = rk9*y(10) - rk13*rk14*y(8) - rk10*y(8)
c     ydot(9)  = rk4*y(4) + rk16*y(5) + rk8*y(6) + rk18*y(7)
c     ydot(10) = rk5*y(3) + rk12*rk14*y(6) + rk20*rk14*y(7)
c    1           + rk13*rk14*y(8) - rk7*y(10)*y(3) - rk17*y(10)*y(12)
c    2           - rk6*y(10) - rk9*y(10)
c     ydot(11) = rk10*y(8)
c     ydot(12) = rk6*y(10) + rk19*rk14*y(5) + rk20*rk14*y(7)
c    1           - rk15*y(2)*y(12) - rk17*y(10)*y(12)
c     return
c     end
c
c     subroutine jex (neq, t, y, j, ia, ja, pdj)
c     real t, y, pdj
c     real rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
c    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
c     dimension y(1), ia(1), ja(1), pdj(1)
c     data rk1/0.1e0/, rk2/10.0e0/, rk3/50.0e0/, rk4/2.5e0/, rk5/0.1e0/,
c    1   rk6/10.0e0/, rk7/50.0e0/, rk8/2.5e0/, rk9/50.0e0/, rk10/5.0e0/,
c    2   rk11/50.0e0/, rk12/50.0e0/, rk13/50.0e0/, rk14/30.0e0/,
c    3   rk15/100.0e0/, rk16/2.5e0/, rk17/100.0e0/, rk18/2.5e0/,
c    4   rk19/50.0e0/, rk20/50.0e0/
c     go to (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), j
c 1   pdj(1) = -rk1
c     pdj(2) = rk1
c     return
c 2   pdj(2) = -rk3*y(3) - rk15*y(12) - rk2
c     pdj(3) = rk2 - rk3*y(3)
c     pdj(4) = rk3*y(3)
c     pdj(5) = rk15*y(12)
c     pdj(12) = -rk15*y(12)
c     return
c 3   pdj(2) = -rk3*y(2)
c     pdj(3) = -rk5 - rk3*y(2) - rk7*y(10)
c     pdj(4) = rk3*y(2)
c     pdj(6) = rk7*y(10)
c     pdj(10) = rk5 - rk7*y(10)
c     return
c 4   pdj(2) = rk11*rk14
c     pdj(3) = rk11*rk14
c     pdj(4) = -rk11*rk14 - rk4
c     pdj(9) = rk4
c     return
c 5   pdj(2) = rk19*rk14
c     pdj(5) = -rk19*rk14 - rk16
c     pdj(9) = rk16
c     pdj(12) = rk19*rk14
c     return
c 6   pdj(3) = rk12*rk14
c     pdj(6) = -rk12*rk14 - rk8
c     pdj(9) = rk8
c     pdj(10) = rk12*rk14
c     return
c 7   pdj(7) = -rk20*rk14 - rk18
c     pdj(9) = rk18
c     pdj(10) = rk20*rk14
c     pdj(12) = rk20*rk14
c     return
c 8   pdj(8) = -rk13*rk14 - rk10
c     pdj(10) = rk13*rk14
c     pdj(11) = rk10
c 9   return
c 10  pdj(3) = -rk7*y(3)
c     pdj(6) = rk7*y(3)
c     pdj(7) = rk17*y(12)
c     pdj(8) = rk9
c     pdj(10) = -rk7*y(3) - rk17*y(12) - rk6 - rk9
c     pdj(12) = rk6 - rk17*y(12)
c 11  return
c 12  pdj(2) = -rk15*y(2)
c     pdj(5) = rk15*y(2)
c     pdj(7) = rk17*y(10)
c     pdj(10) = -rk17*y(10)
c     pdj(12) = -rk15*y(2) - rk17*y(10)
c     return
c     end
c
c the output of this program (on a cray-1 in single precision)
c is as follows..
c
c
c at t =  1.000e-01     no. steps =   12     last step =  1.515e-02
c  y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
c                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
c                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
c
c
c at t =  1.000e+00     no. steps =   33     last step =  7.880e-02
c  y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
c                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
c                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
c
c
c at t =  1.000e+01     no. steps =   48     last step =  1.239e+00
c  y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
c                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
c                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
c
c
c at t =  1.000e+02     no. steps =   91     last step =  3.764e+00
c  y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
c                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
c                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
c
c
c at t =  1.000e+03     no. steps =  111     last step =  4.156e+02
c  y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
c               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
c                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
c
c
c required rwork size = 442   iwork size =  30
c no. steps = 111   no. f-s = 142   no. j-s =   2   no. lu-s =  20
c no. of nonzeros in j =   44   no. of nonzeros in lu =   50
c-----------------------------------------------------------------------
c full description of user interface to lsodes.
c
c the user interface to lsodes consists of the following parts.
c
c i.   the call sequence to subroutine lsodes, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  descriptions of other routines in the lsodes package that may be
c      (optionally) called by the user.  these provide the ability to
c      alter error message handling, save and restore the internal
c      common, and obtain specified derivatives of the solution y(t).
c
c iii. descriptions of common blocks to be declared in overlay
c      or similar environments, or to be saved when doing an interrupt
c      of the problem and continued solution later.
c
c iv.  description of two routines in the lsodes package, either of
c      which the user may replace with his own version, if desired.
c      these relate to the measurement of errors.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
c and those used for both input and output are
c     y, t, istate.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine lsodes to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on the
c initial call for the problem, but not checked thereafter unless a
c change in input parameters is flagged by istate = 3 on input.
c
c the descriptions of the call arguments are as follows.
c
c f      = the name of the user-supplied subroutine defining the
c          ode system.  the system must be put in the first-order
c          form dy/dt = f(t,y), where f is a vector-valued function
c          of the scalar t and the vector y.  subroutine f is to
c          compute the function f.  it is to have the form
c               subroutine f (neq, t, y, ydot)
c               dimension y(1), ydot(1)
c          where neq, t, and y are input, and the array ydot = f(t,y)
c          is output.  y and ydot are arrays of length neq.
c          (in the dimension statement above, 1 is a dummy
c          dimension.. it can be replaced by any value.)
c          subroutine f should not alter y(1),...,y(neq).
c          f must be declared external in the calling program.
c
c          subroutine f may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in f) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y below.
c
c          if quantities computed in the f routine are needed
c          externally to lsodes, an extra call to f should be made
c          for this purpose, for consistent and accurate results.
c          if only the derivative dy/dt is needed, use intdy instead.
c
c neq    = the size of the ode system (number of first order
c          ordinary differential equations).  used only for input.
c          neq may be decreased, but not increased, during the problem.
c          if neq is decreased (with istate = 3 on input), the
c          remaining components of y should be left undisturbed, if
c          these are to be accessed in f and/or jac.
c
c          normally, neq is a scalar, and it is generally referred to
c          as a scalar in this user interface description.  however,
c          neq may be an array, with neq(1) set to the system size.
c          (the lsodes package accesses only neq(1).)  in either case,
c          this parameter is passed as the neq argument in all calls
c          to f and jac.  hence, if it is an array, locations
c          neq(2),... may be used to store other integer data and pass
c          it to f and/or jac.  subroutines f and/or jac must include
c          neq in a dimension statement in that case.
c
c y      = a real array for the vector of dependent variables, of
c          length neq or more.  used for both input and output on the
c          first call (istate = 1), and only for output on other calls.
c          on the first call, y must contain the vector of initial
c          values.  on output, y contains the computed solution vector,
c          evaluated at t.  if desired, the y array may be used
c          for other purposes between calls to the solver.
c
c          this array is passed as the y argument in all calls to
c          f and jac.  hence its length may exceed neq, and locations
c          y(neq+1),... may be used to store other real data and
c          pass it to f and/or jac.  (the lsodes package accesses only
c          y(1),...,y(neq).)
c
c t      = the independent variable.  on input, t is used only on the
c          first call, as the initial point of the integration.
c          on output, after each call, t is the value at which a
c          computed solution y is evaluated (usually the same as tout).
c          on an error return, t is the farthest point reached.
c
c tout   = the next value of t at which a computed solution is desired.
c          used only for input.
c
c          when starting the problem (istate = 1), tout may be equal
c          to t for one call, then should .ne. t for the next call.
c          for the initial t, an input value of tout .ne. t is used
c          in order to determine the direction of the integration
c          (i.e. the algebraic sign of the step sizes) and the rough
c          scale of the problem.  integration in either direction
c          (forward or backward in t) is permitted.
c
c          if itask = 2 or 5 (one-step modes), tout is ignored after
c          the first call (i.e. the first call with tout .ne. t).
c          otherwise, tout is required on every call.
c
c          if itask = 1, 3, or 4, the values of tout need not be
c          monotone, but a value of tout which backs up is limited
c          to the current internal t interval, whose endpoints are
c          tcur - hu and tcur (see optional outputs, below, for
c          tcur and hu).
c
c itol   = an indicator for the type of error control.  see
c          description below under atol.  used only for input.
c
c rtol   = a relative error tolerance parameter, either a scalar or
c          an array of length neq.  see description below under atol.
c          input only.
c
c atol   = an absolute error tolerance parameter, either a scalar or
c          an array of length neq.  input only.
c
c             the input parameters itol, rtol, and atol determine
c          the error control performed by the solver.  the solver will
c          control the vector e = (e(i)) of estimated local errors
c          in y, according to an inequality of the form
c                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
c          and the rms-norm (root-mean-square norm) here is
c          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
c          is a vector of weights which must always be positive, and
c          the values of rtol and atol should all be non-negative.
c          the following table gives the types (scalar/array) of
c          rtol and atol, and the corresponding form of ewt(i).
c
c             itol    rtol       atol          ewt(i)
c              1     scalar     scalar     rtol*abs(y(i)) + atol
c              2     scalar     array      rtol*abs(y(i)) + atol(i)
c              3     array      scalar     rtol(i)*abs(y(i)) + atol
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)
c
c          when either of these parameters is a scalar, it need not
c          be dimensioned in the user-s calling program.
c
c          if none of the above choices (with itol, rtol, and atol
c          fixed throughout the problem) is suitable, more general
c          error controls can be obtained by substituting
c          user-supplied routines for the setting of ewt and/or for
c          the norm calculation.  see part iv below.
c
c          if global errors are to be estimated by making a repeated
c          run on the same problem with smaller tolerances, then all
c          components of rtol and atol (i.e. of ewt) should be scaled
c          down uniformly.
c
c itask  = an index specifying the task to be performed.
c          input only.  itask has the following values and meanings.
c          1  means normal computation of output values of y(t) at
c             t = tout (by overshooting and interpolating).
c          2  means take one step only and return.
c          3  means stop at the first internal mesh point at or
c             beyond t = tout and return.
c          4  means normal computation of output values of y(t) at
c             t = tout but without overshooting t = tcrit.
c             tcrit must be input as rwork(1).  tcrit may be equal to
c             or beyond tout, but not behind it in the direction of
c             integration.  this option is useful if the problem
c             has a singularity at or beyond t = tcrit.
c          5  means take one step, without passing tcrit, and return.
c             tcrit must be input as rwork(1).
c
c          note..  if itask = 4 or 5 and the solver reaches tcrit
c          (within roundoff), it will return t = tcrit (exactly) to
c          indicate this (unless itask = 4 and tout comes before tcrit,
c          in which case answers at t = tout are returned first).
c
c istate = an index used for input and output to specify the
c          the state of the calculation.
c
c          on input, the values of istate are as follows.
c          1  means this is the first call for the problem
c             (initializations will be done).  see note below.
c          2  means this is not the first call, and the calculation
c             is to continue normally, with no change in any input
c             parameters except possibly tout and itask.
c             (if itol, rtol, and/or atol are changed between calls
c             with istate = 2, the new values will be used but not
c             tested for legality.)
c          3  means this is not the first call, and the
c             calculation is to continue normally, but with
c             a change in input parameters other than
c             tout and itask.  changes are allowed in
c             neq, itol, rtol, atol, iopt, lrw, liw, mf,
c             the conditional inputs ia and ja,
c             and any of the optional inputs except h0.
c             in particular, if miter = 1 or 2, a call with istate = 3
c             will cause the sparsity structure of the problem to be
c             recomputed (or reread from ia and ja if moss = 0).
c          note..  a preliminary call with tout = t is not counted
c          as a first call here, as no initialization or checking of
c          input is done.  (such a call is sometimes useful for the
c          purpose of outputting the initial conditions.)
c          thus the first call for which tout .ne. t requires
c          istate = 1 on input.
c
c          on output, istate has the following values and meanings.
c           1  means nothing was done, as tout was equal to t with
c              istate = 1 on input.  (however, an internal counter was
c              set to detect and prevent repeated calls of this type.)
c           2  means the integration was performed successfully.
c          -1  means an excessive amount of work (more than mxstep
c              steps) was done on this call, before completing the
c              requested task, but the integration was otherwise
c              successful as far as t.  (mxstep is an optional input
c              and is normally 500.)  to continue, the user may
c              simply reset istate to a value .gt. 1 and call again
c              (the excess work step counter will be reset to 0).
c              in addition, the user may increase mxstep to avoid
c              this error return (see below on optional inputs).
c          -2  means too much accuracy was requested for the precision
c              of the machine being used.  this was detected before
c              completing the requested task, but the integration
c              was successful as far as t.  to continue, the tolerance
c              parameters must be reset, and istate must be set
c              to 3.  the optional output tolsf may be used for this
c              purpose.  (note.. if this condition is detected before
c              taking any steps, then an illegal input return
c              (istate = -3) occurs instead.)
c          -3  means illegal input was detected, before taking any
c              integration steps.  see written message for details.
c              note..  if the solver detects an infinite loop of calls
c              to the solver with illegal input, it will cause
c              the run to stop.
c          -4  means there were repeated error test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              the problem may have a singularity, or the input
c              may be inappropriate.
c          -5  means there were repeated convergence test failures on
c              one attempted step, before completing the requested
c              task, but the integration was successful as far as t.
c              this may be caused by an inaccurate jacobian matrix,
c              if one is being used.
c          -6  means ewt(i) became zero for some i during the
c              integration.  pure relative error control (atol(i)=0.0)
c              was requested on a variable which has now vanished.
c              the integration was successful as far as t.
c          -7  means a fatal error return flag came from the sparse
c              solver cdrv by way of prjs or slss (numerical
c              factorization or backsolve).  this should never happen.
c              the integration was successful as far as t.
c
c          note.. an error return with istate = -1, -4, or -5 and with
c          miter = 1 or 2 may mean that the sparsity structure of the
c          problem has changed significantly since it was last
c          determined (or input).  in that case, one can attempt to
c          complete the integration by setting istate = 3 on the next
c          call, so that a new structure determination is done.
c
c          note..  since the normal output value of istate is 2,
c          it does not need to be reset for normal continuation.
c          also, since a negative input value of istate will be
c          regarded as illegal, a negative output value requires the
c          user to change it, and possibly other inputs, before
c          calling the solver again.
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used on this call.  input only.
c          the optional inputs are listed separately below.
c          iopt = 0 means no optional inputs are being used.
c                   default values will be used in all cases.
c          iopt = 1 means one or more optional inputs are being used.
c
c rwork  = a work array used for a mixture of real (single precision)
c          and integer work space.
c          the length of rwork (in real words) must be at least
c             20 + nyh*(maxord + 1) + 3*neq + lwm    where
c          nyh    = the initial value of neq,
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
c                   smaller value is given as an optional input),
c          lwm = 0                                    if miter = 0,
c          lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat   if miter = 1,
c          lwm = 2*nnz + 2*neq + (nnz+10*neq)/lenrat  if miter = 2,
c          lwm = neq + 2                              if miter = 3.
c          in the above formulas,
c          nnz    = number of nonzero elements in the jacobian matrix.
c          lenrat = the real to integer wordlength ratio (usually 1 in
c                   single precision and 2 in double precision).
c          (see the mf description for meth and miter.)
c          thus if maxord has its default value and neq is constant,
c          the minimum length of rwork is..
c             20 + 16*neq        for mf = 10,
c             20 + 16*neq + lwm  for mf = 11, 111, 211, 12, 112, 212,
c             22 + 17*neq        for mf = 13,
c             20 +  9*neq        for mf = 20,
c             20 +  9*neq + lwm  for mf = 21, 121, 221, 22, 122, 222,
c             22 + 10*neq        for mf = 23.
c          if miter = 1 or 2, the above formula for lwm is only a
c          crude lower bound.  the required length of rwork cannot
c          be readily predicted in general, as it depends on the
c          sparsity structure of the problem.  some experimentation
c          may be necessary.
c
c          the first 20 words of rwork are reserved for conditional
c          and optional inputs and optional outputs.
c
c          the following word in rwork is a conditional input..
c            rwork(1) = tcrit = critical value of t which the solver
c                       is not to overshoot.  required if itask is
c                       4 or 5, and ignored otherwise.  (see itask.)
c
c lrw    = the length of the array rwork, as declared by the user.
c          (this will be checked by the solver.)
c
c iwork  = an integer work array.  the length of iwork must be at least
c             31 + neq + nnz   if moss = 0 and miter = 1 or 2, or
c             30               otherwise.
c          (nnz is the number of nonzero elements in df/dy.)
c
c          in lsodes, iwork is used only for conditional and
c          optional inputs and optional outputs.
c
c          the following two blocks of words in iwork are conditional
c          inputs, required if moss = 0 and miter = 1 or 2, but not
c          otherwise (see the description of mf for moss).
c            iwork(30+j) = ia(j)     (j=1,...,neq+1)
c            iwork(31+neq+k) = ja(k) (k=1,...,nnz)
c          the two arrays ia and ja describe the sparsity structure
c          to be assumed for the jacobian matrix.  ja contains the row
c          indices where nonzero elements occur, reading in columnwise
c          order, and ia contains the starting locations in ja of the
c          descriptions of columns 1,...,neq, in that order, with
c          ia(1) = 1.  thus, for each column index j = 1,...,neq, the
c          values of the row index i in column j where a nonzero
c          element may occur are given by
c            i = ja(k),  where   ia(j) .le. k .lt. ia(j+1).
c          if nnz is the total number of nonzero locations assumed,
c          then the length of the ja array is nnz, and ia(neq+1) must
c          be nnz + 1.  duplicate entries are not allowed.
c
c liw    = the length of the array iwork, as declared by the user.
c          (this will be checked by the solver.)
c
c note..  the work arrays must not be altered between calls to lsodes
c for the same problem, except possibly for the conditional and
c optional inputs, and except for the last 3*neq words of rwork.
c the latter space is used for internal scratch space, and so is
c available for use by the user outside lsodes between calls, if
c desired (but not for use by f or jac).
c
c jac    = name of user-supplied routine (miter = 1 or moss = 1) to
c          compute the jacobian matrix, df/dy, as a function of
c          the scalar t and the vector y.  it is to have the form
c               subroutine jac (neq, t, y, j, ian, jan, pdj)
c               dimension y(1), ian(1), jan(1), pdj(1)
c          where neq, t, y, j, ian, and jan are input, and the array
c          pdj, of length neq, is to be loaded with column j
c          of the jacobian on output.  thus df(i)/dy(j) is to be
c          loaded into pdj(i) for all relevant values of i.
c          here t and y have the same meaning as in subroutine f,
c          and j is a column index (1 to neq).  ian and jan are
c          undefined in calls to jac for structure determination
c          (moss = 1).  otherwise, ian and jan are structure
c          descriptors, as defined under optional outputs below, and
c          so can be used to determine the relevant row indices i, if
c          desired.  (in the dimension statement above, 1 is a
c          dummy dimension.. it can be replaced by any value.)
c               jac need not provide df/dy exactly.  a crude
c          approximation (possibly with greater sparsity) will do.
c               in any case, pdj is preset to zero by the solver,
c          so that only the nonzero elements need be loaded by jac.
c          calls to jac are made with j = 1,...,neq, in that order, and
c          each such set of calls is preceded by a call to f with the
c          same arguments neq, t, and y.  thus to gain some efficiency,
c          intermediate quantities shared by both calculations may be
c          saved in a user common block by f and not recomputed by jac,
c          if desired.  jac must not alter its input arguments.
c          jac must be declared external in the calling program.
c               subroutine jac may access user-defined quantities in
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array
c          (dimensioned in jac) and/or y has length exceeding neq(1).
c          see the descriptions of neq and y above.
c
c mf     = the method flag.  used only for input.
c          mf has three decimal digits-- moss, meth, miter--
c             mf = 100*moss + 10*meth + miter.
c          moss indicates the method to be used to obtain the sparsity
c          structure of the jacobian matrix if miter = 1 or 2..
c            moss = 0 means the user has supplied ia and ja
c                     (see descriptions under iwork above).
c            moss = 1 means the user has supplied jac (see below)
c                     and the structure will be obtained from neq
c                     initial calls to jac.
c            moss = 2 means the structure will be obtained from neq+1
c                     initial calls to f.
c          meth indicates the basic linear multistep method..
c            meth = 1 means the implicit adams method.
c            meth = 2 means the method based on backward
c                     differentiation formulas (bdf-s).
c          miter indicates the corrector iteration method..
c            miter = 0 means functional iteration (no jacobian matrix
c                      is involved).
c            miter = 1 means chord iteration with a user-supplied
c                      sparse jacobian, given by subroutine jac.
c            miter = 2 means chord iteration with an internally
c                      generated (difference quotient) sparse jacobian
c                      (using ngp extra calls to f per df/dy value,
c                      where ngp is an optional output described below.)
c            miter = 3 means chord iteration with an internally
c                      generated diagonal jacobian approximation.
c                      (using 1 extra call to f per df/dy evaluation).
c          if miter = 1 or moss = 1, the user must supply a subroutine
c          jac (the name is arbitrary) as described above under jac.
c          otherwise, a dummy argument can be used.
c
c          the standard choices for mf are..
c            mf = 10  for a nonstiff problem,
c            mf = 21 or 22 for a stiff problem with ia/ja supplied
c                     (21 if jac is supplied, 22 if not),
c            mf = 121 for a stiff problem with jac supplied,
c                     but not ia/ja,
c            mf = 222 for a stiff problem with neither ia/ja nor
c                     jac supplied.
c          the sparseness structure can be changed during the
c          problem by making a call to lsodes with istate = 3.
c-----------------------------------------------------------------------
c optional inputs.
c
c the following is a list of the optional inputs provided for in the
c call sequence.  (see also part ii.)  for each such input variable,
c this table lists its name as used in this documentation, its
c location in the call sequence, its meaning, and the default value.
c the use of any of these inputs requires iopt = 1, and in that
c case all of these inputs are examined.  a value of zero for any
c of these optional inputs will cause the default value to be used.
c thus to use a subset of the optional inputs, simply preload
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
c then set those of interest to nonzero values.
c
c name    location      meaning and default value
c
c h0      rwork(5)  the step size to be attempted on the first step.
c                   the default value is determined by the solver.
c
c hmax    rwork(6)  the maximum absolute step size allowed.
c                   the default value is infinite.
c
c hmin    rwork(7)  the minimum absolute step size allowed.
c                   the default value is 0.  (this lower bound is not
c                   enforced on the final step before reaching tcrit
c                   when itask = 4 or 5.)
c
c seth    rwork(8)  the element threshhold for sparsity determination
c                   when moss = 1 or 2.  if the absolute value of
c                   an estimated jacobian element is .le. seth, it
c                   will be assumed to be absent in the structure.
c                   the default value of seth is 0.
c
c maxord  iwork(5)  the maximum order to be allowed.  the default
c                   value is 12 if meth = 1, and 5 if meth = 2.
c                   if maxord exceeds the default value, it will
c                   be reduced to the default value.
c                   if maxord is changed during the problem, it may
c                   cause the current order to be reduced.
c
c mxstep  iwork(6)  maximum number of (internally defined) steps
c                   allowed during one call to the solver.
c                   the default value is 500.
c
c mxhnil  iwork(7)  maximum number of messages printed (per problem)
c                   warning that t + h = t on a step (h = step size).
c                   this must be positive to result in a non-default
c                   value.  the default value is 10.
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from lsodes, the variables listed
c below are quantities related to the performance of lsodes
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on any successful return from lsodes, and on any return with
c istate = -1, -2, -4, -5, or -6.  on an illegal input return
c (istate = -3), they will be unchanged from their existing values
c (if any), except possibly for tolsf, lenrw, and leniw.
c on any error return, outputs relevant to the error will be defined,
c as noted below.
c
c name    location      meaning
c
c hu      rwork(11) the step size in t last used (successfully).
c
c hcur    rwork(12) the step size to be attempted on the next step.
c
c tcur    rwork(13) the current value of the independent variable
c                   which the solver has actually reached, i.e. the
c                   current internal mesh point in t.  on output, tcur
c                   will always be at least as far as the argument
c                   t, but may be farther (if interpolation was done).
c
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
c                   computed when a request for too much accuracy was
c                   detected (istate = -3 if detected at the start of
c                   the problem, istate = -2 otherwise).  if itol is
c                   left unaltered but rtol and atol are uniformly
c                   scaled up by a factor of tolsf for the next call,
c                   then the solver is deemed likely to succeed.
c                   (the user may also ignore tolsf and alter the
c                   tolerance parameters in any other way appropriate.)
c
c nst     iwork(11) the number of steps taken for the problem so far.
c
c nfe     iwork(12) the number of f evaluations for the problem so far,
c                   excluding those for structure determination
c                   (moss = 2).
c
c nje     iwork(13) the number of jacobian evaluations for the problem
c                   so far, excluding those for structure determination
c                   (moss = 1).
c
c nqu     iwork(14) the method order last used (successfully).
c
c nqcur   iwork(15) the order to be attempted on the next step.
c
c imxer   iwork(16) the index of the component of largest magnitude in
c                   the weighted local error vector ( e(i)/ewt(i) ),
c                   on an error return with istate = -4 or -5.
c
c lenrw   iwork(17) the length of rwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c leniw   iwork(18) the length of iwork actually required.
c                   this is defined on normal returns and on an illegal
c                   input return for insufficient storage.
c
c nnz     iwork(19) the number of nonzero elements in the jacobian
c                   matrix, including the diagonal (miter = 1 or 2).
c                   (this may differ from that given by ia(neq+1)-1
c                   if moss = 0, because of added diagonal entries.)
c
c ngp     iwork(20) the number of groups of column indices, used in
c                   difference quotient jacobian aproximations if
c                   miter = 2.  this is also the number of extra f
c                   evaluations needed for each jacobian evaluation.
c
c nlu     iwork(21) the number of sparse lu decompositions for the
c                   problem so far.
c
c lyh     iwork(22) the base address in rwork of the history array yh,
c                   described below in this list.
c
c ipian   iwork(23) the base address of the structure descriptor array
c                   ian, described below in this list.
c
c ipjan   iwork(24) the base address of the structure descriptor array
c                   jan, described below in this list.
c
c nzl     iwork(25) the number of nonzero elements in the strict lower
c                   triangle of the lu factorization used in the chord
c                   iteration (miter = 1 or 2).
c
c nzu     iwork(26) the number of nonzero elements in the strict upper
c                   triangle of the lu factorization used in the chord
c                   iteration (miter = 1 or 2).
c                   the total number of nonzeros in the factorization
c                   is therefore nzl + nzu + neq.
c
c the following four arrays are segments of the rwork array which
c may also be of interest to the user as optional outputs.
c for each array, the table below gives its internal name,
c its base address, and its description.
c for yh and acor, the base addresses are in rwork (a real array).
c the integer arrays ian and jan are to be obtained by declaring an
c integer array iwk and identifying iwk(1) with rwork(21), using either
c an equivalence statement or a subroutine call.  then the base
c addresses ipian (of ian) and ipjan (of jan) in iwk are to be obtained
c as optional outputs iwork(23) and iwork(24), respectively.
c thus ian(1) is iwk(ipian), etc.
c
c name    base address      description
c
c ian    ipian (in iwk)  structure descriptor array of size neq + 1.
c jan    ipjan (in iwk)  structure descriptor array of size nnz.
c         (see above)    ian and jan together describe the sparsity
c                        structure of the jacobian matrix, as used by
c                        lsodes when miter = 1 or 2.
c                        jan contains the row indices of the nonzero
c                        locations, reading in columnwise order, and
c                        ian contains the starting locations in jan of
c                        the descriptions of columns 1,...,neq, in
c                        that order, with ian(1) = 1.  thus for each
c                        j = 1,...,neq, the row indices i of the
c                        nonzero locations in column j are
c                        i = jan(k),  ian(j) .le. k .lt. ian(j+1).
c                        note that ian(neq+1) = nnz + 1.
c                        (if moss = 0, ian/jan may differ from the
c                        input ia/ja because of a different ordering
c                        in each column, and added diagonal entries.)
c
c yh      lyh            the nordsieck history array, of size nyh by
c          (optional     (nqcur + 1), where nyh is the initial value
c          output)       of neq.  for j = 0,1,...,nqcur, column j+1
c                        of yh contains hcur**j/factorial(j) times
c                        the j-th derivative of the interpolating
c                        polynomial currently representing the solution,
c                        evaluated at t = tcur.  the base address lyh
c                        is another optional output, listed above.
c
c acor     lenrw-neq+1   array of size neq used for the accumulated
c                        corrections on each step, scaled on output
c                        to represent the estimated local error in y
c                        on the last step.  this is the vector e in
c                        the description of the error control.  it is
c                        defined only on a successful return from
c                        lsodes.
c
c-----------------------------------------------------------------------
c part ii.  other routines callable.
c
c the following are optional calls which the user may make to
c gain additional capabilities in conjunction with lsodes.
c (the routines xsetun and xsetf are designed to conform to the
c slatec error handling package.)
c
c     form of call                  function
c   call xsetun(lun)          set the logical unit number, lun, for
c                             output of messages from lsodes, if
c                             the default is not desired.
c                             the default value of lun is 6.
c
c   call xsetf(mflag)         set a flag to control the printing of
c                             messages by lsodes.
c                             mflag = 0 means do not print. (danger..
c                             this risks losing valuable information.)
c                             mflag = 1 means print (the default).
c
c                             either of the above calls may be made at
c                             any time and will take effect immediately.
c
c   call srcms(rsav,isav,job) saves and restores the contents of
c                             the internal common blocks used by
c                             lsodes (see part iii below).
c                             rsav must be a real array of length 224
c                             or more, and isav must be an integer
c                             array of length 75 or more.
c                             job=1 means save common into rsav/isav.
c                             job=2 means restore common from rsav/isav.
c                                srcms is useful if one is
c                             interrupting a run and restarting
c                             later, or alternating between two or
c                             more problems solved with lsodes.
c
c   call intdy(,,,,,)         provide derivatives of y, of various
c        (see below)          orders, at a specified point t, if
c                             desired.  it may be called only after
c                             a successful return from lsodes.
c
c the detailed instructions for using intdy are as follows.
c the form of the call is..
c
c   lyh = iwork(22)
c   call intdy (t, k, rwork(lyh), nyh, dky, iflag)
c
c the input parameters are..
c
c t         = value of independent variable where answers are desired
c             (normally the same as the t last returned by lsodes).
c             for valid results, t must lie between tcur - hu and tcur.
c             (see optional outputs for tcur and hu.)
c k         = integer order of the derivative desired.  k must satisfy
c             0 .le. k .le. nqcur, where nqcur is the current order
c             (see optional outputs).  the capability corresponding
c             to k = 0, i.e. computing y(t), is already provided
c             by lsodes directly.  since nqcur .ge. 1, the first
c             derivative dy/dt is always available with intdy.
c lyh       = the base address of the history array yh, obtained
c             as an optional output as shown above.
c nyh       = column length of yh, equal to the initial value of neq.
c
c the output parameters are..
c
c dky       = a real array of length neq containing the computed value
c             of the k-th derivative of y(t).
c iflag     = integer flag, returned as 0 if k and t were legal,
c             -1 if k was illegal, and -2 if t was illegal.
c             on an error return, a message is also written.
c-----------------------------------------------------------------------
c part iii.  common blocks.
c
c if lsodes is to be used in an overlay situation, the user
c must declare, in the primary overlay, the variables in..
c   (1) the call sequence to lsodes,
c   (2) the three internal common blocks
c         /ls0001/  of length  257  (218 single precision words
c                         followed by 39 integer words),
c         /lss001/  of length  40    ( 6 single precision words
c                         followed by 34 integer words),
c         /eh0001/  of length  2 (integer words).
c
c if lsodes is used on a system in which the contents of internal
c common blocks are not preserved between calls, the user should
c declare the above three common blocks in his main program to insure
c that their contents are preserved.
c
c if the solution of a given problem by lsodes is to be interrupted
c and then later continued, such as when restarting an interrupted run
c or alternating between two or more problems, the user should save,
c following the return from the last lsodes call prior to the
c interruption, the contents of the call sequence variables and the
c internal common blocks, and later restore these values before the
c next lsodes call for that problem.  to save and restore the common
c blocks, use subroutine srcms (see part ii above).
c
c note.. in this version of lsodes, there are two data statements,
c in subroutines lsodes and xerrwv, which load variables into these
c labeled common blocks.  on some systems, it may be necessary to
c move these to a separate block data subprogram.
c
c-----------------------------------------------------------------------
c part iv.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the lsodes package which
c relate to the measurement of errors.  either routine can be
c replaced by a user-supplied version, if desired.  however, since such
c a replacement may have a major impact on performance, it should be
c done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) ewset.
c the following subroutine is called just before each internal
c integration step, and sets the array of error weights, ewt, as
c described under itol/rtol/atol above..
c     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
c where neq, itol, rtol, and atol are as in the lsodes call sequence,
c ycur contains the current dependent variable vector, and
c ewt is the array of weights set by ewset.
c
c if the user supplies this subroutine, it must return in ewt(i)
c (i = 1,...,neq) a positive quantity suitable for comparing errors
c in y(i) to.  the ewt array returned by ewset is passed to the
c vnorm routine (see below), and also used by lsodes in the computation
c of the optional output imxer, the diagonal jacobian approximation,
c and the increments for difference quotient jacobians.
c
c in the user-supplied version of ewset, it may be desirable to use
c the current values of derivatives of y.  derivatives up to order nq
c are available from the history array yh, described above under
c optional outputs.  in ewset, yh is identical to the ycur array,
c extended to nq + 1 columns with a column length of nyh and scale
c factors of h**j/factorial(j).  on the first call for the problem,
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.
c the quantities nq, nyh, h, and nst can be obtained by including
c in ewset the statements..
c     common /ls0001/ rls(218),ils(39)
c     nq = ils(35)
c     nyh = ils(14)
c     nst = ils(36)
c     h = rls(212)
c thus, for example, the current value of dy/dt can be obtained as
c ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
c unnecessary when nst = 0).
c
c (b) vnorm.
c the following is a real function routine which computes the weighted
c root-mean-square norm of a vector v..
c     d = vnorm (n, v, w)
c where..
c   n = the length of the vector,
c   v = real array of length n containing the vector,
c   w = real array of length n containing weights,
c   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
c vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
c ewt is as set by subroutine ewset.
c
c if the user supplies this function, it should return a non-negative
c value of vnorm suitable for use in the error control in lsodes.
c none of the arguments should be altered by vnorm.
c for example, a user-supplied vnorm routine might..
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
c   -ignore some components of v in the norm, with the effect of
c    suppressing the error control on those components of y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c other routines in the lsodes package.
c
c in addition to subroutine lsodes, the lsodes package includes the
c following subroutines and function routines..
c  iprep    acts as an iterface between lsodes and prep, and also does
c           adjusting of work space pointers and work arrays.
c  prep     is called by iprep to compute sparsity and do sparse matrix
c           preprocessing if miter = 1 or 2.
c  jgroup   is called by prep to compute groups of jacobian column
c           indices for use when miter = 2.
c  adjlr    adjusts the length of required sparse matrix work space.
c           it is called by prep.
c  cntnzu   is called by prep and counts the nonzero elements in the
c           strict upper triangle of j + j-transpose, where j = df/dy.
c  intdy    computes an interpolated value of the y vector at t = tout.
c  stode    is the core integrator, which does one step of the
c           integration and the associated error control.
c  cfode    sets all method coefficients and test constants.
c  prjs     computes and preprocesses the jacobian matrix j = df/dy
c           and the newton iteration matrix p = i - h*l0*j.
c  slss     manages solution of linear system in chord iteration.
c  ewset    sets the error weight vector ewt before each step.
c  vnorm    computes the weighted r.m.s. norm of a vector.
c  srcms    is a user-callable routine to save and restore
c           the contents of the internal common blocks.
c  odrv     constructs a reordering of the rows and columns of
c           a matrix by the minimum degree algorithm.  odrv is a
c           driver routine which calls subroutines md, mdi, mdm,
c           mdp, mdu, and sro.  see ref. 2 for details.  (the odrv
c           module has been modified since ref. 2, however.)
c  cdrv     performs reordering, symbolic factorization, numerical
c           factorization, or linear system solution operations,
c           depending on a path argument ipath.  cdrv is a
c           driver routine which calls subroutines nroc, nsfc,
c           nnfc, nnsc, and nntc.  see ref. 3 for details.
c           lsodes uses cdrv to solve linear systems in which the
c           coefficient matrix is  p = i - con*j, where i is the
c           identity, con is a scalar, and j is an approximation to
c           the jacobian df/dy.  because cdrv deals with rowwise
c           sparsity descriptions, cdrv works with p-transpose, not p.
c  r1mach   computes the unit roundoff in a machine-independent manner.
c  xerrwv, xsetun, and xsetf   handle the printing of all error
c           messages and warnings.  xerrwv is machine-dependent.
c note..  vnorm and r1mach are function routines.
c all the others are subroutines.
c
c the intrinsic and external routines used by lsodes are..
c abs, amax1, amin1, float, max0, min0, mod, sign, sqrt, and write.
c
c-----------------------------------------------------------------------
c the following card is for optimized compilation on lll compilers.
clll. optimize
c-----------------------------------------------------------------------
      external prjs, slss
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem,
     1   j, kgo, lenrat, lenyht, leniw, lenrw, lf0, lia, lja,
     2   lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth
      real atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   r1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following two internal common blocks contain
c (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of each block is as follows..  all real variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsodes first,
c then those local to subroutine stode or subroutine prjs
c (no other routines have own variables), and finally those used
c for communication.  the block ls0001 is declared in subroutines
c lsodes, iprep, prep, intdy, stode, prjs, and slss.  the block lss001
c is declared in subroutines lsodes, iprep, prep, prjs, and slss.
c groups of variables are replaced by dummy arrays in the common
c declarations in routines where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c
      data mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c      data illin/0/, ntrep/0/
c-----------------------------------------------------------------------
c in the data statement below, set lenrat equal to the ratio of
c the wordlength for a real number to that for an integer.  usually,
c lenrat = 1 for single precision and 2 for double precision.  if the
c true ratio is not an integer, use the next smaller integer (.ge. 1).
c-----------------------------------------------------------------------
      data lenrat/1/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c if istate = 1, the final setting of work space pointers, the matrix
c preprocessing, and other initializations are done in block c.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      if (moss .lt. 0 .or. moss .gt. 2) go to 608
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 3) go to 608
      if (miter .eq. 0 .or. miter .eq. 3) moss = 0
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0e0
      hmxi = 0.0e0
      hmin = 0.0e0
      seth = 0.0e0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0e0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0e0) go to 615
      hmxi = 0.0e0
      if (hmax .gt. 0.0e0) hmxi = 1.0e0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0e0) go to 616
      seth = rwork(8)
      if (seth .lt. 0.0e0) go to 609
c check rtol and atol for legality. ------------------------------------
 60   rtoli = rtol(1)
      atoli = atol(1)
      do 65 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0e0) go to 619
        if (atoli .lt. 0.0e0) go to 620
 65     continue
c-----------------------------------------------------------------------
c compute required work array lengths, as far as possible, and test
c these against lrw and liw.  then set tentative pointers for work
c arrays.  pointers to rwork/iwork segments are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  wm, yh, savf, ewt, acor.
c if miter = 1 or 2, the required length of the matrix work space wm
c is not yet known, and so a crude minimum value is used for the
c initial tests of lrw and liw, and yh is temporarily stored as far
c to the right in rwork as possible, to leave the maximum amount
c of space for wm for matrix preprocessing.  thus if miter = 1 or 2
c and moss .ne. 2, some of the segments of rwork are temporarily
c omitted, as they are not needed in the preprocessing.  these
c omitted segments are.. acor if istate = 1, ewt and acor if istate = 3
c and moss = 1, and savf, ewt, and acor if istate = 3 and moss = 0.
c-----------------------------------------------------------------------
      lrat = lenrat
      if (istate .eq. 1) nyh = n
      lwmin = 0
      if (miter .eq. 1) lwmin = 4*n + 10*n/lrat
      if (miter .eq. 2) lwmin = 4*n + 11*n/lrat
      if (miter .eq. 3) lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + n + 1
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
      lia = 31
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)
     1   leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      if (leniw .gt. liw) go to 618
      lja = lia + n + 1
      lia = min0(lia,liw)
      lja = min0(lja,liw)
      lwm = 21
      if (istate .eq. 1) nq = 1
      ncolm = min0(nq+1,maxord+2)
      lenyhm = ncolm*nyh
      lenyht = lenyh
      if (miter .eq. 1 .or. miter .eq. 2) lenyht = lenyhm
      imul = 2
      if (istate .eq. 3) imul = moss
      if (moss .eq. 2) imul = 3
      lrtem = lenyht + imul*n
      lwtem = lwmin
      if (miter .eq. 1 .or. miter .eq. 2) lwtem = lrw - 20 - lrtem
      lenwk = lwtem
      lyhn = lwm + lwtem
      lsavf = lyhn + lenyht
      lewt = lsavf + n
      lacor = lewt + n
      istatc = istate
      if (istate .eq. 1) go to 100
c-----------------------------------------------------------------------
c istate = 3.  move yh to its new location.
c note that only the part of yh needed for the next step, namely
c min(nq+1,maxord+2) columns, is actually moved.
c a temporary error weight array ewt is loaded if moss = 2.
c sparse matrix processing is done in iprep/prep if miter = 1 or 2.
c if maxord was reduced below nq, then the pointers are finally set
c so that savf is identical to yh(*,maxord+2).
c-----------------------------------------------------------------------
      lyhd = lyh - lyhn
      imax = lyhn - 1 + lenyhm
c move yh.  branch for move right, no move, or move left. --------------
      if (lyhd) 70,80,74
 70   do 72 i = lyhn,imax
        j = imax + lyhn - i
 72     rwork(j) = rwork(j+lyhd)
      go to 80
 74   do 76 i = lyhn,imax
 76     rwork(i) = rwork(i+lyhd)
 80   lyh = lyhn
      iwork(22) = lyh
      if (miter .eq. 0 .or. miter .eq. 3) go to 92
      if (moss .ne. 2) go to 85
c temporarily load ewt if miter = 1 or 2 and moss = 2. -----------------
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 82 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 82     rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 85   continue
c iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lsavf = min0(lsavf,lrw)
      lewt = min0(lewt,lrw)
      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (90, 628, 629, 630, 631, 632, 633), ipgo
 90   iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617
c set flag to signal parameter changes to stode. -----------------------
 92   jstart = -1
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0e0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c the sparse matrix preprocessing (miter = 1 or 2), and the
c calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  continue
      lyh = lyhn
      iwork(22) = lyh
      tn = t
      nst = 0
      h = 1.0e0
      nnz = 0
      ngp = 0
      nzl = 0
      nzu = 0
c load the initial value vector in yh. ---------------------------------
      do 105 i = 1,n
 105    rwork(i+lyh-1) = y(i)
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 110 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 110    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
      if (miter .eq. 0 .or. miter .eq. 3) go to 120
c iprep and prep do sparse matrix preprocessing if miter = 1 or 2. -----
      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (115, 628, 629, 630, 631, 632, 633), ipgo
 115  iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617
c check tcrit for legality (itask = 4 or 5). ---------------------------
 120  continue
      if (itask .ne. 4 .and. itask .ne. 5) go to 125
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0e0) go to 625
      if (h0 .ne. 0.0e0 .and. (t + h0 - tcrit)*h0 .gt. 0.0e0)
     1   h0 = tcrit - t
c initialize all remaining parameters. ---------------------------------
 125  uround = r1mach(4)
      jstart = 0
      if (miter .ne. 0) rwork(lwm) = sqrt(uround)
      msbj = 50
      nslj = 0
      ccmxj = 0.2e0
      psmall = 1000.0e0*uround
      rbig = 0.01e0/psmall
      nhnil = 0
      nje = 0
      nlu = 0
      nslast = 0
      hu = 0.0e0
      nqu = 0
      ccmax = 0.3e0
      maxcor = 3
      msbp = 20
      mxncf = 10
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      lf0 = lyh + nyh
      if (h0 .ne. 0.0e0) go to 180
      tdist = abs(tout - t)
      w0 = amax1(abs(t),abs(tout))
      if (tdist .lt. 2.0e0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = amax1(tol,rtol(i))
 140  if (tol .gt. 0.0e0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0e0) tol = amax1(tol,atoli/ayi)
 150    continue
 160  tol = amax1(tol,100.0e0*uround)
      tol = amin1(tol,0.001e0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = amin1(h0,tdist)
      h0 = sign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0e0 + 100.0e0*uround)
      if ((tp - tout)*h .gt. 0.0e0) go to 623
      if ((tn - tout)*h .lt. 0.0e0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
      if ((tcrit - tout)*h .lt. 0.0e0) go to 625
      if ((tn - tout)*h .lt. 0.0e0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
 245  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 510
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0e0) go to 280
      tolsf = tolsf*2.0e0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv(50hlsodes-- warning..internal t (=r1) and h (=r2) are,
     1   50, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h      such that in the machine, t + h = t on the next step  ,
     1   60, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      (h = step size). solver will continue anyway,
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv(50hlsodes-- above warning has been issued i1 times.  ,
     1   50, 102, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      it will not be issued again for this problem,
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
 290  continue
c-----------------------------------------------------------------------
c    call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,wm,f,jac,prjs,slss)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm),
     2   f, jac, prjs, slss)
      kgo = 1 - kflag
      go to (300, 530, 540, 550), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0e0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0e0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsodes.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  60hlsodes-- repeated calls with istate = 1 and tout = t (=r1)  ,
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0e0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xerrwv(50hlsodes-- at current t (=r1), mxstep (=i1) steps   ,
     1   50, 201, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      taken on this call before reaching tout     ,
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0e0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv(50hlsodes-- at t (=r1), ewt(i1) has become r2 .le. 0.,
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv(50hlsodes-- at t (=r1), too much accuracy requested  ,
     1   50, 203, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      for precision of machine..  see tolsf (=r2) ,
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv(50hlsodes-- at t(=r1) and step size h(=r2), the error,
     1   50, 204, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      test failed repeatedly or with abs(h) = hmin,
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), the    ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      corrector convergence failed repeatedly     ,
     1   50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(30h      or with abs(h) = hmin   ,
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 560
c kflag = -3.  fatal error flag returned by prjs or slss (cdrv). -------
 550  call xerrwv(50hlsodes-- at t (=r1) and step size h (=r2), a fatal,
     1   50, 207, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(50h      error flag was returned by cdrv (by way of  ,
     1   50, 207, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(30h      subroutine prjs or slss),
     1   30, 207, 0, 0, 0, 0, 2, tn, h)
      istate = -7
      go to 580
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0e0
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv(30hlsodes-- istate (=i1) illegal ,
     1   30, 1, 0, 1, istate, 0, 0, 0.0e0, 0.0e0)
      go to 700
 602  call xerrwv(30hlsodes-- itask (=i1) illegal  ,
     1   30, 2, 0, 1, itask, 0, 0, 0.0e0, 0.0e0)
      go to 700
 603  call xerrwv(50hlsodes-- istate .gt. 1 but lsodes not initialized ,
     1   50, 3, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      go to 700
 604  call xerrwv(30hlsodes-- neq (=i1) .lt. 1     ,
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0e0, 0.0e0)
      go to 700
 605  call xerrwv(50hlsodes-- istate = 3 and neq increased (i1 to i2)  ,
     1   50, 5, 0, 2, n, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 606  call xerrwv(30hlsodes-- itol (=i1) illegal   ,
     1   30, 6, 0, 1, itol, 0, 0, 0.0e0, 0.0e0)
      go to 700
 607  call xerrwv(30hlsodes-- iopt (=i1) illegal   ,
     1   30, 7, 0, 1, iopt, 0, 0, 0.0e0, 0.0e0)
      go to 700
 608  call xerrwv(30hlsodes-- mf (=i1) illegal     ,
     1   30, 8, 0, 1, mf, 0, 0, 0.0e0, 0.0e0)
      go to 700
 609  call xerrwv(30hlsodes-- seth (=r1) .lt. 0.0  ,
     1   30, 9, 0, 0, 0, 0, 1, seth, 0.0e0)
      go to 700
 611  call xerrwv(30hlsodes-- maxord (=i1) .lt. 0  ,
     1   30, 11, 0, 1, maxord, 0, 0, 0.0e0, 0.0e0)
      go to 700
 612  call xerrwv(30hlsodes-- mxstep (=i1) .lt. 0  ,
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0e0, 0.0e0)
      go to 700
 613  call xerrwv(30hlsodes-- mxhnil (=i1) .lt. 0  ,
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
      go to 700
 614  call xerrwv(40hlsodes-- tout (=r1) behind t (=r2)      ,
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv(50h      integration direction is given by h0 (=r1)  ,
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0e0)
      go to 700
 615  call xerrwv(30hlsodes-- hmax (=r1) .lt. 0.0  ,
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0e0)
      go to 700
 616  call xerrwv(30hlsodes-- hmin (=r1) .lt. 0.0  ,
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0e0)
      go to 700
 617  call xerrwv(50hlsodes-- rwork length is insufficient to proceed. ,
     1   50, 17, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 618  call xerrwv(50hlsodes-- iwork length is insufficient to proceed. ,
     1   50, 18, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. leniw (=i1), exceeds liw (=i2),
     1   60, 18, 0, 2, leniw, liw, 0, 0.0e0, 0.0e0)
      go to 700
 619  call xerrwv(40hlsodes-- rtol(i1) is r1 .lt. 0.0        ,
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0e0)
      go to 700
 620  call xerrwv(40hlsodes-- atol(i1) is r1 .lt. 0.0        ,
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0e0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv(40hlsodes-- ewt(i1) is r1 .le. 0.0         ,
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0e0)
      go to 700
 622  call xerrwv(
     1  60hlsodes-- tout (=r1) too close to t(=r2) to start integration,
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  60hlsodes-- itask = i1 and tout (=r1) behind tcur - hu (= r2)  ,
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ,
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  60hlsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ,
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv(50hlsodes-- at start of problem, too much accuracy   ,
     1   50, 26, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h      requested for precision of machine..  see tolsf (=r1) ,
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0e0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv(50hlsodes-- trouble from intdy. itask = i1, tout = r1,
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0e0)
      go to 700
 628  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine prep).   ,
     1   60, 28, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 28, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 629  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine jgroup). ,
     1   60, 29, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 29, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 630  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine odrv).   ,
     1   60, 30, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 30, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 631  call xerrwv(
     1  60hlsodes-- error from odrv in yale sparse matrix package      ,
     1   60, 31, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(
     1  60h      at t (=r1), odrv returned error flag = i1*neq + i2.   ,
     1   60, 31, 0, 2, imul, irem, 1, tn, 0.0e0)
      go to 700
 632  call xerrwv(
     1  60hlsodes-- rwork length insufficient (for subroutine cdrv).   ,
     1   60, 32, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(
     1  60h        length needed is .ge. lenrw (=i1), exceeds lrw (=i2),
     1   60, 32, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 633  call xerrwv(
     1  60hlsodes-- error from cdrv in yale sparse matrix package      ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(
     1  60h      at t (=r1), cdrv returned error flag = i1*neq + i2.   ,
     1   60, 33, 0, 2, imul, irem, 1, tn, 0.0e0)
      if (imul .eq. 2) call xerrwv(
     1  60h        duplicate entry in sparsity structure descriptors   ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      if (imul .eq. 3 .or. imul .eq. 6) call xerrwv(
     1  60h        insufficient storage for nsfc (called by cdrv)      ,
     1   60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv(50hlsodes-- repeated occurrences of illegal input    ,
     1   50, 302, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
c
 800  call xerrwv(50hlsodes-- run aborted.. apparent infinite loop     ,
     1   50, 303, 2, 0, 0, 0, 0, 0.0e0, 0.0e0)
      return
c----------------------- end of subroutine lsodes ----------------------
      end
      subroutine adjlr (n, isp, ldif)
      integer n, isp, ldif
      dimension isp(1)
c-----------------------------------------------------------------------
c this routine computes an adjustment, ldif, to the required
c integer storage space in iwk (sparse matrix work space).
c it is called only if the word length ratio is lrat = 1.
c this is to account for the possibility that the symbolic lu phase
c may require more storage than the numerical lu and solution phases.
c-----------------------------------------------------------------------
      integer ip, jlmax, jumax, lnfc, lsfc, nzlu
c
      ip = 2*n + 1
c get jlmax = ijl(n) and jumax = iju(n) (sizes of jl and ju). ----------
      jlmax = isp(ip)
      jumax = isp(ip+ip)
c nzlu = (size of l) + (size of u) = (il(n+1)-il(1)) + (iu(n+1)-iu(1)).
      nzlu = isp(n+1) - isp(1) + isp(ip+n+1) - isp(ip+1)
      lsfc = 12*n + 3 + 2*max0(jlmax,jumax)
      lnfc = 9*n + 2 + jlmax + jumax + nzlu
      ldif = max0(0, lsfc - lnfc)
      return
c----------------------- end of subroutine adjlr -----------------------
      end
      subroutine cdrv
     *     (n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
clll. optimize
c*** subroutine cdrv
c*** driver for subroutines for solving sparse nonsymmetric systems of
c       linear equations (compressed pointer storage)
c
c
c    parameters
c    class abbreviations are--
c       n - integer variable
c       f - real variable
c       v - supplies a value to the driver
c       r - returns a result from the driver
c       i - used internally by the driver
c       a - array
c
c class - parameter
c ------+----------
c       -
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c nv    - n     - number of variables/equations.
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c nva   - ia    - pointers to delimit the rows in a.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c fva   - b     - right-hand side b.  b and z can the same array.
c       -           size = n.
c fra   - z     - solution x.  b and z can be the same array.
c       -           size = n.
c
c         the rows and columns of the original matrix m can be
c    reordered (e.g., to reduce fillin or ensure numerical stability)
c    before calling the driver.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
c    in the original order.
c         if the columns have been reordered (i.e.,  c(i).ne.i  for some
c    i), then the driver will call a subroutine (nroc) which rearranges
c    each row of ja and a, leaving the rows in the original order, but
c    placing the elements of each row in increasing order with respect
c    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
c    been called already.
c
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
c       -           ic(c(i)) = i  for i=1,...,n.
c       -           size = n.
c
c         the solution of the system of linear equations is divided into
c    three stages --
c      nsfc -- the matrix m is processed symbolically to determine where
c               fillin will occur during the numeric factorization.
c      nnfc -- the matrix m is factored numerically into the product ldu
c               of a unit lower triangular matrix l, a diagonal matrix
c               d, and a unit upper triangular matrix u, and the system
c               mx = b  is solved.
c      nnsc -- the linear system  mx = b  is solved using the ldu
c  or           factorization from nnfc.
c      nntc -- the transposed linear system  mt x = b  is solved using
c               the ldu factorization from nnf.
c    for several systems whose coefficient matrices have the same
c    nonzero structure, nsfc need be done only once (for the first
c    system).  then nnfc is done once for each additional system.  for
c    several systems with the same coefficient matrix, nsfc and nnfc
c    need be done only once (for the first system).  then nnsc or nntc
c    is done once for each additional right-hand side.
c
c nv    - path  - path specification.  values and their meanings are --
c       -           1  perform nroc, nsfc, and nnfc.
c       -           2  perform nnfc only  (nsfc is assumed to have been
c       -               done in a manner compatible with the storage
c       -               allocation used in the driver).
c       -           3  perform nnsc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           4  perform nntc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           5  perform nroc and nsfc.
c
c         various errors are detected by the driver and the individual
c    subroutines.
c
c nr    - flag  - error flag.  values and their meanings are --
c       -             0     no errors detected
c       -             n+k   null row in a  --  row = k
c       -            2n+k   duplicate entry in a  --  row = k
c       -            3n+k   insufficient storage in nsfc  --  row = k
c       -            4n+1   insufficient storage in nnfc
c       -            5n+k   null pivot  --  row = k
c       -            6n+k   insufficient storage in nsfc  --  row = k
c       -            7n+1   insufficient storage in nnfc
c       -            8n+k   zero pivot  --  row = k
c       -           10n+1   insufficient storage in cdrv
c       -           11n+1   illegal path specification
c
c         working storage is needed for the factored form of the matrix
c    m plus various temporary vectors.  the arrays isp and rsp should be
c    equivalenced.  integer storage is allocated from the beginning of
c    isp and real storage from the end of rsp.
c
c nv    - nsp   - declared dimension of rsp.  nsp generally must
c       -           be larger than  8n+2 + 2k  (where  k = (number of
c       -           nonzero entries in m)).
c nvira - isp   - integer working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = lratio*nsp.
c fvira - rsp   - real working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = nsp.
c nr    - esp   - if sufficient storage was available to perform the
c       -           symbolic factorization (nsfc), then esp is set to
c       -           the amount of excess storage provided (negative if
c       -           insufficient storage was available to perform the
c       -           numeric factorization (nnfc)).
c
c
c  conversion to double precision
c
c    to convert these routines for double precision arrays..
c    (1) use the double precision declarations in place of the real
c    declarations in each subprogram, as given in comment cards.
c    (2) change the data-loaded value of the integer  lratio
c    in subroutine cdrv, as indicated below.
c    (3) change e0 to d0 in the constants in statement number 10
c    in subroutine nnfc and the line following that.
c
      integer  r(1), c(1), ic(1),  ia(1), ja(1),  isp(1), esp,  path,
     *   flag,  d, u, q, row, tmp, ar,  umax
      real  a(1), b(1), z(1), rsp(1)
c     double precision  a(1), b(1), z(1), rsp(1)
c
c  set lratio equal to the ratio between the length of floating point
c  and integer array data.  e. g., lratio = 1 for (real, integer),
c  lratio = 2 for (double precision, integer)
c
      data lratio/1/
c
      if (path.lt.1 .or. 5.lt.path)  go to 111
c******initialize and divide up temporary storage  *******************
      il   = 1
      ijl  = il  + (n+1)
      iu   = ijl +   n
      iju  = iu  + (n+1)
      irl  = iju +   n
      jrl  = irl +   n
      jl   = jrl +   n
c
c  ******  reorder a if necessary, call nsfc if flag is set  ***********
      if ((path-1) * (path-5) .ne. 0)  go to 5
        max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
        jlmax = max/2
        q     = jl   + jlmax
        ira   = q    + (n+1)
        jra   = ira  +   n
        irac  = jra  +   n
        iru   = irac +   n
        jru   = iru  +   n
        jutmp = jru  +   n
        jumax = lratio*nsp  + 1 - jutmp
        esp = max/lratio
        if (jlmax.le.0 .or. jumax.le.0)  go to 110
c
        do 1 i=1,n
          if (c(i).ne.i)  go to 2
   1      continue
        go to 3
   2    ar = nsp + 1 - n
        call  nroc
     *     (n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
        if (flag.ne.0)  go to 100
c
   3    call  nsfc
     *     (n, r, ic, ia,ja,
     *      jlmax, isp(il), isp(jl), isp(ijl),
     *      jumax, isp(iu), isp(jutmp), isp(iju),
     *      isp(q), isp(ira), isp(jra), isp(irac),
     *      isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
        if(flag .ne. 0)  go to 100
c  ******  move ju next to jl  *****************************************
        jlmax = isp(ijl+n-1)
        ju    = jl + jlmax
        jumax = isp(iju+n-1)
        if (jumax.le.0)  go to 5
        do 4 j=1,jumax
   4      isp(ju+j-1) = isp(jutmp+j-1)
c
c  ******  call remaining subroutines  *********************************
   5  jlmax = isp(ijl+n-1)
      ju    = jl  + jlmax
      jumax = isp(iju+n-1)
      l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
      lmax  = isp(il+n) - 1
      d     = l   + lmax
      u     = d   + n
      row   = nsp + 1 - n
      tmp   = row - n
      umax  = tmp - u
      esp   = umax - (isp(iu+n) - 1)
c
      if ((path-1) * (path-2) .ne. 0)  go to 6
        if (umax.lt.0)  go to 110
        call nnfc
     *     (n,  r, c, ic,  ia, ja, a, z, b,
     *      lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),
     *      umax, isp(iu), isp(ju), isp(iju), rsp(u),
     *      rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
        if(flag .ne. 0)  go to 100
c
   6  if ((path-3) .ne. 0)  go to 7
        call nnsc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
c
   7  if ((path-4) .ne. 0)  go to 8
        call nntc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
   8  return
c
c ** error.. error detected in nroc, nsfc, nnfc, or nnsc
 100  return
c ** error.. insufficient storage
 110  flag = 10*n + 1
      return
c ** error.. illegal path specification
 111  flag = 11*n + 1
      return
      end
      subroutine cfode (meth, elco, tesco)
clll. optimize
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
c-----------------------------------------------------------------------
c cfode is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for all orders and saved.
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfode is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array contains the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the implicit adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array contains test constants used for the
c local error test and the selection of step size and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c nq + 1 if k = 3.
c-----------------------------------------------------------------------
      dimension pc(12)
c
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0e0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0e0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0e0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = float(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0e0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
c----------------------- end of subroutine cfode -----------------------
      end
      subroutine cntnzu (n, ia, ja, nzsut)
      integer n, ia, ja, nzsut
      dimension ia(1), ja(1)
c-----------------------------------------------------------------------
c this routine counts the number of nonzero elements in the strict
c upper triangle of the matrix m + m(transpose), where the sparsity
c structure of m is given by pointer arrays ia and ja.
c this is needed to compute the storage requirements for the
c sparse matrix reordering operation in odrv.
c-----------------------------------------------------------------------
      integer ii, jj, j, jmin, jmax, k, kmin, kmax, num
c
      num = 0
      do 50 ii = 1,n
        jmin = ia(ii)
        jmax = ia(ii+1) - 1
        if (jmin .gt. jmax) go to 50
        do 40 j = jmin,jmax
          if (ja(j) - ii) 10, 40, 30
 10       jj =ja(j)
          kmin = ia(jj)
          kmax = ia(jj+1) - 1
          if (kmin .gt. kmax) go to 30
          do 20 k = kmin,kmax
            if (ja(k) .eq. ii) go to 40
 20         continue
 30       num = num + 1
 40       continue
 50     continue
      nzsut = num
      return
c----------------------- end of subroutine cntnzu ----------------------
      end
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
clll. optimize
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
      integer n, itol
      integer i
      real rtol, atol, ycur, ewt
      dimension rtol(1), atol(1), ycur(n), ewt(n)
c
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine ewset -----------------------
      end
      subroutine intdy (t, k, yh, nyh, dky, iflag)
clll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp
      dimension yh(nyh,1), dky(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c intdy computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0e0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,
     1   30, 51, 0, 1, k, 0, 0, 0.0e0, 0.0e0)
      iflag = -1
      return
 90   call xerrwv(30hintdy--  t (=r1) illegal      ,
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0e0)
      call xerrwv(
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c----------------------- end of subroutine intdy -----------------------
      end
      subroutine iprep (neq, y, rwork, ia, ja, ipflag, f, jac)
clll. optimize
      external f, jac
      integer neq, ia, ja, ipflag
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imax, lewtn, lyhd, lyhn
      real y, rwork
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real rlss
      dimension neq(1), y(1), rwork(1), ia(1), ja(1)
      common /ls0001/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c this routine serves as an interface between the driver and
c subroutine prep.  it is called only if miter is 1 or 2.
c tasks performed here are..
c  * call prep,
c  * reset the required wm segment length lenwk,
c  * move yh back to its final location (following wm in rwork),
c  * reset pointers for yh, savf, ewt, and acor, and
c  * move ewt to its new position if istate = 1.
c ipflag is an output error indication flag.  ipflag = 0 if there was
c no trouble, and ipflag is the value of the prep error flag ipper
c if there was trouble in subroutine prep.
c-----------------------------------------------------------------------
      ipflag = 0
c call prep to do matrix preprocessing operations. ---------------------
      call prep (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt),
     1   rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag, f, jac)
      lenwk = max0(lreq,lwmin)
      if (ipflag .lt. 0) return
c if prep was successful, move yh to end of required space for wm. -----
      lyhn = lwm + lenwk
      if (lyhn .gt. lyh) return
      lyhd = lyh - lyhn
      if (lyhd .eq. 0) go to 20
      imax = lyhn - 1 + lenyhm
      do 10 i = lyhn,imax
 10     rwork(i) = rwork(i+lyhd)
      lyh = lyhn
c reset pointers for savf, ewt, and acor. ------------------------------
 20   lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
      if (istatc .eq. 3) go to 40
c if istate = 1, move ewt (left) to its new position. ------------------
      if (lewtn .gt. lewt) return
      do 30 i = 1,n
 30     rwork(i+lewtn-1) = rwork(i+lewt-1)
 40   lewt = lewtn
      return
c----------------------- end of subroutine iprep -----------------------
      end
      subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)
clll. optimize
      integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier
      dimension ia(1), ja(1), igp(1), jgp(n), incl(n), jdone(n)
c-----------------------------------------------------------------------
c this subroutine constructs groupings of the column indices of
c the jacobian matrix, used in the numerical evaluation of the
c jacobian by finite differences.
c
c input..
c n      = the order of the matrix.
c ia,ja  = sparse structure descriptors of the matrix by rows.
c maxg   = length of available storate in the igp array.
c
c output..
c ngrp   = number of groups.
c jgp    = array of length n containing the column indices by groups.
c igp    = pointer array of length ngrp + 1 to the locations in jgp
c          of the beginning of each group.
c ier    = error indicator.  ier = 0 if no error occurred, or 1 if
c          maxg was insufficient.
c
c incl and jdone are working arrays of length n.
c-----------------------------------------------------------------------
      integer i, j, k, kmin, kmax, ncol, ng
c
      ier = 0
      do 10 j = 1,n
 10     jdone(j) = 0
      ncol = 1
      do 60 ng = 1,maxg
        igp(ng) = ncol
        do 20 i = 1,n
 20       incl(i) = 0
        do 50 j = 1,n
c reject column j if it is already in a group.--------------------------
          if (jdone(j) .eq. 1) go to 50
          kmin = ia(j)
          kmax = ia(j+1) - 1
          do 30 k = kmin,kmax
c reject column j if it overlaps any column already in this group.------
            i = ja(k)
            if (incl(i) .eq. 1) go to 50
 30         continue
c accept column j into group ng.----------------------------------------
          jgp(ncol) = j
          ncol = ncol + 1
          jdone(j) = 1
          do 40 k = kmin,kmax
            i = ja(k)
 40         incl(i) = 1
 50       continue
c stop if this group is empty (grouping is complete).-------------------
        if (ncol .eq. igp(ng)) go to 70
 60     continue
c error return if not all columns were chosen (maxg too small).---------
      if (ncol .le. n) go to 80
      ng = maxg
 70   ngrp = ng - 1
      return
 80   ier = 1
      return
c----------------------- end of subroutine jgroup ----------------------
      end
      subroutine md
     *     (n, ia,ja, max, v,l, head,last,next, mark, flag)
clll. optimize
c***********************************************************************
c  md -- minimum degree algorithm (based on element model)
c***********************************************************************
c
c  description
c
c    md finds a minimum degree ordering of the rows and columns of a
c    general sparse matrix m stored in (ia,ja,a) format.
c    when the structure of m is nonsymmetric, the ordering is that
c    obtained for the symmetric matrix  m + m-transpose.
c
c
c  additional parameters
c
c    max  - declared dimension of the one-dimensional arrays v and l.
c           max must be at least  n+2k,  where k is the number of
c           nonzeroes in the strict upper triangle of m + m-transpose
c
c    v    - integer one-dimensional work array.  dimension = max
c
c    l    - integer one-dimensional work array.  dimension = max
c
c    head - integer one-dimensional work array.  dimension = n
c
c    last - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    next - integer one-dimensional array used to return the inverse of
c           the permutation returned in last.  dimension = n
c
c    mark - integer one-dimensional work array (may be the same as v).
c           dimension = n
c
c    flag - integer error flag.  values and their meanings are -
c             0     no errors detected
c             9n+k  insufficient storage in md
c
c
c  definitions of internal parameters
c
c    ---------+---------------------------------------------------------
c    v(s)     - value field of list entry
c    ---------+---------------------------------------------------------
c    l(s)     - link field of list entry  (0 =) end of list)
c    ---------+---------------------------------------------------------
c    l(vi)    - pointer to element list of uneliminated vertex vi
c    ---------+---------------------------------------------------------
c    l(ej)    - pointer to boundary list of active element ej
c    ---------+---------------------------------------------------------
c    head(d)  - vj =) vj head of d-list d
c             -  0 =) no vertex in d-list d
c
c
c             -                  vi uneliminated vertex
c             -          vi in ek           -       vi not in ek
c    ---------+-----------------------------+---------------------------
c    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
c             -                             -  0 =) vi tail of d-list
c    ---------+-----------------------------+---------------------------
c    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
c             --vk =) compute degree        - vj =) vj last in d-list
c             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
c             -  0 =) do not compute degree -
c    ---------+-----------------------------+---------------------------
c    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
c
c
c             -                   vi eliminated vertex
c             -      ei active element      -           otherwise
c    ---------+-----------------------------+---------------------------
c    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
c             -       to be eliminated      -       to be eliminated
c    ---------+-----------------------------+---------------------------
c    last(vi) -  m =) size of ei = m        - undefined
c    ---------+-----------------------------+---------------------------
c    mark(vi) - -m =) overlap count of ei   - undefined
c             -       with ek = m           -
c             - otherwise nonnegative tag   -
c             -       .lt. mark(vk)         -
c
c-----------------------------------------------------------------------
c
      integer  ia(1), ja(1),  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1),  flag,  tag, dmin, vk,ek, tail
      equivalence  (vk,ek)
c
c----initialization
      tag = 0
      call  mdi
     *   (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return
c
      k = 0
      dmin = 1
c
c----while  k .lt. n  do
   1  if (k.ge.n)  go to 4
c
c------search for vertex of minimum degree
   2    if (head(dmin).gt.0)  go to 3
          dmin = dmin + 1
          go to 2
c
c------remove vertex vk of minimum degree from degree list
   3    vk = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin
c
c------number vertex vk, adjust tag, and tag vk
        k = k+1
        next(vk) = -k
        last(ek) = dmin - 1
        tag = tag + last(ek)
        mark(vk) = tag
c
c------form element ek from uneliminated neighbors of vk
        call  mdm
     *     (vk,tail, v,l, last,next, mark)
c
c------purge inactive elements and do mass elimination
        call  mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c
c------update degrees of uneliminated vertices in ek
        call  mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c
        go to 1
c
c----generate inverse permutation from permutation
   4  do 5 k=1,n
        next(k) = -next(k)
   5    last(next(k)) = k
c
      return
      end
      subroutine mdi
     *     (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
clll. optimize
c***********************************************************************
c  mdi -- initialization
c***********************************************************************
      integer  ia(1), ja(1),  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1), tag,  flag,  sfs, vi,dvi, vj
c
c----initialize degrees, element lists, and degree lists
      do 1 vi=1,n
        mark(vi) = 1
        l(vi) = 0
   1    head(vi) = 0
      sfs = n+1
c
c----create nonzero structure
c----for each nonzero entry a(vi,vj)
      do 6 vi=1,n
        jmin = ia(vi)
        jmax = ia(vi+1) - 1
        if (jmin.gt.jmax)  go to 6
        do 5 j=jmin,jmax
          vj = ja(j)
          if (vj-vi) 2, 5, 4
c
c------if a(vi,vj) is in strict lower triangle
c------check for previous occurrence of a(vj,vi)
   2      lvk = vi
          kmax = mark(vi) - 1
          if (kmax .eq. 0) go to 4
          do 3 k=1,kmax
            lvk = l(lvk)
            if (v(lvk).eq.vj) go to 5
   3        continue
c----for unentered entries a(vi,vj)
   4        if (sfs.ge.max)  go to 101
c
c------enter vj in element list for vi
            mark(vi) = mark(vi) + 1
            v(sfs) = vj
            l(sfs) = l(vi)
            l(vi) = sfs
            sfs = sfs+1
c
c------enter vi in element list for vj
            mark(vj) = mark(vj) + 1
            v(sfs) = vi
            l(sfs) = l(vj)
            l(vj) = sfs
            sfs = sfs+1
   5      continue
   6    continue
c
c----create degree lists and initialize mark vector
      do 7 vi=1,n
        dvi = mark(vi)
        next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        nextvi = next(vi)
        if (nextvi.gt.0)  last(nextvi) = vi
   7    mark(vi) = tag
c
      return
c
c ** error-  insufficient storage
 101  flag = 9*n + vi
      return
      end
      subroutine mdm
     *     (vk,tail, v,l, last,next, mark)
clll. optimize
c***********************************************************************
c  mdm -- form element from uneliminated neighbors of vk
c***********************************************************************
      integer  vk, tail,  v(1), l(1),   last(1), next(1),   mark(1),
     *   tag, s,ls,vs,es, b,lb,vb, blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag and list of uneliminated neighbors
      tag = mark(vk)
      tail = vk
c
c----for each vertex/element vs/es in element list of vk
      ls = l(vk)
   1  s = ls
      if (s.eq.0)  go to 5
        ls = l(s)
        vs = v(s)
        if (next(vs).lt.0)  go to 2
c
c------if vs is uneliminated vertex, then tag and append to list of
c------uneliminated neighbors
          mark(vs) = tag
          l(tail) = s
          tail = s
          go to 4
c
c------if es is active element, then ...
c--------for each vertex vb in boundary list of element es
   2      lb = l(es)
          blpmax = last(es)
          do 3 blp=1,blpmax
            b = lb
            lb = l(b)
            vb = v(b)
c
c----------if vb is untagged vertex, then tag and append to list of
c----------uneliminated neighbors
            if (mark(vb).ge.tag)  go to 3
              mark(vb) = tag
              l(tail) = b
              tail = b
   3        continue
c
c--------mark es inactive
          mark(es) = tag
c
   4    go to 1
c
c----terminate list of uneliminated neighbors
   5  l(tail) = 0
c
      return
      end
      subroutine mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
clll. optimize
c***********************************************************************
c  mdp -- purge inactive elements and do mass elimination
c***********************************************************************
      integer  ek, tail,  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax
c
c----initialize tag
      tag = mark(ek)
c
c----for each vertex vi in ek
      li = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 12
      do 11 ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)
c
c------remove vi from degree list
        if (last(vi).eq.0)  go to 3
          if (last(vi).gt.0)  go to 1
            head(-last(vi)) = next(vi)
            go to 2
   1        next(last(vi)) = next(vi)
   2      if (next(vi).gt.0)  last(next(vi)) = last(vi)
c
c------remove inactive items from element list of vi
   3    ls = vi
   4    s = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).lt.tag)  go to 5
            free = ls
            l(s) = l(ls)
            ls = s
   5      go to 4
c
c------if vi is interior vertex, then remove from list and eliminate
   6    lvi = l(vi)
        if (lvi.ne.0)  go to 7
          l(i) = l(li)
          li = i
c
          k = k+1
          next(vi) = -k
          last(ek) = last(ek) - 1
          go to 11
c
c------else ...
c--------classify vertex vi
   7      if (l(lvi).ne.0)  go to 9
            evi = v(lvi)
            if (next(evi).ge.0)  go to 9
              if (mark(evi).lt.0)  go to 8
c
c----------if vi is prototype vertex, then mark as such, initialize
c----------overlap count for corresponding element, and move vi to end
c----------of boundary list
                last(vi) = evi
                mark(evi) = -1
                l(tail) = li
                tail = li
                l(i) = l(li)
                li = i
                go to 10
c
c----------else if vi is duplicate vertex, then mark as such and adjust
c----------overlap count for corresponding element
   8            last(vi) = 0
                mark(evi) = mark(evi) - 1
                go to 10
c
c----------else mark vi to compute degree
   9            last(vi) = -ek
c
c--------insert ek in element list of vi
  10      v(free) = ek
          l(free) = l(vi)
          l(vi) = free
  11    continue
c
c----terminate boundary list
  12  l(tail) = 0
c
      return
      end
      subroutine mdu
     *     (ek,dmin, v,l, head,last,next, mark)
clll. optimize
c***********************************************************************
c  mdu -- update degrees of uneliminated vertices in ek
c***********************************************************************
      integer  ek, dmin,  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,
     *   blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag
      tag = mark(ek) - last(ek)
c
c----for each vertex vi in ek
      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 11
      do 10 ilp=1,ilpmax
        i = l(i)
        vi = v(i)
        if (last(vi))  1, 10, 8
c
c------if vi neither prototype nor duplicate vertex, then merge elements
c------to compute degree
   1      tag = tag + 1
          dvi = last(ek)
c
c--------for each vertex/element vs/es in element list of vi
          s = l(vi)
   2      s = l(s)
          if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3
c
c----------if vs is uneliminated vertex, then tag and adjust degree
              mark(vs) = tag
              dvi = dvi + 1
              go to 5
c
c----------if es is active element, then expand
c------------check for outmatched vertex
   3          if (mark(es).lt.0)  go to 6
c
c------------for each vertex vb in es
              b = es
              blpmax = last(es)
              do 4 blp=1,blpmax
                b = l(b)
                vb = v(b)
c
c--------------if vb is untagged, then tag and adjust degree
                if (mark(vb).ge.tag)  go to 4
                  mark(vb) = tag
                  dvi = dvi + 1
   4            continue
c
   5        go to 2
c
c------else if vi is outmatched vertex, then adjust overlaps but do not
c------compute degree
   6      last(vi) = 0
          mark(es) = mark(es) - 1
   7      s = l(s)
          if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7
c
c------else if vi is prototype vertex, then calculate degree by
c------inclusion/exclusion and reset overlap count
   8      evi = last(vi)
          dvi = last(ek) + last(evi) + mark(evi)
          mark(evi) = 0
c
c------insert vi in appropriate degree list
   9    next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        if (dvi.lt.dmin)  dmin = dvi
c
  10    continue
c
  11  return
      end
      subroutine nnfc
     *     (n, r,c,ic, ia,ja,a, z, b,
     *      lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,
     *      row, tmp, irl,jrl, flag)
clll. optimize
c*** subroutine nnfc
c*** numerical ldu-factorization of sparse nonsymmetric matrix and
c      solution of system of linear equations (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
c
      integer rk,umax
      integer  r(1), c(1), ic(1), ia(1), ja(1), il(1), jl(1), ijl(1)
      integer  iu(1), ju(1), iju(1), irl(1), jrl(1), flag
      real  a(1), l(1), d(1), u(1), z(1), b(1), row(1)
      real tmp(1), lki, sum, dk
c     double precision  a(1), l(1), d(1), u(1), z(1), b(1), row(1)
c     double precision  tmp(1), lki, sum, dk
c
c  ******  initialize pointers and test storage  ***********************
      if(il(n+1)-1 .gt. lmax) go to 104
      if(iu(n+1)-1 .gt. umax) go to 107
      do 1 k=1,n
        irl(k) = il(k)
        jrl(k) = 0
   1    continue
c
c  ******  for each row  ***********************************************
      do 19 k=1,n
c  ******  reverse jrl and zero row where kth row of l will fill in  ***
        row(k) = 0
        i1 = 0
        if (jrl(k) .eq. 0) go to 3
        i = jrl(k)
   2    i2 = jrl(i)
        jrl(i) = i1
        i1 = i
        row(i) = 0
        i = i2
        if (i .ne. 0) go to 2
c  ******  set row to zero where u will fill in  ***********************
   3    jmin = iju(k)
        jmax = jmin + iu(k+1) - iu(k) - 1
        if (jmin .gt. jmax) go to 5
        do 4 j=jmin,jmax
   4      row(ju(j)) = 0
c  ******  place kth row of a in row  **********************************
   5    rk = r(k)
        jmin = ia(rk)
        jmax = ia(rk+1) - 1
        do 6 j=jmin,jmax
          row(ic(ja(j))) = a(j)
   6      continue
c  ******  initialize sum, and link through jrl  ***********************
        sum = b(rk)
        i = i1
        if (i .eq. 0) go to 10
c  ******  assign the kth row of l and adjust row, sum  ****************
   7      lki = -row(i)
c  ******  if l is not required, then comment out the following line  **
          l(irl(i)) = -lki
          sum = sum + lki * tmp(i)
          jmin = iu(i)
          jmax = iu(i+1) - 1
          if (jmin .gt. jmax) go to 9
          mu = iju(i) - jmin
          do 8 j=jmin,jmax
   8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
   9      i = jrl(i)
          if (i .ne. 0) go to 7
c
c  ******  assign kth row of u and diagonal d, set tmp(k)  *************
  10    if (row(k) .eq. 0.0e0) go to 108
        dk = 1.0e0 / row(k)
        d(k) = dk
        tmp(k) = sum * dk
        if (k .eq. n) go to 19
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 12
        mu = iju(k) - jmin
        do 11 j=jmin,jmax
  11      u(j) = row(ju(mu+j)) * dk
  12    continue
c
c  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1
        if (i .eq. 0) go to 18
  14    irl(i) = irl(i) + 1
        i1 = jrl(i)
        if (irl(i) .ge. il(i+1)) go to 17
        ijlb = irl(i) - il(i) + ijl(i)
        j = jl(ijlb)
  15    if (i .gt. jrl(j)) go to 16
          j = jrl(j)
          go to 15
  16    jrl(i) = jrl(j)
        jrl(j) = i
  17    i = i1
        if (i .ne. 0) go to 14
  18    if (irl(k) .ge. il(k+1)) go to 19
        j = jl(ijl(k))
        jrl(k) = jrl(j)
        jrl(j) = k
  19    continue
c
c  ******  solve  ux = tmp  by back substitution  **********************
      k = n
      do 22 i=1,n
        sum =  tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 21
        mu = iju(k) - jmin
        do 20 j=jmin,jmax
  20      sum = sum - u(j) * tmp(ju(mu+j))
  21    tmp(k) =  sum
        z(c(k)) =  sum
  22    k = k-1
      flag = 0
      return
c
c ** error.. insufficient storage for l
 104  flag = 4*n + 1
      return
c ** error.. insufficient storage for u
 107  flag = 7*n + 1
      return
c ** error.. zero pivot
 108  flag = 8*n + k
      return
      end
      subroutine nnsc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
clll. optimize
c*** subroutine nnsc
c*** numerical solution of sparse nonsymmetric system of linear
c      equations given ldu-factorization (compressed pointer storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving  ly = b.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(1), c(1), il(1), jl(1), ijl(1), iu(1), ju(1), iju(1)
      real l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk, sum
c     double precision  l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(r(k))
c  ******  solve  ly = b  by forward substitution  *********************
      do 3 k=1,n
        jmin = il(k)
        jmax = il(k+1) - 1
        tmpk = -d(k) * tmp(k)
        tmp(k) = -tmpk
        if (jmin .gt. jmax) go to 3
        ml = ijl(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
   3    continue
c  ******  solve  ux = y  by back substitution  ************************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax) go to 5
        mu = iju(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + u(j) * tmp(ju(mu+j))
   5    tmp(k) = -sum
        z(c(k)) = -sum
        k = k - 1
   6    continue
      return
      end
      subroutine nntc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
clll. optimize
c*** subroutine nntc
c*** numeric solution of the transpose of a sparse nonsymmetric system
c      of linear equations given lu-factorization (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving ut y = b
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(1), c(1), il(1), jl(1), ijl(1), iu(1), ju(1), iju(1)
      real l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum
c     double precision l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(c(k))
c  ******  solve  ut y = b  by forward substitution  *******************
      do 3 k=1,n
        jmin = iu(k)
        jmax = iu(k+1) - 1
        tmpk = -tmp(k)
        if (jmin .gt. jmax) go to 3
        mu = iju(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
   3    continue
c  ******  solve  lt x = y  by back substitution  **********************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = il(k)
        jmax = il(k+1) - 1
        if (jmin .gt. jmax) go to 5
        ml = ijl(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + l(j) * tmp(jl(ml+j))
   5    tmp(k) = -sum * d(k)
        z(r(k)) = tmp(k)
        k = k - 1
   6    continue
      return
      end
      subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
clll. optimize
c
c       ----------------------------------------------------------------
c
c               yale sparse matrix package - nonsymmetric codes
c                    solving the system of equations mx = b
c
c    i.   calling sequences
c         the coefficient matrix can be processed by an ordering routine
c    (e.g., to reduce fillin or ensure numerical stability) before using
c    the remaining subroutines.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
c    is used, then nroc should be used to reorder the coefficient matrix
c    the calling sequence is --
c        (       (matrix ordering))
c        (nroc   (matrix reordering))
c         nsfc   (symbolic factorization to determine where fillin will
c                  occur during numeric factorization)
c         nnfc   (numeric factorization into product ldu of unit lower
c                  triangular matrix l, diagonal matrix d, and unit
c                  upper triangular matrix u, and solution of linear
c                  system)
c         nnsc   (solution of linear system for additional right-hand
c                  side using ldu factorization from nnfc)
c    (if only one system of equations is to be solved, then the
c    subroutine trk should be used.)
c
c    ii.  storage of sparse matrices
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c         the strict upper (lower) triangular portion of the matrix
c    u (l) is stored in a similar fashion using the arrays  iu, ju, u
c    (il, jl, l)  except that an additional array iju (ijl) is used to
c    compress storage of ju (jl) by allowing some sequences of column
c    (row) indices to used for more than one row (column)  (n.b., l is
c    stored by columns).  iju(k) (ijl(k)) points to the starting
c    location in ju (jl) of entries for the kth row (column).
c    compression in ju (jl) occurs in two ways.  first, if a row
c    (column) i was merged into the current row (column) k, and the
c    number of elements merged in from (the tail portion of) row
c    (column) i is the same as the final length of row (column) k, then
c    the kth row (column) and the tail of row (column) i are identical
c    and iju(k) (ijl(k)) points to the start of the tail.  second, if
c    some tail portion of the (k-1)st row (column) is identical to the
c    head of the kth row (column), then iju(k) (ijl(k)) points to the
c    start of that tail portion.  for example, the nonzero structure of
c    the strict upper triangular part of the matrix
c            d 0 x x x
c            0 d 0 x x
c            0 0 d x 0
c            0 0 0 d x
c            0 0 0 0 d
c    would be represented as
c                - 1 2 3 4 5 6
c            ----+------------
c             iu - 1 4 6 7 8 8
c             ju - 3 4 5 4
c            iju - 1 2 4 3           .
c    the diagonal entries of l and u are assumed to be equal to one and
c    are not stored.  the array d contains the reciprocals of the
c    diagonal entries of the matrix d.
c
c    iii. additional storage savings
c         in nsfc, r and ic can be the same array in the calling
c    sequence if no reordering of the coefficient matrix has been done.
c         in nnfc, r, c, and ic can all be the same array if no
c    reordering has been done.  if only the rows have been reordered,
c    then c and ic can be the same array.  if the row and column
c    orderings are the same, then r and c can be the same array.  z and
c    row can be the same array.
c         in nnsc or nntc, r and c can be the same array if no
c    reordering has been done or if the row and column orderings are the
c    same.  z and b can be the same array.  however, then b will be
c    destroyed.
c
c    iv.  parameters
c         following is a list of parameters to the programs.  names are
c    uniform among the various subroutines.  class abbreviations are --
c       n - integer variable
c       f - real variable
c       v - supplies a value to a subroutine
c       r - returns a result from a subroutine
c       i - used internally by a subroutine
c       a - array
c
c class - parameter
c ------+----------
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c fva   - b     - right-hand side b.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c fvra  - d     - reciprocals of the diagonal entries of the matrix d.
c       -           size = n.
c nr    - flag  - error flag.  values and their meanings are --
c       -            0     no errors detected
c       -            n+k   null row in a  --  row = k
c       -           2n+k   duplicate entry in a  --  row = k
c       -           3n+k   insufficient storage for jl  --  row = k
c       -           4n+1   insufficient storage for l
c       -           5n+k   null pivot  --  row = k
c       -           6n+k   insufficient storage for ju  --  row = k
c       -           7n+1   insufficient storage for u
c       -           8n+k   zero pivot  --  row = k
c nva   - ia    - pointers to delimit the rows of a.
c       -           size = n+1.
c nvra  - ijl   - pointers to the first element in each column in jl,
c       -           used to compress storage in jl.
c       -           size = n.
c nvra  - iju   - pointers to the first element in each row in ju, used
c       -           to compress storage in ju.
c       -           size = n.
c nvra  - il    - pointers to delimit the columns of l.
c       -           size = n+1.
c nvra  - iu    - pointers to delimit the rows of u.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c nvra  - jl    - row numbers corresponding to the elements of l.
c       -           size = jlmax.
c nv    - jlmax - declared dimension of jl.  jlmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin minus compression.
c nvra  - ju    - column numbers corresponding to the elements of u.
c       -           size = jumax.
c nv    - jumax - declared dimension of ju.  jumax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin minus compression.
c fvra  - l     - nonzero entries in the strict lower triangular portion
c       -           of the matrix l, stored by columns.
c       -           size = lmax.
c nv    - lmax  - declared dimension of l.  lmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin  (il(n+1)-1 after nsfc).
c nv    - n     - number of variables/equations.
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c fvra  - u     - nonzero entries in the strict upper triangular portion
c       -           of the matrix u, stored by rows.
c       -           size = umax.
c nv    - umax  - declared dimension of u.  umax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin  (iu(n+1)-1 after nsfc).
c fra   - z     - solution x.
c       -           size = n.
c
c       ----------------------------------------------------------------
c
c*** subroutine nroc
c*** reorders rows of a, leaving row order unchanged
c
c
c       input parameters.. n, ic, ia, ja, a
c       output parameters.. ja, a, flag
c
c       parameters used internally..
c nia   - p     - at the kth step, p is a linked list of the reordered
c       -           column indices of the kth row of a.  p(n+1) points
c       -           to the first entry in the list.
c       -           size = n+1.
c nia   - jar   - at the kth step,jar contains the elements of the
c       -           reordered column indices of a.
c       -           size = n.
c fia   - ar    - at the kth step, ar contains the elements of the
c       -           reordered row of a.
c       -           size = n.
c
      integer  ic(1), ia(1), ja(1), jar(1), p(1), flag
      real  a(1), ar(1)
c     double precision  a(1), ar(1)
c
c  ******  for each nonempty row  *******************************
      do 5 k=1,n
        jmin = ia(k)
        jmax = ia(k+1) - 1
        if(jmin .gt. jmax) go to 5
        p(n+1) = n + 1
c  ******  insert each element in the list  *********************
        do 3 j=jmin,jmax
          newj = ic(ja(j))
          i = n + 1
   1      if(p(i) .ge. newj) go to 2
            i = p(i)
            go to 1
   2      if(p(i) .eq. newj) go to 102
          p(newj) = p(i)
          p(i) = newj
          jar(newj) = ja(j)
          ar(newj) = a(j)
   3      continue
c  ******  replace old row in ja and a  *************************
        i = n + 1
        do 4 j=jmin,jmax
          i = p(i)
          ja(j) = jar(i)
   4      a(j) = ar(i)
   5    continue
      flag = 0
      return
c
c ** error.. duplicate entry in a
 102  flag = n + k
      return
      end
      subroutine nsfc
     *      (n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,
     *       q, ira,jra, irac, irl,jrl, iru,jru, flag)
clll. optimize
c*** subroutine nsfc
c*** symbolic ldu-factorization of nonsymmetric sparse matrix
c      (compressed pointer storage)
c
c
c       input variables.. n, r, ic, ia, ja, jlmax, jumax.
c       output variables.. il, jl, ijl, iu, ju, iju, flag.
c
c       parameters used internally..
c nia   - q     - suppose  m*  is the result of reordering  m.  if
c       -           processing of the ith row of  m*  (hence the ith
c       -           row of  u) is being done,  q(j)  is initially
c       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
c       -           values need not be stored, each entry points to the
c       -           next nonzero and  q(n+1)  points to the first.  n+1
c       -           indicates the end of the list.  for example, if n=9
c       -           and the 5th row of  m*  is
c       -              0 x x 0 x 0 0 x 0
c       -           then  q  will initially be
c       -              a a a a 8 a a 10 5           (a - arbitrary).
c       -           as the algorithm proceeds, other elements of  q
c       -           are inserted in the list because of fillin.
c       -           q  is used in an analogous manner to compute the
c       -           ith column of  l.
c       -           size = n+1.
c nia   - ira,  - vectors used to find the columns of  m.  at the kth
c nia   - jra,      step of the factorization,  irac(k)  points to the
c nia   - irac      head of a linked list in  jra  of row indices i
c       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
c       -           indicates the end of the list.  ira(i)  (i.ge.k)
c       -           points to the smallest j such that j .ge. k and
c       -           m(i,j)  is nonzero.
c       -           size of each = n.
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c nia   - iru,  - vectors used in a manner analogous to  irl and jrl
c nia   - jru       to find the columns of  u.
c       -           size of each = n.
c
c  internal variables..
c    jlptr - points to the last position used in  jl.
c    juptr - points to the last position used in  ju.
c    jmin,jmax - are the indices in  a or u  of the first and last
c                elements to be examined in a given row.
c                for example,  jmin=ia(k), jmax=ia(k+1)-1.
c
      integer cend, qm, rend, rk, vj
      integer ia(1), ja(1), ira(1), jra(1), il(1), jl(1), ijl(1)
      integer iu(1), ju(1), iju(1), irl(1), jrl(1), iru(1), jru(1)
      integer r(1), ic(1), q(1), irac(1), flag
c
c  ******  initialize pointers  ****************************************
      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do 1 k=1,n
        irac(k) = 0
        jra(k) = 0
        jrl(k) = 0
   1    jru(k) = 0
c  ******  initialize column pointers for a  ***************************
      do 2 k=1,n
        rk = r(k)
        iak = ia(rk)
        if (iak .ge. ia(rk+1))  go to 101
        jaiak = ic(ja(iak))
        if (jaiak .gt. k)  go to 105
        jra(k) = irac(jaiak)
        irac(jaiak) = k
   2    ira(k) = iak
c
c  ******  for each column of l and row of u  **************************
      do 41 k=1,n
c
c  ******  initialize q for computing kth column of l  *****************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth column of a  ******************************
        vj = irac(k)
        if (vj .eq. 0)  go to 5
   3      qm = np1
   4      m = qm
          qm =  q(m)
          if (qm .lt. vj)  go to 4
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if (vj .ne. 0)  go to 3
c  ******  link through jru  *******************************************
   5    lastid = 0
        lasti = 0
        ijl(k) = jlptr
        i = k
   6      i = jru(i)
          if (i .eq. 0)  go to 10
          qm = np1
          jmin = irl(i)
          jmax = ijl(i) + il(i+1) - il(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 6
          jtmp = jl(jmin)
          if (jtmp .ne. k)  long = long + 1
          if (jtmp .eq. k)  r(i) = -r(i)
          if (lastid .ge. long)  go to 7
            lasti = i
            lastid = long
c  ******  and merge the corresponding columns into the kth column  ****
   7      do 9 j=jmin,jmax
            vj = jl(j)
   8        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 8
            if (qm .eq. vj)  go to 9
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
   9        continue
            go to 6
c  ******  lasti is the longest column merged into the kth  ************
c  ******  see if it equals the entire kth column  *********************
  10    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 17
        if (lastid .ne. luk)  go to 11
c  ******  if so, jl can be compressed  ********************************
        irll = irl(lasti)
        ijl(k) = irll + 1
        if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
        go to 17
c  ******  if not, see if kth column can overlap the previous one  *****
  11    if (jlmin .gt. jlptr)  go to 15
        qm = q(qm)
        do 12 j=jlmin,jlptr
          if (jl(j) - qm)  12, 13, 15
  12      continue
        go to 15
  13    ijl(k) = j
        do 14 i=j,jlptr
          if (jl(i) .ne. qm)  go to 15
          qm = q(qm)
          if (qm .gt. n)  go to 17
  14      continue
        jlptr = j - 1
c  ******  move column indices from q to jl, update vectors  ***********
  15    jlmin = jlptr + 1
        ijl(k) = jlmin
        if (luk .eq. 0)  go to 17
        jlptr = jlptr + luk
        if (jlptr .gt. jlmax)  go to 103
          qm = q(np1)
          do 16 j=jlmin,jlptr
            qm = q(qm)
  16        jl(j) = qm
  17    irl(k) = ijl(k)
        il(k+1) = il(k) + luk
c
c  ******  initialize q for computing kth row of u  ********************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth row of reordered a  ***********************
        rk = r(k)
        jmin = ira(k)
        jmax = ia(rk+1) - 1
        if (jmin .gt. jmax)  go to 20
        do 19 j=jmin,jmax
          vj = ic(ja(j))
          qm = np1
  18      m = qm
          qm = q(m)
          if (qm .lt. vj)  go to 18
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
  19      continue
c  ******  link through jrl,  ******************************************
  20    lastid = 0
        lasti = 0
        iju(k) = juptr
        i = k
        i1 = jrl(k)
  21      i = i1
          if (i .eq. 0)  go to 26
          i1 = jrl(i)
          qm = np1
          jmin = iru(i)
          jmax = iju(i) + iu(i+1) - iu(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 21
          jtmp = ju(jmin)
          if (jtmp .eq. k)  go to 22
c  ******  update irl and jrl, *****************************************
            long = long + 1
            cend = ijl(i) + il(i+1) - il(i)
            irl(i) = irl(i) + 1
            if (irl(i) .ge. cend)  go to 22
              j = jl(irl(i))
              jrl(i) = jrl(j)
              jrl(j) = i
  22      if (lastid .ge. long)  go to 23
            lasti = i
            lastid = long
c  ******  and merge the corresponding rows into the kth row  **********
  23      do 25 j=jmin,jmax
            vj = ju(j)
  24        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 24
            if (qm .eq. vj)  go to 25
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
  25        continue
          go to 21
c  ******  update jrl(k) and irl(k)  ***********************************
  26    if (il(k+1) .le. il(k))  go to 27
          j = jl(irl(k))
          jrl(k) = jrl(j)
          jrl(j) = k
c  ******  lasti is the longest row merged into the kth  ***************
c  ******  see if it equals the entire kth row  ************************
  27    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 34
        if (lastid .ne. luk)  go to 28
c  ******  if so, ju can be compressed  ********************************
        irul = iru(lasti)
        iju(k) = irul + 1
        if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
        go to 34
c  ******  if not, see if kth row can overlap the previous one  ********
  28    if (jumin .gt. juptr)  go to 32
        qm = q(qm)
        do 29 j=jumin,juptr
          if (ju(j) - qm)  29, 30, 32
  29      continue
        go to 32
  30    iju(k) = j
        do 31 i=j,juptr
          if (ju(i) .ne. qm)  go to 32
          qm = q(qm)
          if (qm .gt. n)  go to 34
  31      continue
        juptr = j - 1
c  ******  move row indices from q to ju, update vectors  **************
  32    jumin = juptr + 1
        iju(k) = jumin
        if (luk .eq. 0)  go to 34
        juptr = juptr + luk
        if (juptr .gt. jumax)  go to 106
          qm = q(np1)
          do 33 j=jumin,juptr
            qm = q(qm)
  33        ju(j) = qm
  34    iru(k) = iju(k)
        iu(k+1) = iu(k) + luk
c
c  ******  update iru, jru  ********************************************
        i = k
  35      i1 = jru(i)
          if (r(i) .lt. 0)  go to 36
          rend = iju(i) + iu(i+1) - iu(i)
          if (iru(i) .ge. rend)  go to 37
            j = ju(iru(i))
            jru(i) = jru(j)
            jru(j) = i
            go to 37
  36      r(i) = -r(i)
  37      i = i1
          if (i .eq. 0)  go to 38
          iru(i) = iru(i) + 1
          go to 35
c
c  ******  update ira, jra, irac  **************************************
  38    i = irac(k)
        if (i .eq. 0)  go to 41
  39      i1 = jra(i)
          ira(i) = ira(i) + 1
          if (ira(i) .ge. ia(r(i)+1))  go to 40
          irai = ira(i)
          jairai = ic(ja(irai))
          if (jairai .gt. i)  go to 40
          jra(i) = irac(jairai)
          irac(jairai) = i
  40      i = i1
          if (i .ne. 0)  go to 39
  41    continue
c
      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return
c
c ** error.. null row in a
 101  flag = n + rk
      return
c ** error.. duplicate entry in a
 102  flag = 2*n + rk
      return
c ** error.. insufficient storage for jl
 103  flag = 3*n + k
      return
c ** error.. null pivot
 105  flag = 5*n + k
      return
c ** error.. insufficient storage for ju
 106  flag = 6*n + k
      return
      end
      subroutine odrv
     *     (n, ia,ja,a, p,ip, nsp,isp, path, flag)
clll. optimize
c                                                                 5/2/83
c***********************************************************************
c  odrv -- driver for sparse matrix reordering routines
c***********************************************************************
c
c  description
c
c    odrv finds a minimum degree ordering of the rows and columns
c    of a matrix m stored in (ia,ja,a) format (see below).  for the
c    reordered matrix, the work and storage required to perform
c    gaussian elimination is (usually) significantly less.
c
c    note.. odrv and its subordinate routines have been modified to
c    compute orderings for general matrices, not necessarily having any
c    symmetry.  the miminum degree ordering is computed for the
c    structure of the symmetric matrix  m + m-transpose.
c    modifications to the original odrv module have been made in
c    the coding in subroutine mdi, and in the initial comments in
c    subroutines odrv and md.
c
c    if only the nonzero entries in the upper triangle of m are being
c    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c    with the diagonal entries placed first in each row.  this is to
c    ensure that if m(i,j) will be in the upper triangle of m with
c    respect to the new ordering, then m(i,j) is stored in row i (and
c    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
c    strict lower triangle of m, then m(j,i) is stored in row j (and
c    thus m(i,j) is not stored).
c
c
c  storage of sparse matrices
c
c    the nonzero entries of the matrix m are stored row-by-row in the
c    array a.  to identify the individual nonzero entries in each row,
c    we need to know in which column each entry lies.  these column
c    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  to identify the individual rows, we need to know where
c    each row starts.  these row pointers are stored in the array ia.
c    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
c    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
c    the first location following the last element in the last row.
c    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
c    the nonzero entries in the i-th row are stored consecutively in
c
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c
c    and the corresponding column indices are stored consecutively in
c
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c
c    when the coefficient matrix is symmetric, only the nonzero entries
c    in the upper triangle need be stored.  for example, the matrix
c
c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )
c
c    could be stored as
c
c            - 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia - 1  4  5  8 12 14
c         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
c
c    or (symmetrically) as
c
c            - 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia - 1  4  5  7  9 10
c         ja - 1  3  4  2  3  4  4  5  5
c          a - 1  2  3  4  5  6  7  8  9          .
c
c
c  parameters
c
c    n    - order of the matrix
c
c    ia   - integer one-dimensional array containing pointers to delimit
c           rows in ja and a.  dimension = n+1
c
c    ja   - integer one-dimensional array containing the column indices
c           corresponding to the elements of a.  dimension = number of
c           nonzero entries in (the upper triangle of) m
c
c    a    - real one-dimensional array containing the nonzero entries in
c           (the upper triangle of) m, stored by rows.  dimension =
c           number of nonzero entries in (the upper triangle of) m
c
c    p    - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    ip   - integer one-dimensional array used to return the inverse of
c           the permutation returned in p.  dimension = n
c
c    nsp  - declared dimension of the one-dimensional array isp.  nsp
c           must be at least  3n+4k,  where k is the number of nonzeroes
c           in the strict upper triangle of m
c
c    isp  - integer one-dimensional array used for working storage.
c           dimension = nsp
c
c    path - integer path specification.  values and their meanings are -
c             1  find minimum degree ordering only
c             2  find minimum degree ordering and reorder symmetrically
c                  stored matrix (used when only the nonzero entries in
c                  the upper triangle of m are being stored)
c             3  reorder symmetrically stored matrix as specified by
c                  input permutation (used when an ordering has already
c                  been determined and only the nonzero entries in the
c                  upper triangle of m are being stored)
c             4  same as 2 but put diagonal entries at start of each row
c             5  same as 3 but put diagonal entries at start of each row
c
c    flag - integer error flag.  values and their meanings are -
c               0    no errors detected
c              9n+k  insufficient storage in md
c             10n+1  insufficient storage in odrv
c             11n+1  illegal path specification
c
c
c  conversion from real to double precision
c
c    change the real declarations in odrv and sro to double precision
c    declarations.
c
c-----------------------------------------------------------------------
c
      integer  ia(1), ja(1),  p(1), ip(1),  isp(1),  path,  flag,
     *   v, l, head,  tmp, q
      real  a(1)
c...  double precision  a(1)
      logical  dflag
c
c----initialize error flag and validate path specification
      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111
c
c----allocate storage and find minimum degree ordering
      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110
c
        call  md
     *     (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100
c
c----allocate storage and symmetrically reorder matrix
   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110
c
        dflag = path.eq.4 .or. path.eq.5
        call sro
     *     (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
c
   2  return
c
c ** error -- error detected in md
 100  return
c ** error -- insufficient storage
 110  flag = 10*n + 1
      return
c ** error -- illegal path specified
 111  flag = 11*n + 1
      return
      end
      subroutine prep (neq, y, yh, savf, ewt, ftem, ia, ja,
     1                     wk, iwk, ipper, f, jac)
clll. optimize
      external f,jac
      integer neq, ia, ja, iwk, ipper
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k,
     1   knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
      real y, yh, savf, ewt, ftem, wk
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth
      real dq, dyj, erwt, fac, yj
      dimension neq(1), y(1), yh(1), savf(1), ewt(1), ftem(1),
     1   ia(1), ja(1), wk(1), iwk(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c this routine performs preprocessing related to the sparse linear
c systems that must be solved if miter = 1 or 2.
c the operations that are performed here are..
c  * compute sparseness structure of jacobian according to moss,
c  * compute grouping of column indices (miter = 2),
c  * compute a new ordering of rows and columns of the matrix,
c  * reorder ja corresponding to the new ordering,
c  * perform a symbolic lu factorization of the matrix, and
c  * set pointers for segments of the iwk/wk array.
c in addition to variables described previously, prep uses the
c following for communication..
c yh     = the history array.  only the first column, containing the
c          current y vector, is used.  used only if moss .ne. 0.
c savf   = a work array of length neq, used only if moss .ne. 0.
c ewt    = array of length neq containing (inverted) error weights.
c          used only if moss = 2 or if istate = moss = 1.
c ftem   = a work array of length neq, identical to acor in the driver,
c          used only if moss = 2.
c wk     = a real work array of length lenwk, identical to wm in
c          the driver.
c iwk    = integer work array, assumed to occupy the same space as wk.
c lenwk  = the length of the work arrays wk and iwk.
c istatc = a copy of the driver input argument istate (= 1 on the
c          first call, = 3 on a continuation call).
c iys    = flag value from odrv or cdrv.
c ipper  = output error flag with the following values and meanings..
c          0  no error.
c         -1  insufficient storage for internal structure pointers.
c         -2  insufficient storage for jgroup.
c         -3  insufficient storage for odrv.
c         -4  other error flag from odrv (should never occur).
c         -5  insufficient storage for cdrv.
c         -6  other error flag from cdrv.
c-----------------------------------------------------------------------
      ibian = lrat*2
      ipian = ibian + 1
      np1 = n + 1
      ipjan = ipian + np1
      ibjan = ipjan - 1
      liwk = lenwk*lrat
      if (ipjan+n-1 .gt. liwk) go to 210
      if (moss .eq. 0) go to 30
c
      if (istatc .eq. 3) go to 20
c istate = 1 and moss .ne. 0.  perturb y for structure determination. --
      do 10 i = 1,n
        erwt = 1.0e0/ewt(i)
        fac = 1.0e0 + 1.0e0/(float(i)+1.0e0)
        y(i) = y(i) + fac*sign(erwt,y(i))
 10     continue
      go to (70, 100), moss
c
 20   continue
c istate = 3 and moss .ne. 0.  load y from yh(*,1). --------------------
      do 25 i = 1,n
 25     y(i) = yh(i)
      go to (70, 100), moss
c
c moss = 0.  process user-s ia,ja.  add diagonal entries if necessary. -
 30   knew = ipjan
      kmin = ia(1)
      iwk(ipian) = 1
      do 60 j = 1,n
        jfound = 0
        kmax = ia(j+1) - 1
        if (kmin .gt. kmax) go to 45
        do 40 k = kmin,kmax
          i = ja(k)
          if (i .eq. j) jfound = 1
          if (knew .gt. liwk) go to 210
          iwk(knew) = i
          knew = knew + 1
 40       continue
        if (jfound .eq. 1) go to 50
 45     if (knew .gt. liwk) go to 210
        iwk(knew) = j
        knew = knew + 1
 50     iwk(ipian+j) = knew + 1 - ipjan
        kmin = kmax + 1
 60     continue
      go to 140
c
c moss = 1.  compute structure from user-supplied jacobian routine jac.
 70   continue
c a dummy call to f allows user to create temporaries for use in jac. --
      call f (neq, tn, y, savf)
      k = ipjan
      iwk(ipian) = 1
      do 90 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        do 75 i = 1,n
 75       savf(i) = 0.0e0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), savf)
        do 80 i = 1,n
          if (abs(savf(i)) .le. seth) go to 80
          if (i .eq. j) go to 80
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 80       continue
        iwk(ipian+j) = k + 1 - ipjan
 90     continue
      go to 140
c
c moss = 2.  compute structure from results of n + 1 calls to f. -------
 100  k = ipjan
      iwk(ipian) = 1
      call f (neq, tn, y, savf)
      do 120 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        yj = y(j)
        erwt = 1.0e0/ewt(j)
        dyj = sign(erwt,yj)
        y(j) = yj + dyj
        call f (neq, tn, y, ftem)
        y(j) = yj
        do 110 i = 1,n
          dq = (ftem(i) - savf(i))/dyj
          if (abs(dq) .le. seth) go to 110
          if (i .eq. j) go to 110
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 110      continue
        iwk(ipian+j) = k + 1 - ipjan
 120    continue
c
 140  continue
      if (moss .eq. 0 .or. istatc .ne. 1) go to 150
c if istate = 1 and moss .ne. 0, restore y from yh. --------------------
      do 145 i = 1,n
 145    y(i) = yh(i)
 150  nnz = iwk(ipian+n) - 1
      lenigp = 0
      ipigp = ipjan + nnz
      if (miter .ne. 2) go to 160
c
c compute grouping of column indices (miter = 2). ----------------------
      maxg = np1
      ipjgp = ipjan + nnz
      ibjgp = ipjgp - 1
      ipigp = ipjgp + n
      iptt1 = ipigp + np1
      iptt2 = iptt1 + n
      lreq = iptt2 + n - 1
      if (lreq .gt. liwk) go to 220
      call jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp),
     1   iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
      if (ier .ne. 0) go to 220
      lenigp = ngp + 1
c
c compute new ordering of rows/columns of jacobian. --------------------
 160  ipr = ipigp + lenigp
      ipc = ipr
      ipic = ipc + n
      ipisp = ipic + n
      iprsp = (ipisp - 2)/lrat + 2
      iesp = lenwk + 1 - iprsp
      if (iesp .lt. 0) go to 230
      ibr = ipr - 1
      do 170 i = 1,n
 170    iwk(ibr+i) = i
      nsp = liwk + 1 - ipisp
      call odrv (n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic),
     1   nsp, iwk(ipisp), 1, iys)
      if (iys .eq. 11*n+1) go to 240
      if (iys .ne. 0) go to 230
c
c reorder jan and do symbolic lu factorization of matrix. --------------
      ipa = lenwk + 1 - nnz
      nsp = ipa - iprsp
      lreq = max0(12*n/lrat, 6*n/lrat+2*n+nnz) + 3
      lreq = lreq + iprsp - 1 + nnz
      if (lreq .gt. lenwk) go to 250
      iba = ipa - 1
      do 180 i = 1,nnz
 180    wk(iba+i) = 0.0e0
      ipisp = lrat*(iprsp - 1) + 1
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
      lreq = lenwk - iesp
      if (iys .eq. 10*n+1) go to 250
      if (iys .ne. 0) go to 260
      ipil = ipisp
      ipiu = ipil + 2*n + 1
      nzu = iwk(ipil+n) - iwk(ipil)
      nzl = iwk(ipiu+n) - iwk(ipiu)
      if (lrat .gt. 1) go to 190
      call adjlr (n, iwk(ipisp), ldif)
      lreq = lreq + ldif
 190  continue
      if (lrat .eq. 2 .and. nnz .eq. n) lreq = lreq + 1
      nsp = nsp + lreq - lenwk
      ipa = lreq + 1 - nnz
      iba = ipa - 1
      ipper = 0
      return
c
 210  ipper = -1
      lreq = 2 + (2*n + 1)/lrat
      lreq = max0(lenwk+1,lreq)
      return
c
 220  ipper = -2
      lreq = (lreq - 1)/lrat + 1
      return
c
 230  ipper = -3
      call cntnzu (n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
      return
c
 240  ipper = -4
      return
c
 250  ipper = -5
      return
c
 260  ipper = -6
      lreq = lenwk
      return
c----------------------- end of subroutine prep ------------------------
      end
      subroutine prjs (neq,y,yh,nyh,ewt,ftem,savf,wk,iwk,f,jac)
clll. optimize
      external f,jac
      integer neq, nyh, iwk
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
      real y, yh, ewt, ftem, savf, wk
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth
      real con, di, fac, hl0, pij, r, r0, rcon, rcont,
     1   srur, vnorm
      dimension neq(1), y(1), yh(nyh,1), ewt(1), ftem(1), savf(1),
     1   wk(1), iwk(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c prjs is called to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c j is computed by columns, either by the user-supplied routine jac
c if miter = 1, or by finite differencing if miter = 2.
c if miter = 3, a diagonal approximation to j is used.
c if miter = 1 or 2, and if the existing value of the jacobian
c (as contained in p) is considered acceptable, then a new value of
c p is reconstructed from the old value.  in any case, when miter
c is 1 or 2, the p matrix is subjected to lu decomposition in cdrv.
c p and its lu decomposition are stored (separately) in wk.
c
c in addition to variables described previously, communication
c with prjs uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stode).
c savf  = array containing f evaluated at predicted y.
c wk    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3, and p and its sparse
c         lu decomposition if miter is 1 or 2.
c         storage of matrix elements starts at wk(3).
c         wk also contains the following matrix-related data..
c         wk(1) = sqrt(uround), used in numerical jacobian increments.
c         wk(2) = h*el0, saved for later use if miter = 3.
c iwk   = integer work space for matrix-related data, assumed to
c         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
c         are assumed to have identical locations.
c el0   = el(1) (input).
c ierpj = output error flag (in common).
c       = 0 if no error.
c       = 1  if zero pivot found in cdrv.
c       = 2  if a singular matrix arose with miter = 3.
c       = -1 if insufficient storage for cdrv (should not occur here).
c       = -2 if other error found in cdrv (should not occur here).
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses other variables in common.
c-----------------------------------------------------------------------
      hl0 = h*el0
      con = -hl0
      if (miter .eq. 3) go to 300
c see whether j should be reevaluated (jok = 0) or not (jok = 1). ------
      jok = 1
      if (nst .eq. 0 .or. nst .ge. nslj+msbj) jok = 0
      if (icf .eq. 1 .and. abs(rc - 1.0e0) .lt. ccmxj) jok = 0
      if (icf .eq. 2) jok = 0
      if (jok .eq. 1) go to 250
c
c miter = 1 or 2, and the jacobian is to be reevaluated. ---------------
 20   jcur = 1
      nje = nje + 1
      nslj = nst
      iplost = 0
      conmin = abs(con)
      go to (100, 200), miter
c
c if miter = 1, call jac, multiply by scalar, and add identity. --------
 100  continue
      kmin = iwk(ipian)
      do 130 j = 1, n
        kmax = iwk(ipian+j) - 1
        do 110 i = 1,n
 110      ftem(i) = 0.0e0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem)
        do 120 k = kmin, kmax
          i = iwk(ibjan+k)
          wk(iba+k) = ftem(i)*con
          if (i .eq. j) wk(iba+k) = wk(iba+k) + 1.0e0
 120      continue
        kmin = kmax + 1
 130    continue
      go to 290
c
c if miter = 2, make ngp calls to f to approximate j and p. ------------
 200  continue
      fac = vnorm(n, savf, ewt)
      r0 = 1000.0e0 * abs(h) * uround * float(n) * fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      srur = wk(1)
      jmin = iwk(ipigp)
      do 240 ng = 1,ngp
        jmax = iwk(ipigp+ng) - 1
        do 210 j = jmin,jmax
          jj = iwk(ibjgp+j)
          r = amax1(srur*abs(y(jj)),r0/ewt(jj))
 210      y(jj) = y(jj) + r
        call f (neq, tn, y, ftem)
        do 230 j = jmin,jmax
          jj = iwk(ibjgp+j)
          y(jj) = yh(jj,1)
          r = amax1(srur*abs(y(jj)),r0/ewt(jj))
          fac = -hl0/r
          kmin =iwk(ibian+jj)
          kmax =iwk(ibian+jj+1) - 1
          do 220 k = kmin,kmax
            i = iwk(ibjan+k)
            wk(iba+k) = (ftem(i) - savf(i))*fac
            if (i .eq. jj) wk(iba+k) = wk(iba+k) + 1.0e0
 220        continue
 230      continue
        jmin = jmax + 1
 240    continue
      nfe = nfe + ngp
      go to 290
c
c if jok = 1, reconstruct new p from old p. ----------------------------
 250  jcur = 0
      rcon = con/con0
      rcont = abs(con)/conmin
      if (rcont .gt. rbig .and. iplost .eq. 1) go to 20
      kmin = iwk(ipian)
      do 275 j = 1,n
        kmax = iwk(ipian+j) - 1
        do 270 k = kmin,kmax
          i = iwk(ibjan+k)
          pij = wk(iba+k)
          if (i .ne. j) go to 260
          pij = pij - 1.0e0
          if (abs(pij) .ge. psmall) go to 260
            iplost = 1
            conmin = amin1(abs(con0),conmin)
 260      pij = pij*rcon
          if (i .eq. j) pij = pij + 1.0e0
          wk(iba+k) = pij
 270      continue
        kmin = kmax + 1
 275    continue
c
c do numerical factorization of p matrix. ------------------------------
 290  nlu = nlu + 1
      con0 = con
      ierpj = 0
      do 295 i = 1,n
 295    ftem(i) = 0.0e0
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
      if (iys .eq. 0) return
      imul = (iys - 1)/n
      ierpj = -2
      if (imul .eq. 8) ierpj = 1
      if (imul .eq. 10) ierpj = -1
      return
c
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  continue
      jcur = 1
      nje = nje + 1
      wk(2) = hl0
      ierpj = 0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wk(3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1e0*r0 - h*(wk(i+2) - savf(i))
        wk(i+2) = 1.0e0
        if (abs(r0) .lt. uround/ewt(i)) go to 320
        if (abs(di) .eq. 0.0e0) go to 330
        wk(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ierpj = 2
      return
c----------------------- end of subroutine prjs ------------------------
      end
      subroutine slss (wk, iwk, x, tem)
clll. optimize
      integer iwk
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     3   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i
      real wk, x, tem
      real rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real rlss
      real di, hl0, phl0, r
      dimension wk(1), iwk(1), x(1), tem(1)
      common /ls0001/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),
     1   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     4   nslj, ngp, nlu, nnz, nsp, nzl, nzu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls cdrv to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c communication with slss uses the following variables..
c wk    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wk(3).
c         wk also contains the following matrix-related data..
c         wk(1) = sqrt(uround) (not used here),
c         wk(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwk   = integer work space for matrix-related data, assumed to
c         be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp)
c         are assumed to have identical locations.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).
c         iersl = 0  if no trouble occurred.
c         iersl = -1 if cdrv returned an error flag (miter = 1 or 2).
c                    this should never occur and is considered fatal.
c         iersl = 1  if a singular matrix arose with miter = 3.
c this routine also uses other variables in common.
c-----------------------------------------------------------------------
      iersl = 0
      go to (100, 100, 300), miter
 100  call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),
     1   wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
      if (iersl .ne. 0) iersl = -1
      return
c
 300  phl0 = wk(2)
      hl0 = h*el0
      wk(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0e0 - r*(1.0e0 - 1.0e0/wk(i+2))
        if (abs(di) .eq. 0.0e0) go to 390
 320    wk(i+2) = 1.0e0/di
 330  do 340 i = 1,n
 340    x(i) = wk(i+2)*x(i)
      return
 390  iersl = 1
      return
c
c----------------------- end of subroutine slss ------------------------
      end
      subroutine sro
     *     (n, ip, ia,ja,a, q, r, dflag)
clll. optimize
c***********************************************************************
c  sro -- symmetric reordering of sparse symmetric matrix
c***********************************************************************
c
c  description
c
c    the nonzero entries of the matrix m are assumed to be stored
c    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
c    are stored if i ne j).
c
c    sro does not rearrange the order of the rows, but does move
c    nonzeroes from one row to another to ensure that if m(i,j) will be
c    in the upper triangle of m with respect to the new ordering, then
c    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
c    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
c    stored in row j (and thus m(i,j) is not stored).
c
c
c  additional parameters
c
c    q     - integer one-dimensional work array.  dimension = n
c
c    r     - integer one-dimensional work array.  dimension = number of
c            nonzero entries in the upper triangle of m
c
c    dflag - logical variable.  if dflag = .true., then store nonzero
c            diagonal elements at the beginning of the row
c
c-----------------------------------------------------------------------
c
      integer  ip(1),  ia(1), ja(1),  q(1), r(1)
      real  a(1),  ak
c...  double precision  a(1),  ak
      logical  dflag
c
c
c--phase 1 -- find row in which to store each nonzero
c----initialize count of nonzeroes to be stored in each row
      do 1 i=1,n
  1     q(i) = 0
c
c----for each nonzero element a(j)
      do 3 i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if (jmin.gt.jmax)  go to 3
        do 2 j=jmin,jmax
c
c--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
          k = ja(j)
          if (ip(k).lt.ip(i))  ja(j) = i
          if (ip(k).ge.ip(i))  k = i
          r(j) = k
c
c--------... and increment count of nonzeroes (=q(r(j)) in that row
  2       q(k) = q(k) + 1
  3     continue
c
c
c--phase 2 -- find new ia and permutation to apply to (ja,a)
c----determine pointers to delimit rows in permuted (ja,a)
      do 4 i=1,n
        ia(i+1) = ia(i) + q(i)
  4     q(i) = ia(i+1)
c
c----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
c----for each nonzero element (in reverse order)
      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do 6 jdummy=jmin,jmax
        i = r(j)
        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  go to 5
c
c------if dflag, then put diagonal nonzero at beginning of row
          r(j) = ia(i)
          ilast = i
          go to 6
c
c------put (off-diagonal) nonzero in last unused location in row
  5       q(i) = q(i) - 1
          r(j) = q(i)
c
  6     j = j-1
c
c
c--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
      do 8 j=jmin,jmax
  7     if (r(j).eq.j)  go to 8
          k = r(j)
          r(j) = r(k)
          r(k) = k
          jak = ja(k)
          ja(k) = ja(j)
          ja(j) = jak
          ak = a(k)
          a(k) = a(j)
          a(j) = ak
          go to 7
  8     continue
c
      return
      end
      subroutine stode (neq, y, yh, nyh, yh1, ewt, savf, acor,
     1   wm, iwm, f, jac, pjac, slvs)
clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real y, yh, yh1, ewt, savf, acor, wm
      real conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnorm
      dimension neq(1), y(1), yh(nyh,1), yh1(1), ewt(1), savf(1),
     1   acor(1), wm(1), iwm(1)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c-----------------------------------------------------------------------
c stode performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stode is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the type of chord method used, or the jacobian structure.
c communication with stode is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, if a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0e0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0e0
      rc = 0.0e0
      el0 = 1.0e0
      crate = 0.7e0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/float(l)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = amin1(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = amin1(rh,abs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = amax1(rh,hmin/abs(h))
 175  rh = amin1(rh,rmax)
      rh = rh/amax1(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = amax1(0.2e0*crate,del/delp)
      dcon = del*amin1(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = amin1(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = amax1(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0e0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c----------------------- end of subroutine stode -----------------------
      end
      real function vnorm (n, v, w)
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
c-----------------------------------------------------------------------
      integer n,   i
      real v, w,   sum
      dimension v(n), w(n)
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm = sqrt(sum/float(n))
      return
c----------------------- end of function vnorm -------------------------
      end
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      integer msg, nmes, nerr, level, ni, i1, i2, nr,
     1   i, lun, lunit, mesflg, ncpw, nch, nwds
      real r1, r2
      dimension msg(nmes)
c-----------------------------------------------------------------------
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in single precision and are printed
c    in e21.13 format.
c 5. the common block /eh0001/ below is data-loaded (a machine-
c    dependent feature) with default values.
c    this block is needed for proper retention of parameters used by
c    this routine which the user can reset by calling xsetf or xsetun.
c    the variables in this block are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement below.
c
c for some systems, the data statement below must be replaced
c by a separate block data subprogram.
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh0001/ mesflg, lunit
c
c      data mesflg/1/, lunit/6/
c-----------------------------------------------------------------------
c the following data-loaded value of ncpw is valid for the cdc-6600
c and cdc-7600 computers.
c     data ncpw/10/
c the following is valid for the cray-1 computer.
c     data ncpw/8/
c the following is valid for the burroughs 6700 and 7800 computers.
c     data ncpw/6/
c the following is valid for the pdp-10 computer.
c     data ncpw/5/
c the following is valid for the vax computer with 4 bytes per integer,
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.
      data ncpw/4/
c the following is valid for the pdp-11, or vax with 2-byte integers.
c     data ncpw/2/
c-----------------------------------------------------------------------
c
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, 10) (msg(i),i=1,nwds)
c-----------------------------------------------------------------------
c the following format statement is to have the form
c 10  format(1x,mmann)
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.
c the following is valid for ncpw = 10.
c 10  format(1x,6a10)
c the following is valid for ncpw = 8.
c 10  format(1x,8a8)
c the following is valid for ncpw = 6.
c 10  format(1x,10a6)
c the following is valid for ncpw = 5.
c 10  format(1x,12a5)
c the following is valid for ncpw = 4.
  10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,e21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,e21.13,3x,4hr2 =,e21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c----------------------- end of subroutine xerrwv ----------------------
      end
c-----------------------------------------------------------------------
      real function r1mach(i)
c
c  single-precision machine constants
c
c  r1mach(1) = b**(emin-1), the smallest positive magnitude.
c
c  r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  r1mach(3) = b**(-t), the smallest relative spacing.
c
c  r1mach(4) = b**(1-t), the largest relative spacing.
c
c  r1mach(5) = log10(b)
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.
c  on rare machines a static statement may need to be added.
c  (but probably more systems prohibit it than require it.)
c
c  for ieee-arithmetic machines (binary standard), the first
c  set of constants below should be appropriate.
c
c  where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  sometimes this requires using
c  equivalent integer arrays.  if your compiler uses half-word
c  integers by default (sometimes called integer*2), you may need to
c  change integer to integer*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc
c
      real rmach(5)
c
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
c
c     machine constants for ieee arithmetic machines, such as the at&t
c     3b series, motorola 68000 based machines (e.g. sun 3 and at&t
c     pc 7300), and 8087 based micros (e.g. ibm pc and at&t 6300).
c
       data small(1) /     8388608 /
       data large(1) /  2139095039 /
       data right(1) /   864026624 /
       data diver(1) /   872415232 /
       data log10(1) /  1050288283 /, sc/987/
c
c     machine constants for amdahl machines.
c
c      data small(1) /    1048576 /
c      data large(1) / 2147483647 /
c      data right(1) /  990904320 /
c      data diver(1) / 1007681536 /
c      data log10(1) / 1091781651 /, sc/987/
c
c     machine constants for the burroughs 1700 system.
c
c      data rmach(1) / z400800000 /
c      data rmach(2) / z5ffffffff /
c      data rmach(3) / z4e9800000 /
c      data rmach(4) / z4ea800000 /
c      data rmach(5) / z500e730e8 /, sc/987/
c
c     machine constants for the burroughs 5700/6700/7700 systems.
c
c      data rmach(1) / o1771000000000000 /
c      data rmach(2) / o0777777777777777 /
c      data rmach(3) / o1311000000000000 /
c      data rmach(4) / o1301000000000000 /
c      data rmach(5) / o1157163034761675 /, sc/987/
c
c     machine constants for ftn4 on the cdc 6000/7000 series.
c
c      data rmach(1) / 00564000000000000000b /
c      data rmach(2) / 37767777777777777776b /
c      data rmach(3) / 16414000000000000000b /
c      data rmach(4) / 16424000000000000000b /
c      data rmach(5) / 17164642023241175720b /, sc/987/
c
c     machine constants for ftn5 on the cdc 6000/7000 series.
c
c      data rmach(1) / o"00564000000000000000" /
c      data rmach(2) / o"37767777777777777776" /
c      data rmach(3) / o"16414000000000000000" /
c      data rmach(4) / o"16424000000000000000" /
c      data rmach(5) / o"17164642023241175720" /, sc/987/
c
c     machine constants for convex c-1.
c
c      data rmach(1) / '00800000'x /
c      data rmach(2) / '7fffffff'x /
c      data rmach(3) / '34800000'x /
c      data rmach(4) / '35000000'x /
c      data rmach(5) / '3f9a209b'x /, sc/987/
c
c     machine constants for the cray 1, xmp, 2, and 3.
c
c      data rmach(1) / 200034000000000000000b /
c      data rmach(2) / 577767777777777777776b /
c      data rmach(3) / 377224000000000000000b /
c      data rmach(4) / 377234000000000000000b /
c      data rmach(5) / 377774642023241175720b /, sc/987/
c
c     machine constants for the data general eclipse s/200.
c
c     note - it may be appropriate to include the following line -
c     static rmach(5)
c
c      data small/20k,0/,large/77777k,177777k/
c      data right/35420k,0/,diver/36020k,0/
c      data log10/40423k,42023k/, sc/987/
c
c     machine constants for the harris slash 6 and slash 7.
c
c      data small(1),small(2) / '20000000, '00000201 /
c      data large(1),large(2) / '37777777, '00000177 /
c      data right(1),right(2) / '20000000, '00000352 /
c      data diver(1),diver(2) / '20000000, '00000353 /
c      data log10(1),log10(2) / '23210115, '00000377 /, sc/987/
c
c     machine constants for the honeywell dps 8/70 series.
c
c      data rmach(1) / o402400000000 /
c      data rmach(2) / o376777777777 /
c      data rmach(3) / o714400000000 /
c      data rmach(4) / o716400000000 /
c      data rmach(5) / o776464202324 /, sc/987/
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c      data rmach(1) / z00100000 /
c      data rmach(2) / z7fffffff /
c      data rmach(3) / z3b100000 /
c      data rmach(4) / z3c100000 /
c      data rmach(5) / z41134413 /, sc/987/
c
c     machine constants for the interdata 8/32
c     with the unix system fortran 77 compiler.
c
c     for the interdata fortran vii compiler replace
c     the z's specifying hex constants with y's.
c
c      data rmach(1) / z'00100000' /
c      data rmach(2) / z'7effffff' /
c      data rmach(3) / z'3b100000' /
c      data rmach(4) / z'3c100000' /
c      data rmach(5) / z'41134413' /, sc/987/
c
c     machine constants for the pdp-10 (ka or ki processor).
c
c      data rmach(1) / "000400000000 /
c      data rmach(2) / "377777777777 /
c      data rmach(3) / "146400000000 /
c      data rmach(4) / "147400000000 /
c      data rmach(5) / "177464202324 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers (expressed in integer and octal).
c
c      data small(1) /    8388608 /
c      data large(1) / 2147483647 /
c      data right(1) /  880803840 /
c      data diver(1) /  889192448 /
c      data log10(1) / 1067065499 /, sc/987/
c
c      data rmach(1) / o00040000000 /
c      data rmach(2) / o17777777777 /
c      data rmach(3) / o06440000000 /
c      data rmach(4) / o06500000000 /
c      data rmach(5) / o07746420233 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     16-bit integers  (expressed in integer and octal).
c
c      data small(1),small(2) /   128,     0 /
c      data large(1),large(2) / 32767,    -1 /
c      data right(1),right(2) / 13440,     0 /
c      data diver(1),diver(2) / 13568,     0 /
c      data log10(1),log10(2) / 16282,  8347 /, sc/987/
c
c      data small(1),small(2) / o000200, o000000 /
c      data large(1),large(2) / o077777, o177777 /
c      data right(1),right(2) / o032200, o000000 /
c      data diver(1),diver(2) / o032400, o000000 /
c      data log10(1),log10(2) / o037632, o020233 /, sc/987/
c
c     machine constants for the sequent balance 8000.
c
c      data small(1) / $00800000 /
c      data large(1) / $7f7fffff /
c      data right(1) / $33800000 /
c      data diver(1) / $34000000 /
c      data log10(1) / $3e9a209b /, sc/987/
c
c     machine constants for the univac 1100 series.
c
c      data rmach(1) / o000400000000 /
c      data rmach(2) / o377777777777 /
c      data rmach(3) / o146400000000 /
c      data rmach(4) / o147400000000 /
c      data rmach(5) / o177464202324 /, sc/987/
c
c     machine constants for the vax unix f77 compiler.
c
c      data small(1) /       128 /
c      data large(1) /    -32769 /
c      data right(1) /     13440 /
c      data diver(1) /     13568 /
c      data log10(1) / 547045274 /, sc/987/
c
c     machine constants for the vax-11 with
c     fortran iv-plus compiler.
c
c      data rmach(1) / z00000080 /
c      data rmach(2) / zffff7fff /
c      data rmach(3) / z00003480 /
c      data rmach(4) / z00003500 /
c      data rmach(5) / z209b3f9a /, sc/987/
c
c     machine constants for vax/vms version 2.2.
c
c      data rmach(1) /       '80'x /
c      data rmach(2) / 'ffff7fff'x /
c      data rmach(3) /     '3480'x /
c      data rmach(4) /     '3500'x /
c      data rmach(5) / '209b3f9a'x /, sc/987/
c
c  ***  issue stop 778 if all data statements are commented...
      if (sc .ne. 987) stop 778
      if (i .lt. 1  .or.  i .gt. 5) goto 999
      r1mach = rmach(i)
      return
  999 write(*,1999) i
 1999 format(' r1mach - i out of bounds',i10)
      stop
      end
c
c subroutine xsetf

      subroutine xsetf (mflag)
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsetf -----------------------
      end

      subroutine srcms (rsav, isav, job)
c-----------------------------------------------------------------------
c this routine saves or restores (depending on job) the contents of
c the common blocks ls0001, lss001, and eh0001, which are used
c internally by one or more odepack solvers.
c
c rsav = real array of length 224 or more.
c isav = integer array of length 75 or more.
c job  = flag indicating to save or restore the common blocks..
c        job  = 1 if common is to be saved (written to rsav/isav)
c        job  = 2 if common is to be restored (read from rsav/isav)
c        a call with job = 2 presumes a prior call with job = 1.
c-----------------------------------------------------------------------
      integer isav, job
      integer ieh, ils, ilss
      integer i, lenil, leniss, lenrl, lenrss
      real rsav,   rls, rlss
      dimension rsav(1), isav(1)
      common /ls0001/ rls(218), ils(39)
      common /lss001/ rlss(6), ilss(34)
      common /eh0001/ ieh(2)
      data lenrl/218/, lenil/39/, lenrss/6/, leniss/34/
c
      if (job .eq. 2) go to 100
      do 10 i = 1,lenrl
 10     rsav(i) = rls(i)
      do 15 i = 1,lenrss
 15     rsav(lenrl+i) = rlss(i)
c
      do 20 i = 1,lenil
 20     isav(i) = ils(i)
      do 25 i = 1,leniss
 25     isav(lenil+i) = ilss(i)
c
      isav(lenil+leniss+1) = ieh(1)
      isav(lenil+leniss+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrl
 110     rls(i) = rsav(i)
      do 115 i = 1,lenrss
 115     rlss(i) = rsav(lenrl+i)
c
      do 120 i = 1,lenil
 120     ils(i) = isav(i)
      do 125 i = 1,leniss
 125     ilss(i) = isav(lenil+i)
c
      ieh(1) = isav(lenil+leniss+1)
      ieh(2) = isav(lenil+leniss+2)
      return
c----------------------- end of subroutine srcms -----------------------
      end

c      block data
c      common /eh0001/ mesflg, lunit
c      common /ls0001/ rowns(209),
c     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
c     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
c     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
c     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
c     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
c      data illin/0/, ntrep/0/
c      data mesflg/1/, lunit/6/
c      end



