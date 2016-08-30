module chemprod_mark
contains
    subroutine chemprod(ig,i,j,k,delta,dt,ne)  ! calculate the production of species
      include 'chm1.inc'
      include 'gas1.inc'
        integer :: ig,i,j,k,ne
        real :: delta,delta0,dt
!ozone and its the production occurs in daytime         
!  F(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
     delta=0.0
   IF(ig==11.and.rk_photo(jphoto_no2)>0.0) THEN
     delta0=60.*dt*(rk_com(33)*cnn(kno)*cnn(kho2)& !k4[NO][HO2] in molec/cc
              +rk_com(57)*cnn(kch3o2)*cnn(kno)& !k5[NO][CH3O2]) 
              +rk_urb(17)*cnn(kto2)*cnn(kno)&   !k[TO2][NO]
              +rk_com(58)*cnn(kethp)*cnn(kno)&  !K[ETHP][NO] 
              +rk_urb(28)*cnn(kro2)*cnn(kno)&   !K[RO2][NO]
              +rk_com(71)*cnn(kc2o3)*cnn(kno)&  !K[c2o3][no] 
              +rk_urb(29)*cnn(kano2)*cnn(kno)&  !k[ano2][no]
              +rk_urb(30)*cnn(knap)*cnn(kno)&  !k[nap][no]
              +rk_urb(31)*cnn(kxo2)*cnn(kno)&  !k[xo2][no]
              +rk_bio(8)*cnn(kisopp)*cnn(kno)& !k[isopp][no]
              +rk_bio(9)*cnn(kisopn)*cnn(kno)& !k[isopn][no]
              +rk_bio(10)*cnn(kisopo2)*cnn(kno)) !k[isopo2][no]
!-----------convert molec/cc to ppbv-----
    delta = delta0/cair_mlc*1.e+9     
   ENDIF
!      if(i==46.and.j==30.and.k==3.and.ig==11.and.ne==1)   print*,delta,'555'

    return
    end subroutine
   
    subroutine chemprodloss(myid,delta1,delta2,delta3,delta4,dt,i,j,k,ne) 
      include 'chm1.inc'
      include 'gas1.inc'
      integer myid,i,j,k,ne
      real delta1,delta2,delta3,delta4,dt
!ozone production and loss,and loss due to nox,radicals from radical-radical reactions
!F(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
!L(o3)=k3[O1D][H2O]+k8[O3][HO2]+k7[O3][OH]+(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
!Ln(radical)=k[no2][oh]
!Lr(ridical)=k[RO2][HO2]+2k[HO2]^2
! Ln/(Ln+Lr)<0.4 is nox limited, >0,6 VOCs limited
! refrence :Kleinman, L. I., and Coauthors, 1997: Dependence of O3 production on
! NO and hydrocarbons in the troposphere. Geophys. Res. Lett., 24, 2299-2302.
delta1=60.*dt*(rk_com(33)*cnn(kno)*cnn(kho2)& !k4[NO][HO2] in molec/cc
            +rk_com(57)*cnn(kch3o2)*cnn(kno)& !k5[NO][CH3O2])
            +rk_urb(17)*cnn(kto2)*cnn(kno)&   !k[TO2][NO]
            +rk_com(58)*cnn(kethp)*cnn(kno)&  !K[ETHP][NO]
            +rk_urb(28)*cnn(kro2)*cnn(kno)&   !K[RO2][NO]
            +rk_com(71)*cnn(kc2o3)*cnn(kno)&  !K[c2o3][no]
            +rk_urb(29)*cnn(kano2)*cnn(kno)&  !k[ano2][no]
            +rk_urb(30)*cnn(knap)*cnn(kno)&   !k[nap][no]
            +rk_urb(31)*cnn(kxo2)*cnn(kno)&   !k[xo2][no]
            +rk_bio(8)*cnn(kisopp)*cnn(kno)&  !k[isopp][no]
            +rk_bio(9)*cnn(kisopn)*cnn(kno)&  !k[isopn][no]
            +rk_bio(10)*cnn(kisopo2)*cnn(kno)&!k[isopo2][no]
            )/cair_mlc*1.e+9

!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta1,'111'

      delta2= 60.*dt*(rk_com(12)*cnn(ko1d)*h2o& ! k3[O1D][H2O]
            +rk_com(21)*cnn(ko3)*cnn(kho2)&  ! k8[O3][HO2]
            +rk_com(20)*cnn(ko3)*cnn(koh)&   ! k7[O3][OH]
            +(rk_com(18)*cnn(kno)*cnn(ko3)*& !k10[NO][O3]
             rk_com(24)*cnn(kno2)*cnn(koh)/& !*k12[NO2][OH]
            (rk_com(1)*cnn(kno2)+rk_com(24)*cnn(kno2)*cnn(koh))&! k11[NO2]+k12[NO2][OH]
            )&      !(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
            )/cair_mlc*1.e+9            
            
!             if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta2,'2222'

     delta3 = 60.*dt*(rk_com(24)*cnn(kno2)*cnn(koh))/cair_mlc*1.e+9 !k[no2][oh]
!            if(i==46.and.j==30.and.k==3.and.ne==1) print*,delta3,'333'
     
     delta4 = 60.*dt*(rk_com(61)*cnn(kch3o2)*cnn(kho2)& !k[HO2][ch3o2]
             +rk_com(62)*cnn(kethp)*cnn(kho2)&!k[HO2][ETHP]
             +rk_urb(36)*cnn(kro2)*cnn(kho2)&!k[HO2][RO2]
             +rk_com(73)*cnn(kc2o3)*cnn(kho2)&!k[HO2][C2O3]
             +rk_urb(37)*cnn(kano2)*cnn(kho2)&!k[HO2][ANO2]
             +rk_urb(38)*cnn(knap)*cnn(kho2)&!k[HO2][NAP]
             +rk_bio(11)*cnn(kisopp)*cnn(kho2)&!k[HO2][ISOPP]
             +rk_bio(12)*cnn(kisopn)*cnn(kho2)&!k[HO2][ISOPN]
             +rk_bio(13)*cnn(kisopo2)*cnn(kho2)&!k[HO2][ISOPO2]
             +rk_urb(39)*cnn(kxo2)*cnn(kho2)&!k[HO2][XO2]
             +2.*rk_com(31)*cnn(kho2)*cnn(kho2)& !2k[HO2]^2
             )/cair_mlc*1.e+9

       d1=rk_com(61)*cnn(kch3o2)*cnn(kho2)
       d2=rk_com(62)*cnn(kethp)*cnn(kho2)
       d3=rk_urb(36)*cnn(kro2)*cnn(kho2)
       d4=rk_com(73)*cnn(kc2o3)*cnn(kho2)
       d5=rk_urb(37)*cnn(kano2)*cnn(kho2)
       d6=rk_urb(38)*cnn(knap)*cnn(kho2)
       d7=rk_bio(11)*cnn(kisopp)*cnn(kho2)
       d8=rk_bio(12)*cnn(kisopn)*cnn(kho2)
       d9=rk_bio(13)*cnn(kisopo2)*cnn(kho2)
       d10=rk_urb(39)*cnn(kxo2)*cnn(kho2)
       d11=2.*rk_com(31)*cnn(kho2)*cnn(kho2)

             if(i==46.and.j==30.and.k==3.and.ne==1)&
     !        print*,d1,d2,d3,d4, '111'
     !        print*,d5,d6,d7,d8, '222'
     !        print*,d9,d10,d11 ,'333' 
     !         print*,delta3,delta4,60.*dt*d11/cair_mlc*1.e+9,60.*dt*(d1+d2+d3+d4)/cair_mlc*1.e+9,60.*dt*(d5+d6+d7+d8)/cair_mlc*1.e+9,&
     !                  60.*dt*(d9+d10)/cair_mlc*1.e+9
     return
     end subroutine
end module chemprod_mark
