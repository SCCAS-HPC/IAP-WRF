module getsmark
contains
subroutine GetSMarkChange(DeltSpeMark,OrgCon,TmpSM, &
                          MapS,ISrcDefined,IHgtLev,ismMax)  
real DeltSpeMark,OrgCon
real TmpSM(ismMax)
real,allocatable,dimension(:) ::  NewCon
integer MapS

allocate(NewCon(ismMax))
do i=1,ismMax
 NewCon(i)=OrgCon*TmpSM(i)
enddo
 ! 1 initial condition
 ! 2 Stratopshere
 ! 3 Boudary
 ! definde source * 3 + 3
 if(MapS.gt.ISrcDefined) then
    print *,'MapS has problem',MapS,ISrcDefined
    stop 
 endif
 ic= 3 + IHgtLev*ISrcDefined + MapS
 NewCon(ic)=NewCon(ic)+DeltSpeMark
! ConNew = DeltSpeMark+OrgCon 
  ConNew=0
  do i=1,ismMax
     ConNew=ConNew+NewCon(i)
  enddo
do i=1,ismMax
 TmpSM(i)=NewCon(i)/ConNew
enddo
deallocate(NewCon)
return
end subroutine 

subroutine SMmixing(ConOrg,TmpSMthis,deltc,TmpSMother,ismMax)
integer ismMax,is
real ConOrg, TmpSmthis(ismMax),deltc,TmpSMother(ismMax)
!real,allocatable,dimension(:) ::  NewCon
real cc
real NewCon(100)

!allocate(NewCon(ismMax))
!if( TmpSmthis(7).gt.0)then
!    print *, TmpSmthis(7),'ddddd',TmpSMother(7)
!    print *,'---------------------',Conorg,deltc
!   endif
cc=0
ConOrg = MAX(ConOrg, 1.E-20)

IF(deltc>1.E-5) THEN
do is=1,ismMax
  NewCon(is)=ConOrg*TmpSMthis(is)+deltc*TmpSMother(is)
  cc=cc+NewCon(is)
enddo

do is=1,ismMax
  if(cc.gt.0)TmpSMthis(is)=NewCon(is)/cc
  !TmpSMthis(is)=TmpSMthis(is)
enddo

ENDIF
!if( TmpSmthis(7).gt.0)print *, TmpSmthis(7),'eeeee',TmpSMother(7)

!deallocate(NewCon)
return
end subroutine

subroutine GetBounChange(DeltSpeMark,OrgCon,TmpSM, &
                         Iset,ismMax)
real DeltSpeMark,OrgCon
real TmpSM(ismMax)
real,allocatable,dimension(:) ::  NewCon

allocate(NewCon(ismMax))

do i=1,ismMax
 NewCon(i)=OrgCon*TmpSM(i)
enddo
 ! 3 initial condition
 ! 2 Stratopshere
 ! 1  Boudary
 ic= Iset

 NewCon(ic)=NewCon(ic)+DeltSpeMark
! ConNew = DeltSpeMark+OrgCon
  ConNew=0
  do i=1,ismMax
     ConNew=ConNew+NewCon(i)
  enddo
do i=1,ismMax
 TmpSM(i)=NewCon(i)/ConNew
enddo
deallocate(NewCon)
return
end subroutine

end module getsmark
