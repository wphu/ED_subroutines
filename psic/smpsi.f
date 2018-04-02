c  Main program implementing a simple model for plasma surface interactions

      implicit none

      real*8 Te, Ti, nei, Phi, Tp, theta, z1, z2, am1, am2, es,
     &       tgdns, gamma, Ei, Cs, qi, qe, TL, L, KA, KB, A, B,
     &       alpha, fr, KBoltz, MassP, Pii, MassE
      real*8 YieldPhyH, YieldChm, YieldRESH, YieldPhyC, YieldRESC,
     &       RnC, YieldH, YieldC, YEffGross, YEffNet, EroRateGross,
     &       EroRateNet, GammaInc, GammaIncC, GammaEvap, RnH
      integer i
      real*8 rnion
      parameter ( KBoltz=1.3804d-16, MassP=1.6726d-24,
     &            Pii=3.141592654, MassE=9.1094d-28 )

       open ( 2, file='smpsi.inp', status='old' )
       open ( 3, file='smpsi.out', status='unknown' )

c  Inputs to the subroutine
c  For a table of Surface binding energies and Densities of various
c  targets look in the "incatm" file provided with PSIC or run PSIC,
c  click on target atom and click on element of interest.
       read(2,*)Ti         ! Plasma Ion Temperature at the Plasma-
                           !    Sheath interface (eV).
       read(2,*)Te         ! Plasma Electron Temperature at the
                           !    Plasma-Sheath interface (eV).
       read(2,*)nei        ! Plasma Density at the Plasma-Sheath
                           !    interface (ions/cm^3).
       read(2,*)Phi        ! Applied bias voltage (Volts).
       read(2,*)theta      ! Angle the projectile makes with normal
                           !    to target (degrees).
       read(2,*)z1         ! Atomic number of incident atom.
       read(2,*)z2         ! Atomic number of target atom.
       read(2,*)am1        ! Atomic mass of incident atom (amu).
       read(2,*)am2        ! Atomic mass of target atom (amu).
       read(2,*)es         ! Surface binding energy (heat of
                           !    sublimation) of target (eV).
       read(2,*)tgdns      ! Density of the target (gms/cm^3).
       read(2,*)L          ! Thickness (cms) of the target.
       read(2,*)TL         ! Temperature at the opposite end (as opposed
                           ! to plasma facing end) of the target (Kelvins).
       read(2,*)KA         ! Parameter "a" in the temperature dependent
                           ! thermal conductivity profile k(T)= 1/(aT+B).
       read(2,*)KB         ! parameter "b" in the thermal conductivity
                           ! profile.
       read(2,*)A          ! Parameter in the vapor pressure curve
                           ! (see Eqn.2, and subroutine thermev for data
                           !  on A abd B for many plasma facing materials)
       read(2,*)B          ! Parameter in the vapor pressure curve
       read(2,*)alpha      ! Sticking coeff (0.6-0.9 for metals and
                           ! 0.05 for C_3 and 0.4 for C).
       read(2,*)fr         ! Returning fraction of outgoing atoms.
       read(2,*)gamma      ! (depends on the type of flow assumed in
                           !  the pre-sheath; =1 for isothermal flow,
                           !  = 5/2 for adiabatic flow with isotropic
                           !  pressure, = 3 fpr adiabatic flow) Ref.20
                           !  of paper for more detailed info.
c  Converting L into meters for subroutine heatdiff
      L = L / 100.0d0

c  Finding the sheath potential drop (See Eqn.3)
      if (Phi.eq.0.0d0) Phi = 0.5d0 * Te * log( 2.0d0 * Pii * MassE /
     &   MassP / am1 * ( 1.0d0 + Ti / Te ) )
      print*,"Sheath potential drop Phi = ",Phi," eV"

c  Bombarding ion energy, Ei (in eV) (Eqn.5)
      Ei = Te + ( 2.0d0 + gamma ) * Ti - Phi
      print*,"Bombarding ion energy = ",Ei," eV"

c  Sound Speed (in cms/sec)
      Cs = ( KBoltz / MassP * ( Te + gamma * Ti ))**0.5
      print*,"Sound Speed = ",Cs," cms/sec"

c  Incident particle flux, GammaInc (in prt/cm^2/sec)
c    Eqn.6 of paper for incident particle flux.
      GammaInc = nei * Cs
      print*,"Incident particle flux = ",GammaInc," prt/cm^2/sec"

c  Heat flux deposited by ions and electrons, qi and qe,
c  (in ergs/cm^2/sec) (Eqns.7,8)
      qi = GammaInc * Ei * 1.6022d-12
      qe = 2.0d0 * Te * GammaInc * dexp( Phi / Te ) * 1.6022d-12
c      print*,"Ion energy flux = ",qi," ergs/cm^2/sec"
c      print*,"Electron energy flux = ",qe," ergs/cm^2/sec"
      qi = qi * 1.0d-3
      qe = qe * 1.0d-3
      print*,"Ion energy flux = ",qi," Watts/m^2"
      print*,"Electron energy flux = ",qe," Watts/m^2"

c  Finding the temperature of the target plate, Tp (in Kelvins)
c  Call the heat diffusion subroutine. (Eqn.13)
      call heatdiff(qi,qe,TL,L,KA,KB,Tp)
      print*,"Surface temperature = ",Tp," K"

c  Finding the evaporated flux, GammaEvap (in atoms/cm^2/sec)
c  (Eqn.1)
      call thermev(Tp, am2, A, B, alpha, GammaEvap)
      GammaEvap = GammaEvap * 1.0d-4
      print*,"Evaporated flux = ",GammaEvap," atms/cm^2/sec"

c  Finding the gross and net erosion yields.
c    Physical sputtering yield due to H ion incident on Graphite.
      call physput(theta,Ei,z1,z2,am1,am2,es,tgdns,YieldPhyH)
      print*,"Phys sput yield due to H ion = ",YieldPhyH
c    Physical self sputtering yield.
      call physput(theta,Ei,z2,z2,am2,am2,es,tgdns,YieldPhyC)
      print*,"Phys sput yield due to C ion = ",YieldPhyC
c    RES yield due to H ion incident on Graphite.
      call ressput(Ei,z1,z2,am1,am2,es,Tp,GammaInc,YieldRESH)
      print*,"RES yield due to H ion = ",YieldRESH
c    Chemical sputtering yield of H incident on graphite
      call chmsput(Ei,z1,z2,am1,am2,Tp,GammaInc,YieldChm)
      print*,"Chem Sput yield due to H ion = ",YieldChm
c    Backscattering of C incident on C
      RnC = rnion(theta,Ei,nint(z2),nint(am2),1,nint(am2),1)
      print*,"Backscattering coeff of C ion = ",RnC

c    Backscattering of H incident on target
      RnH = rnion(theta,Ei,nint(z1),nint(am1),1,nint(am2),1)

c    Total yield due to H incidence
      YieldH = YieldPhyH + YieldRESH + YieldChm

c    RES yield due to self sputtering.
      GammaIncC = fr * (YieldH * GammaInc + GammaEvap) ! roughly ...
      call ressput(Ei,z2,z2,am2,am2,es,Tp,GammaIncC,YieldRESC)
      print*,"RES yield due to C ion = ",YieldRESC

c    Total yield due to C incidence
      YieldC = YieldPhyC + YieldRESC

c    Efective gross erosion yield, (Eqn.19)
      YEffGross = ( YieldH + GammaEvap/GammaInc ) * (1.0d0 - fr*RnC )
     &     / (1.0d0 - fr * ( YieldC + RnC ))
c    Efective net erosion yield, (Eqn.20)
      YEffNet = ( YieldH + GammaEvap/GammaInc ) * ( 1.0d0 - fr )
     &     / (1.0d0 - fr * ( YieldC + RnC ))
      print*,"Effective Gross erosion yield = ",YEffGross
      print*,"Effective Net erosion yield = ",YEffNet

c  Erosion rate of target (Eqn.21)
c    Gross
      EroRateGross = YEffGross*GammaInc*am2/tgdns/6.0221d23
c    Net
      EroRateNet = YEffNet*GammaInc*am2/tgdns/6.0221d23
      print*,"Effective Gross erosion rate = ",EroRateGross," cms/sec"
      print*,"Effective Net erosion rate = ",EroRateNet," cms/sec"

      stop
      end
