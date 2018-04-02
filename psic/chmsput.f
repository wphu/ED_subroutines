c  SUBROUTINE  chmsput

c  Reference: J. Roth, "Chemical erosion of carbon-based materials in
              fusion devices", J. Nucl. Mater., Vol.266-269 (1999) 51-57.

c  Inputs to the subroutine :
c  --------------------------
c  eo -> energy of incident projectile (eV)
c  z1	-> nuclear charge of incident atom.
c  z2	-> nuclear charge of target atom.
c  am1	-> atomic mass of incident atom (amu).
c  am2	-> atomic mass of target atom (amu).
c  tp	-> temperature of the target plate (K).
c  gamma-> incident ion flux (/cm**2/sec).

c  Outputs from the subroutine :
c  -----------------------------
c  yldchm -> Chemical sputtering yeild.

c  Make the following change yourselves (ive used erel = 1.8 eV).
c  erel = 1.8 eV for pure carbon,
c       = 1.5 eV for Si, Ti, W doped carbon, and
c       = 1.2 eV for B doped carbon.

c  Subroutine to evaluate sputtering yeilds. 

      subroutine chmsput(eo,z1,z2,am1,am2,tp,gamma,yldchm)

      implicit real*8(a-h,o-z)

      real*8 eo, z1, z2, am1, am2, tp, gamma, etf,
     @    stoppwr, yeildpar, yldchm
 
c  Input checks
      if (eo.lt.0.0) then
        write(*,*)"Input error; Correct so that 0 < eo"
        return
      endif
      if (z1.lt.1.0d0.or.z1.ge.2.0d0) then
        write(*,*)"Input error z1 = ", z1
        write(*,*)"chmsput valid only for hydrogen isotope incident"
        write(*,*)"on graphite; Correct so that 1.0 <= z1 < 2.0"
        return
      endif
      if (am1.lt.1.0d0.or.am1.ge.4.0d0) then
        write(*,*)"Input error"
        write(*,*)"chmsput valid only for hydrogen isotope incident"
        write(*,*)"on graphite; Correct so that 1.0 <= am1 < 4.0"
        return
      endif
      if (z2.lt.6.0d0.or.z2.ge.7.0d0) then
        write(*,*)"Input error"
        write(*,*)"chmsput valid only for graphite target"
        write(*,*)"Correct so that 6.0 <= z2 < 7.0"
        return
      endif
      if (am2.lt.12.0d0.or.am2.ge.13.0d0) then
        write(*,*)"Input error"
        write(*,*)"chmsput valid only for graphite target"
        write(*,*)"Correct so that 12.0 <= am2 < 13.0"
        return
      endif
      if (tp.gt.3000.0d0) then
        write(*,*)"Input error; Correct so that tp <= 3000 K"
        return
      endif

c  Initialising
      pwr2by3 = 2.0d0 / 3.0d0
      z123z223sum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

c  The Thomas Fermi energy etf.
      etf = 30.74d0 * (am1+am2)/am2  * z1*z2*z123z223sum**0.5d0

c  Converting target temperature from Kelvin to eV
      tpev = tp / 1.1604d4
	
      eobyetf = eo / etf
      stoppwr1 = (0.5d0*dlog(1.0d0+1.2288d0*eobyetf))
      stoppwr2 = (eobyetf + 0.1728d0*eobyetf**0.5d0 +
     @    0.008d0*eobyetf**0.1504d0)
      stoppwr = stoppwr1 / stoppwr2

      flxm2 = gamma*1.0d4
      flxlmt = 1.0d30*exp(-1.4/tpev)
      if (flxm2.gt.flxlmt) then  ! Assuming high flux
          C = 1.0d0/(1.0d0+3.0d-23*flxm2)
      else
          C = 1.0d0/(1+3.0d7*dexp(-1.4/tpev))
      endif
      etherm = 1.7d0
      erel = 1.8d0
      edam = 15.0d0
      edes = 2.0d0
      Csp3num = C*(2.0d-32*flxm2+dexp(-etherm/tpev))
      Csp3denom = 2.0d-32*flxm2+
     @    (1.0d0+2.0d29/flxm2*dexp(-erel/tpev))*dexp(-etherm/tpev)
      Csp3 = Csp3num/Csp3denom
      ytherm = Csp3*0.033*dexp(-etherm/tpev) /
     @  (2.0d-32*flxm2+dexp(-etherm/tpev))
      if (nint(am1).eq.1) then
          yeildpar = 0.035
      elseif (nint(am1).eq.2) then
          yeildpar = 0.1
      elseif (nint(am1).eq.3) then
          yeildpar = 0.12
      endif
      D = 250.0d0/am1
      ydam = yeildpar*stoppwr*(1.0d0-(edam/eo)**pwr2by3)*
     @    (1.0d0-edam/eo)**2.0d0
      ydes = yeildpar*stoppwr*(1.0d0-(edes/eo)**pwr2by3)*
     @    (1.0d0-edes/eo)**2.0d0
      ysurf = csp3*ydes/(1.0d0+dexp((eo-65.0d0)/40.0d0))
      yldthmprt = ytherm*(1.0d0+D*ydam)
      yldchm = yldthmprt+ysurf

      if (yldchm.le.0.0d0) yldchm = 0.0d0

      return
      end
