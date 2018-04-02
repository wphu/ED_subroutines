c  Subroutine to evaluate sputtering yields of
c  Graphite due to Radiation Enhanced Sublimation.

c  Reference: J. Roth, E. Vietzke and A.A. Haasz, "Erosion of graphite",
              Supplement to the journal Nuclear Fusion, Vol. 1 (1991) 63-78.

      subroutine ressput(eo,z1,z2,am1,am2,es,tp,gamma,yldres)

      implicit real*8(a-h,o-z)

      real*8 eo, z1, z2, am1, am2, es, tp, gamma, eth, etf,
     @stoppwr, ethbyeo, yldparres
 
c  Inputs to the code

c  eo	-> incident ion energy (in eV).
c  z1	-> nuclear charge of incident atom.
c  z2	-> nuclear charge of target atom.
c  am1	-> atomic mass of incident atom (amu).
c  am2	-> atomic mass of target atom (amu).
c  es	-> surface binding energy (heat of sublimation) of target (eV).
c  tp	-> temperature of the target plate (eV).
c  gamma -> incident ion flux (/cm^2/sec).
                                                                                
c  Input checks
      yldres = 0.0
      if (eo.lt.0.0) then
        write(*,*)"Input error; Correct so that 0 < eo"
        return
      endif
      if (z2.lt.6.0d0.or.z2.ge.7.0d0) then
        write(*,*)"Input error"
        write(*,*)"ressput valid only for graphite target"
        write(*,*)"Correct so that 6.0 <= z2 < 7.0"
        return
      endif
      if (am2.lt.12.0d0.or.am2.ge.13.0d0) then
        write(*,*)"Input error"
        write(*,*)"ressput valid only for graphite target"
        write(*,*)"Correct so that 12.0 <= am2 < 13.0"
        return
      endif
      if (tp.gt.3000.0d0) then
        write(*,*)"Input error; Correct so that tp <= 3000 K"
        return
      endif
      if (gamma.le.0.0) then
        write(*,*)"Input error; Correct so that gamma > 0.0"
        return
      endif

c  Initialising 
      pwr1by6 = 1.0d0 / 6.0d0
      pwr5by6 = 5.0d0 / 6.0d0
      pwr1by3 = 1.0d0 / 3.0d0
      pwr2by3 = 2.0d0 / 3.0d0

      z123xz223 = z1**pwr2by3 * z2**pwr2by3
      z123z223sum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

      eth = ( 7.0d0/am2byam1**0.54d0 + 0.15d0*am2byam1**1.12d0) * es

      etf = 30.74d0 * (am1+am2)/am2  * z1*z2*z123z223sum**0.5d0
      tpev = tp / 11600.0d0

      yldparres = 54.0d0 * am1**1.18d0 * dexp((-0.78d0)/tpev)
     @  * (gamma/1.0d16)**(-0.1)
	
      ethbyeo = eth / eo
      eobyetf = eo / etf
      stoppwr1 = (0.5d0*dlog(1.0d0+1.2288d0*eobyetf))
      stoppwr2 = (eobyetf + 0.1728d0*eobyetf**0.5d0 +
     @  0.008d0*eobyetf**0.1504d0)
      stoppwr = stoppwr1 / stoppwr2

c  Calculation of the RES sputtering yield

      yldres = yldparres * stoppwr * (1.0d0-ethbyeo**pwr2by3)
     @  * ( 1.0d0 - ethbyeo ) ** 2

      if (yldres.le.0.0d0) yldres = 0.0d0

      return
      end
