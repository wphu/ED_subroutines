c  Subroutine to evaluate sputtering yeilds of any monoatomic
c  target material due to Physical Sputtering.

c  Reference: "Revised formulae for sputtering data", C. Garc\'ia
c              Rosales, W. Eckstein, J. Roth, J. Nucl. Mat. 218
c              (1994) 8-17.

      subroutine physput(theta,eo,z1,z2,am1,am2,es,tgdns,yldphy)
      implicit real*8(a-h,o-z)

      real*8 eo, z1, z2, am1, am2, es, eth, etf, E_L,
     @stoppwr, yeildpar, ethbyeo, yldphy
 
c  Inputs to the subroutine
c  theta -> angle (in degrees) with normal to target of incident particle.
c  eo	-> energy of the incident particle (eV)
c  z1	-> nuclear charge of incident atom.
c  z2	-> nuclear charge of target atom.
c  am1	-> atomic mass of incident atom (amu).
c  am2	-> atomic mass of target atom (amu).
c  es	-> surface binding energy (heat of sublimation) of target (eV).
c  tgdns -> target density (gms/cc)

c  Output from the subroutine
c  yldphy -> The physical sputtering yield.

c  Checking inputs
      if (theta.lt.0.0.or.theta.ge.90.0d0) then
        write(*,*)"Input error; Correct so that 0 <= theta < 90"
        return
      endif
      if (eo.lt.0.0) then
        write(*,*)"Input error; Correct so that 0 < eo"
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

      ypar1 = 1.633 * z123xz223 * z123z223sum**pwr1by3 /
     @    es**pwr2by3
      ypar2 = am1**pwr5by6 * am2**pwr1by6 / (am1 + am2)
      ypar3 = (0.15d0 + 0.05d0*am2byam1) /
     @    (1.0d0 + 0.05d0*am2byam1**1.6d0)
      yeildpar = ypar1 * ypar2 * ypar3

      eth = ( 7.0d0/am2byam1**0.54d0 + 0.15d0*am2byam1**1.12d0) * es
      ethbyeo = eth / eo

      etf = 30.736d0 * (am1+am2)/am2  * z1*z2*z123z223sum**0.5d0

      eobyetf = eo / etf
      stoppwr1 = (0.5d0*dlog(1.0d0+1.2288d0*eobyetf))
      stoppwr2 = (eobyetf + 0.1728d0*eobyetf**0.5d0 +
     @    0.008d0*eobyetf**0.1504d0)
      stoppwr = stoppwr1 / stoppwr2

c  iionflg -> flag for light/heavy ion.
c  iionflg = 0 => light ion sputtering.
c  iionflg = 1 => heavy ion sputtering.
      if (am1.le.4.0) then
          iionflg = 0
      else
          iionflg = 1
      endif

      Ro = (am2/tgdns/6.0221d23)**(1.0d0/3.0d0)*1.0d8
      a = 0.4685d0*dsqrt(1.0d0/z123z223sum)
      if (iionflg.eq.0) then
          f = (0.94 - 1.33d-3 * am2byam1) * dsqrt(es)
          gam = 4.0d0*am1*am2/(am1+am2)**2
          q = dsqrt(es/gam/eo)
          anu = a/Ro*dsqrt(1.0d0/2.0d0/eobyetf/q)
          ang_opt = 90.0d0 - 57.29*anu
      else
c       For a really good quantitative study you have to
c       input the value of fs from the above reference. I am
c       circumventing the need of making a big database
c       by using fig-3 from the following reference:
c       Institute of Plasma Physics, Nagoya University
c       report - IPPJ-AM-26, Japan.
          if (am2byam1.le.2.0d0) fs = 1.75d0
          if (2.0d0.le.am2byam1.and.am2byam1.lt.4.0d0) fs = 1.5d0
          if (4.0d0.le.am2byam1.and.am2byam1.lt.6.0d0) fs = 1.25d0
          if (6.0d0.le.am2byam1.and.am2byam1.lt.8.0d0) fs = 1.0d0
          if (am2byam1.ge.8.0d0) fs = 0.8d0
c       The if .. then ... else below is introduced to avoid
c       the blowing up of f when eta = 0.0d0 (Xavier Bonnin suggested
c       this during inclusion of this routine in SOLPS-B2.5).   
          eta = 1.0d0 - dsqrt(eth/eo)
          if (eta.gt.0.0d0) then
              f = min(10.0d0,fs*(1.0d0+2.5d0*(1-eta)/eta))
          else
              f = fs
          endif
          psi = (a/Ro)**(3.0d0/2.0d0)*dsqrt(z1*z2/dsqrt(z123z223sum)/eo)
          ang_opt = 90.0d0-286.0*psi**0.45
      endif
      Sigma = dcos(ang_opt*1.7453293d-2)
      t = 1.0d0/dcos(theta*1.7453293d-2)
      angcntrb = dexp(f*(Sigma*(1-t)+dlog(t)))

c  Calculation of the physical sputtering yeild
      yldphy = yeildpar * stoppwr * (1.0d0-ethbyeo**pwr2by3)
     @    * ( 1.0d0 - ethbyeo ) ** 2 * angcntrb
      if (yldphy.lt.0.0d0) yldphy = 0.0d0

      return
      end
