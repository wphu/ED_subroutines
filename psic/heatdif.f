c  Subroutine to find the surface temperature of a plasma facing
c  surface assuming a 1-D. steady state, heat diffusion with a
c  thermal conductivity that varies as 1/AT+B

c  Reference: R. Behrisch, "Contribution of the different erosion
c             processes to material release from the vessel walls
c             of fusion devices during plasma operation", Contrib.
c             Plasma Phys. Vol. 42 (2002) 431-444

      subroutine heatdiff(qi,qe,TL,L,KA,KB,Tp)

      implicit real*8(a-h,o-z)

      real*8 qi, qe, TL, L, KA, KB, Tp

c  Inputs:
c  qi -> incident ion heat flux (joules/m^2/sec)
c  qe -> incident electron heat flux (joules/m^2/sec)
c  TL -> Temperature in K at x=L (x=0 is plasma facing surface,
c                            x=L is the cold end).
c  L  -> Thickness of the target plate (m).
c  KA -> parameter in the expression for conductivity (see above)
c  KB -> parameter in the expression for conductivity (see above)

c  Outputs:
c  Tp -> Surface temperature (at x=0) (in Kelvins).

      Tp = (TL + KB/KA) * dexp(KA*(qi+qe)*L) - KB/KA

      if (Tp.ge.3000) Tp = 3000.0d0
      if (Tp.le.TL) Tp = TL

      return
      end      
