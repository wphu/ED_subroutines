C   Subroutine to evaluate the thermal evaporation of Graphite.

c   Ref: 1) C. Garc'ia Rosales, Erosion processes in plasma-wall
c           interactions, J. Nucl. Mater. 211 (1994) Pg.202-214.
c        2) R.A. Langley, Evaporation, Data Compendium for Plasma
c           Surface Interactions, Special Issue, NUCLEAR FUSION,
c           IAEA, Vienna, (1984) 55-61.

c   Written by Manoj Warrier on 17-05-96.

*********************************************************************************

c  This routine takes the temperature of the graphite plate ( tp ), and the
c  atomic mass (am2 ) of the target as inputs and outputs the flux of carbon
c  atoms evaporated from the graphite plate ( flxthev ) per meter**2 per sec.
c
c  Table of A and B for some first wall materials in Fusion research:
c  (Ref (2) above).
c  (note that A and B are provided such that the vappres is in Pascals.
c  element    |  A       |  B      |
c  Aluminum   | -16165   | 10.91   |
c  Beryllium  | -16720   | 11.61   |
c  Graphite   | -40181   | 14.80   |
c  Copper     | -17079   | 11.26   |
c  Molybdenum | -33069   | 12.02   |
c  Titanium   | -23340   | 9.79    |
c  Tungsten   | -44485   | 12.74   |
c  Titanium - |          |         |
c     Carbide | -39117.9 | 15.23   |
c  Stainless -|          |         |
c     Steel   | -14321.9 | 9.44    |

	subroutine thermev(tp,am2,A,B,alpha,flxthev)

	implicit real*8(a-h,o-z)
	real*8 tp, am2, flxthev, vappres, alpha

c Inputs:
c  tp       ->	temperature of the Graphite plate in Kelvin.
c  am2      ->	atomic mass of the target ( in amu )
c  A        ->  Parameter in typical Vapor pressure curve.
c  B        ->  Parameter in typical Vapor pressure curve.
c  alpha    ->	Sticking Coefficient ( = 0.6-0.9 for metals and
c               0.05 for C_3 and 0.4 for C )

c Output:
c  flxthev  ->	Flux of target atoms evaporated ( /m**2/sec )

c  vappres  ->	Vapour Pressure of the target in Pascals.

c  Checking inputs
      if (tp.gt.3000.0d0) then
        write(*,*)"Input error; Correct so that tp <= 3000 K"
        return
      endif

      vappres = 10.0d0**(A/tp + B)
      flxthev = alpha * 2.62521d24 * vappres / dsqrt(am2 * tp)
      if (flxthev.le.0.0) flxthev = 0.0

      return
      end
