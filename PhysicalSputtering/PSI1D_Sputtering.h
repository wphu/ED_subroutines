/* ==============================================================
Subroutine to evaluate sputtering yeilds of any
monoatomic target material due to Physical Sputtering
Ref: Subroutines for some plasma surface interaction processes:
     hpysical sputtering, chemical erosion, radiation enhanced
     sublimation, backscattering and thermal evaporation.

!!! 戴舒宇和桑超峰的程序，靶板密度的单位不一样，桑超峰的是g/cm^3，戴书宇的是atoms/Ai^3 （Ai = 10^-10 m）
!!! 此处是改自psic/physput.f
================================================================*/



#ifndef PSI1D_SPUTTERING_H
#define PSI1D_SPUTTERING_H

#include <vector>



class PSI1D_Sputtering
{

public:
    //! Constructor for Collisions between two species
    PSI1D_Sputtering(
    double an1_in,
    double am1_in,
    double an2_in,
    double am2_in,
    double es_in,
    double density_in );


    ~PSI1D_Sputtering();


    // parameters for sputtering
    double an1;     // atomic number of incident atomic
    double an2;     // atomic number of target atomic
    double am1;     // atomic mass of incident atomic (amu)
    double am2;     // atomic mass of target atomic (amu)
    double es;      //surface binding energy (heat of sublimation) of target (eV).
                    // W: 11.75(old value 8.68) C: 7.41
    double density; // target density, unit: g/cm^3 （Ai = 10^-10 m）
                    // C: 2.248 g/cm^3, W: 19.35 g/cm^3
    double n;       // target number density, unit: atoms/Ai^3 （Ai = 10^-10 m)  n = 1.0e-24 * density * NA / am2
    double Ro;      // n^(-1/3)
    int ionflag;    //ionflag -> flag for light/heavy ion.
                    //ionflag = 0 => light ion sputtering.
                    //ionflag = 1 => heavy ion sputtering.
    double Q, eth,eth1, mu, etf, aL, Mratio;

    void init();
    double phy_sput_yield(double ke, double theta);

private:


};


#endif
