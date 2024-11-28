//////////////////////////
// background.hpp
//////////////////////////
// 
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: September 2018
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include <gsl/gsl_integration.h>


double FermiDiracIntegrand(double q, void * w)
{
	return q * q * sqrt(q * q + *(double *)w) / (exp(q) + 1.0l);
}

//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
// 
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
// 
//////////////////////////

double FermiDiracIntegral(double &w)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	
	f.function = &FermiDiracIntegrand;
	f.params = &w;
		
	gsl_integration_qng(&f, 0.0l, 24.0l, 5.0e-7, 1.0e-7, &result, &err, &n);
	
	return result;
}


//////////////////////////
// bg_ncdm (1)
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//
// Returns: value for the background model
// 
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo, const int p)
{
	if (p < 0 || p >= cosmo.num_ncdm)
		return 0;
	else
	{
		double w = a * cosmo.m_ncdm[p] / (pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST);
		w *= w;
		
		return FermiDiracIntegral(w) * cosmo.Omega_ncdm[p] * pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p] / C_FD_NORM / a;
	}
}


//////////////////////////
// bg_ncdm (2)
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such that
//   multiple calls at the same value of a will not result in multiple integrations
//   being carried out. This assumes that the cosmological model should not change!
//
// Returns: value for the background model
// 
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo)
{
	double w;
	static double result = -1.0;
	static double a_prev = -1.0;
	
	if (a != a_prev)
	{
		result = 0.0;
		a_prev = a;
		
		for (int p = 0; p < cosmo.num_ncdm; p++)
			result += bg_ncdm(a, cosmo, p);
	}
	
	return result;
}


// Default Omega_Lambda_now functor
struct OmegaLambdaNowDefault
{
    inline double operator()(double a, const cosmology &cosmo) const
    {
        if (a < 0.341672)
            return -cosmo.Omega_Lambda_0;
        else
            return cosmo.Omega_Lambda_0;
    }
};

// Omega_Lambda_now functor for before the discontinuity (always negative)
struct OmegaLambdaNegative
{
    inline double operator()(double a, const cosmology &cosmo) const
    {
        return -cosmo.Omega_Lambda_0;
    }
};

// Omega_Lambda_now functor for after the discontinuity (always positive)
struct OmegaLambdaPositive
{
    inline double operator()(double a, const cosmology &cosmo) const
    {
        return cosmo.Omega_Lambda_0;
    }
};

double Omega_Lambda_now(const double a, const cosmology cosmo) {
    return OmegaLambdaNowDefault()(a, cosmo);
}

//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
// 
//////////////////////////

template <typename TOmegaLambdaFunc = OmegaLambdaNowDefault>
double Hconf(const double a, const double fourpiG, const cosmology cosmo, TOmegaLambdaFunc omegaLambdaNowFunc = TOmegaLambdaFunc())
{
    double Omega_Lambda_eff = omegaLambdaNowFunc(a, cosmo);
	return   sqrt((2. * fourpiG / 3.)
           * (  ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a)
              + (Omega_Lambda_eff * a * a)
              + (cosmo.Omega_rad / a / a)
              + (cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 1. + 3. * (cosmo.w0_fld + cosmo.wa_fld)))
             ));
}


double Omega_m(const double a, const cosmology cosmo) { return cosmo.Omega_m / (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + Omega_Lambda_now(a, cosmo) * a * a * a + cosmo.Omega_rad / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld))); }

double Omega_rad(const double a, const cosmology cosmo) { return (cosmo.Omega_rad + (bg_ncdm(a, cosmo) + cosmo.Omega_cdm + cosmo.Omega_b - cosmo.Omega_m) * a) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + Omega_Lambda_now(a, cosmo) * a * a * a * a + cosmo.Omega_rad + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld) - 1.)); }

double Omega_Lambda(const double a, const cosmology cosmo) { return Omega_Lambda_now(a, cosmo) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a / a / a + Omega_Lambda_now(a, cosmo) + cosmo.Omega_rad / a / a / a / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. + 3. * (cosmo.w0_fld + cosmo.wa_fld))); }


//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a fourth-order
//   Runge-Kutta method
// 
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
// 
//////////////////////////


// Template for RK4 step
template <typename TOmegaLambdaFunc>
double rk4_step(double a, const double fourpiG, const cosmology &cosmo, const double dtau)
{
    double k1a, k2a, k3a, k4a;

    k1a = a * Hconf<TOmegaLambdaFunc>(a, fourpiG, cosmo);
    k2a = (a + k1a * dtau / 2.) * Hconf<TOmegaLambdaFunc>(a + k1a * dtau / 2., fourpiG, cosmo);
    k3a = (a + k2a * dtau / 2.) * Hconf<TOmegaLambdaFunc>(a + k2a * dtau / 2., fourpiG, cosmo);
    k4a = (a + k3a * dtau) * Hconf<TOmegaLambdaFunc>(a + k3a * dtau, fourpiG, cosmo);

    return a + dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}

void rungekutta4bg(double &a, const double fourpiG, const cosmology cosmo, const double dtau)
{
    const double a_discontinuity = 0.341672;
    double a_new;

    // Step 1: Perform the full integration step
    a_new = rk4_step<OmegaLambdaNowDefault>(a, fourpiG, cosmo, dtau);

    // Step 2: Check if we crossed the discontinuity
    if (a < a_discontinuity && a_new > a_discontinuity)
    {
        // Step 3: Estimate dtau to reach the discontinuity
        double H_at_a = Hconf<OmegaLambdaNowDefault>(a, fourpiG, cosmo);
        double dtau_to_discontinuity = (a_discontinuity - a) / (a * H_at_a);
        double dtau_remaining = dtau - dtau_to_discontinuity;

        // Step 4: Integrate the first part using OmegaLambdaNegative
        a = rk4_step<OmegaLambdaNegative>(a, fourpiG, cosmo, dtau_to_discontinuity);

        // Step 5: Integrate the second part using OmegaLambdaPositive
        a = rk4_step<OmegaLambdaPositive>(a, fourpiG, cosmo, dtau_remaining);

        // Step 6: Check if a < a_discontinuity
        if (a < a_discontinuity)
        {
            // Handle the unexpected result
            fprintf(stderr, "Error: Integration resulted in a < a_discontinuity.\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        // No discontinuity crossed; update 'a'
        a = a_new;
    }
}


double particleHorizonIntegrand(double sqrta, void * cosmo)
{
	return 2. / (sqrta * Hconf(sqrta*sqrta, 1., *(cosmology *)cosmo));
}

//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: particle horizon (tau)
// 
//////////////////////////

double particleHorizon(const double a, const double fourpiG, cosmology & cosmo)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	
	f.function = &particleHorizonIntegrand;
	f.params = &cosmo;
	
	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);
	
	return result / sqrt(fourpiG);
}

#endif

