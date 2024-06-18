#include "LeptonInjector_HNL/DecayWidths.h"

/*
Notes:
-All decay widths are without the leading |U_alpha|^2, so they can be used for several mixing values (by re-weighting during the analysis).
-The paper says that formulae for decay rates are identical for dirac and majorana, but for majorana the total width is 2x the sum of all, to account for charge-conjugated decay modes, how do we apply this, what type do we use? (Compare to what happens in HNLDecayEnergy.cxx)
*/

// Two-body decays:
// match results from https://arxiv.org/abs/2007.03701v4 - LF

// neutrino + meson (pi0/eta/eta_prime)
double gamma_meson_nu(double G_F, double f, double m_meson, double m_HNL)
{
	if (m_HNL < m_meson)
		return 0;
	double mass_ratio = pow(m_meson / m_HNL, 2);
	double factor = pow(G_F, 2) * pow(f, 2) * pow(m_HNL, 3) / (32 * M_PI);
	return factor * pow(1 - mass_ratio, 2);
}

// H (+) + lepton (-)
double gamma_H_l(double G_F, double V, double f, double m_H, double m_l, double m_HNL)
{
	if (m_HNL < m_H + m_l)
		return 0;
	double mass_ratio1 = pow(m_l / m_HNL, 2);
	double mass_ratio2 = pow(m_H / m_HNL, 2);
	double mass_diff_minus = pow((m_H - m_l) / m_HNL, 2);
	double mass_diff_plus = pow((m_H + m_l) / m_HNL, 2);
	double factor = pow(G_F, 2) * pow(V, 2) * pow(f, 2) * pow(m_HNL, 3) / (16 * M_PI);
	return factor * (pow(1 - mass_ratio1, 2) - mass_ratio2 * (1 + mass_ratio1)) * sqrt((1 - mass_diff_minus) * (1 - mass_diff_minus));
}

// rho (+) + lepton (-)
double gamma_rho_l(double f_v, double G_F, double V, double m_rho, double m_l, double m_HNL)
{
	if (m_HNL < m_rho + m_l)
		return 0;
	double mass_ratio1 = pow(m_l / m_HNL, 2);
	double mass_ratio2 = pow(m_rho / m_HNL, 2);
	double mass_diff_minus = pow((m_rho - m_l) / m_HNL, 2);
	double mass_diff_plus = pow((m_rho + m_l) / m_HNL, 2);
	double mass_diff_minus_2 = (pow(m_l, 2) - 2 * pow(m_rho, 2)) / pow(m_HNL, 2);
	double factor = pow(f_v / m_rho, 2) * pow(G_F, 2) * pow(V, 2) * pow(m_HNL, 3) / (8 * M_PI);
	return factor * (pow(1 - mass_ratio1, 2) + mass_ratio2 * (1 + mass_diff_minus_2)) * sqrt((1 - mass_diff_minus) * (1 - mass_diff_minus)); // found bug here!
}

// (neutral) vector meson + neutrino
double gamma_V_nu(double f_v, double g_v, double G_F, double m_V, double m_HNL)
{
	if (m_HNL < m_V)
		return 0;
	double mass_ratio = pow(m_V / m_HNL, 2);
	double factor = pow(G_F, 2) * pow(m_HNL, 3) * pow(f_v / m_V, 2) * pow(g_v, 2) / (32 * M_PI);
	return factor * (1 + 2 * mass_ratio) * pow(1 - mass_ratio, 2);
}

// Three-body decays:
// match results from https://arxiv.org/abs/2007.03701v4 - LF

// neutrino a + neutrino b (+) + neutrino b (-)
double gamma_nu_nu_nu(double G_F, double m_HNL)
{
	double factor = pow(G_F, 2) * pow(m_HNL, 5) / (192 * pow(M_PI, 3)); // common factor shared by all 3-body decay widths
	return factor;
}

//  lepton a (-) + neutrino b + lepton b (+)
double gamma_nu1_l1_l2(double G_F, double m_l1, double m_l2, double m_HNL)
{
	if (m_HNL < m_l1 + m_l2)
		return 0;
	double factor = pow(G_F, 2) * pow(m_HNL, 5) / (192 * pow(M_PI, 3)); // common factor shared by all 3-body decay widths
	double mass_ratio = std::max(m_l1, m_l2) / m_HNL;
	double x_factor = 1 - (8 * pow(mass_ratio, 2)) + (8 * pow(mass_ratio, 6)) - pow(mass_ratio, 8) - (12 * pow(mass_ratio, 4) * log(pow(mass_ratio, 2)));
	return factor * x_factor;
}

// neutrino a + lepton b (+) + lepton b (-)
// gamma_same_flavor, C3 and C4 are not used (only need a =! b case, because nutau+tau(+)+tau(-) is not accessible below 3 GeV masses) - commented out
// added return 0 statements, to make sure the decaywidth is really 0 if squareroot/fraction/log are ill defined - might this be the bug that breaks the mu mu decay channel kinematics?
double gamma_nu_ll(double G_F, double sin2thetaW, double m_l, double m_HNL)
{
	if (m_HNL < 2 * m_l)
		return 0;
	double x = pow(m_l / m_HNL, 2);
	double L, num, den;
	if (x <= 0.25)
	{ // otherwise denominator breaks
		num = 1 - 3 * x - (1 - x) * sqrt(1 - 4 * x);
		den = x * (1 + sqrt(1 - 4 * x));
	}
	if (den != 0 && num / den > 0)
		L = log(num / den); // otherwise fraction or log breaks
	double C1 = 0.25 * (1 - 4 * sin2thetaW + 8 * pow(sin2thetaW, 2));
	double C2 = 0.5 * sin2thetaW * (2 * sin2thetaW - 1);
	// double C3 = 0.25*(1+4*sin2thetaW+8*pow(sin2thetaW,2));
	// double C4 = 0.5*sin2thetaW*(2*sin2thetaW+1);
	double factor = pow(G_F, 2) * pow(m_HNL, 5) / (192 * pow(M_PI, 3)); // common factor shared by all 3-body decay widths
	// double gamma_same_flavor = (C3*((1-14*x-2*pow(x,2)-12*pow(x,3))*sqrt(1-4*x)+12*pow(x,2)*(pow(x,2)-1)*L)+(4*C4*(x*(2+10*x-12*pow(x,2))*sqrt(1-4*x)+6*pow(x,2)*(1-2*x+2*pow(x,2))*L)));
	double gamma_opposite_flavor = (C1 * ((1 - 14 * x - 2 * pow(x, 2) - 12 * pow(x, 3)) * sqrt(1 - 4 * x) + 12 * pow(x, 2) * (pow(x, 2) - 1) * L) + (4 * C2 * (x * (2 + 10 * x - 12 * pow(x, 2)) * sqrt(1 - 4 * x) + 6 * pow(x, 2) * (1 - 2 * x + 2 * pow(x, 2)) * L)));
	return factor * gamma_opposite_flavor;
}
