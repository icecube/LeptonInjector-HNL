#ifndef DECAYWIDTHS_H_INCLUDED  // To make sure you don't declare the function more than once by including the header multiple times.
#define DECAYWIDTHS_H_INCLUDED

#include <iostream>
#include <math.h>

double gamma_meson_nu(double G_F, double f, double m_meson, double m_HNL);
double gamma_H_l(double G_F, double V, double f, double m_H, double m_l, double m_HNL);
double gamma_rho_l(double g_p, double G_F, double V, double m_rho, double m_l, double m_HNL);
// double gamma_rho0_nu(double f_rho, double G_F, double m_rho0, double m_HNL);
double gamma_V_nu(double f_v, double g_v, double G_F, double m_V, double m_HNL);
double gamma_nu_nu_nu(double G_F, double m_HNL);
double gamma_nu1_l1_l2(double G_F, double m_l1, double m_l2, double m_HNL);
double gamma_nu_ll(double G_F, double sin2thetaW, double m_l, double m_HNL);

#endif
