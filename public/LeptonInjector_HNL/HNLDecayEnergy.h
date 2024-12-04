#ifndef HNLDECAYENERGY_H_INCLUDED    // To make sure you don't declare the function more than once by including the header multiple times.
#define HNLDECAYENERGY_H_INCLUDED

#include <iostream>
#include "math.h"
#include "string.h"
#include <typeinfo> 
#include <vector>

namespace HNLDecayEnergy{

double TwoBodyRestE1(double m_HNL, double m1, double m2);
double GammaFromE(double E, double m);
double GammaFromBeta(double Beta);
double Beta(double E, double m);
std::vector<double> Lorentz(double E, std::vector<double> p, double boost[3], double gamma);
std::vector<double> RotateVector(std::vector<double> p, double alpha, int dir);

// Define all explicit decay modes, so they can be used as (external) functions
// nu_nu_nu
double gamma_nu_nu_nu_overload(double m_HNL);
// nu_e_e
double gamma_e_e(double m_HNL);
// pi0_nu
double gamma_p0_nu(double m_HNL);
// nu_mu_mu
double gamma_nu_mu_mu(double m_HNL);
// eta_nu
double gamma_eta_nu(double m_HNL);
// rho0_nu
double gamma_rho0_nu(double m_HNL);
// omega_nu
double gamma_omega_nu(double m_HNL);

// K*0_nu and K*0bar_nu are not accessible at tree level, apparently.

// etaprime_nu
double gamma_etaprime_nu(double m_HNL);
// phi_nu
double gamma_phi_nu(double m_HNL);
// nue_e_tau
double gamma_nue_e_tau(double m_HNL);
// numu_mu_tau
double gamma_numu_mu_tau(double m_HNL);
// tau_pi
double gamma_tau_pi(double m_HNL);
// // D*0_nu  FIXME: This is missing
// double gamma_D*0_nu(double m_HNL);
// // D*0bar_nu  FIXME: This is missing
// double gamma_D*0bar_nu(double m_HNL);
// tau_K
double gamma_tau_K(double m_HNL);
// tau_rho
double gamma_tau_rho(double m_HNL);
// // tau_K*+  FIXME: This is missing
// double gamma_tau_K*+(double m_HNL);
// // etac_nu  FIXME: This is missing
// double gamma_etac_nu(double m_HNL);

double FullWidth(double m_HNL);
int DecayChannel(double rand, double m_HNL);
void TwoBodyDecay(std::vector<double> &decay_products, bool isParticle, double pdgId_daughter, double m_daughter1, double m_daughter2, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double coszenith_daughter, double azimuth_daughter);
void ThreeBodyDecay(std::vector<double> &decay_products, bool isParticle, std::string file, double m_1, double m_2, double m_3, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, int rand_event, int channel);
// void ThreeBodyDecay_ThreeNeutrino(std::vector<double> &decay_products, bool isParticle, double m_1, double m_2, double m_3, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, int rand_event, int channel);
std::vector<double> GetDaughterVector(bool isParticle, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_channel, double rand_coszenith, double rand_azimuth, int rand_event);

};  // namespace HNLDecayEnergy

#endif
