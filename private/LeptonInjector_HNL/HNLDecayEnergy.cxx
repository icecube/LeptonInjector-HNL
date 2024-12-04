#include "LeptonInjector_HNL/DecayWidths.h"
#include "LeptonInjector_HNL/HNLDecayEnergy.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
//Borrowed include statements
#include <LeptonInjector_HNL/LeptonInjector_HNL.h>
#include <cassert>
#include <tuple>
#include <utility> 
#include <boost/math/constants/constants.hpp>
#include <boost/make_shared.hpp>

#include <icetray/open.h>
#include <icetray/I3Units.h>
#include <dataclasses/physics/I3MCTree.h>

static const double c = 1.0;

// Lepton masses in GeV
static const double m_electron = 511e-06;
static const double m_muon     = 0.10566;
static const double m_tau      = 1.777;
static const double m_neutrino = 0.0;
// Chosen to agree with values used in MadGraph

// Meson masses in GeV - taken from https://en.wikipedia.org/wiki/List_of_mesons, with the given digits (rounded)
static const double m_pi0       = 0.135;
static const double m_pi_plus   = 0.140;
static const double m_k_plus    = 0.494;
static const double m_rho       = 0.775;
static const double m_omega		= 0.783;
static const double m_eta       = 0.548;
static const double m_eta_prime = 0.958;
static const double m_phi		= 1.019;

// Decay constants in GeV (switched to the ones from https://arxiv.org/abs/2007.03701)
static const double f_pi        = 0.130;
static const double f_kaon      = 0.156;
static const double f_eta       = 0.0816;
static const double f_etaprime  = -0.0946;
static const double f_rho       = 0.171;
static const double f_omega		= 0.155;
static const double f_phi		= 0.232;

// CKM mixing elements (from https://doi.org/10.1093/ptep/ptaa104)
static const double V_ud = 0.97370;
static const double V_us = 0.2245;

// Declare global physical constants
static const double G_F         = 1.17E-05;  // Fermi constant in GeV^-2
static const double sin2thetaW  = 0.2229;  // sin^2(weak_angle) from CODATA 2018

static int counter = 0;

double HNLDecayEnergy::TwoBodyRestE1(double m_HNL, double m1, double m2){
    return (pow(m_HNL,2)+pow(m1,2)-pow(m2,2))/(2*m_HNL);
}

double HNLDecayEnergy::GammaFromE(double E, double m){
    return E/m;
}

double HNLDecayEnergy::GammaFromBeta(double Beta){
    return 1./sqrt(1-Beta*Beta);
}

double HNLDecayEnergy::Beta(double E, double m){
    return sqrt(1-pow(m/E,2));
}

std::vector<double> HNLDecayEnergy::Lorentz(double E, std::vector<double> p, double boost[3], double gamma){

	std::vector<double> p_prime;
    
	double boost_normsq = pow(boost[0],2)+pow(boost[1],2)+pow(boost[2],2);
    
    if (sqrt(boost_normsq)>1){
        log_fatal_stream("Relativity Error: passed boost with velocity > c , boost magnitude = " << sqrt(boost_normsq) );
    }
    
	// active transformation 
    double lorentz_matrix[4][4] = {
        { gamma, gamma*boost[0]/c, gamma*boost[1]/c, gamma*boost[2]/c },
        { gamma*boost[0]/c, 1 + (gamma - 1)*boost[0]*boost[0]/boost_normsq, (gamma - 1)*boost[0]*boost[1]/boost_normsq, (gamma - 1)*boost[0]*boost[2]/boost_normsq },
        { gamma*boost[1]/c, (gamma - 1)*boost[0]*boost[1]/boost_normsq, 1 + (gamma - 1)*boost[1]*boost[1]/boost_normsq, (gamma - 1)*boost[1]*boost[2]/boost_normsq },
        { gamma*boost[2]/c, (gamma - 1)*boost[0]*boost[2]/boost_normsq, (gamma - 1)*boost[1]*boost[2]/boost_normsq, 1 + (gamma - 1)*boost[2]*boost[2]/boost_normsq }
    };

    p_prime.push_back(E*lorentz_matrix[1][0] + p[0]*lorentz_matrix[1][1] + p[1]*lorentz_matrix[1][2] + p[2]*lorentz_matrix[1][3]);
    p_prime.push_back(E*lorentz_matrix[2][0] + p[0]*lorentz_matrix[2][1] + p[1]*lorentz_matrix[2][2] + p[2]*lorentz_matrix[2][3]);
    p_prime.push_back(E*lorentz_matrix[3][0] + p[0]*lorentz_matrix[3][1] + p[1]*lorentz_matrix[3][2] + p[2]*lorentz_matrix[3][3]);
    double E_prime =  E*lorentz_matrix[0][0] + p[0]*lorentz_matrix[0][1] + p[1]*lorentz_matrix[0][2] + p[2]*lorentz_matrix[0][3];
    p_prime.push_back(E_prime);
    
    /********
    Be aware!
    Components here are p_prime = (p'_x, p'_y, p'_z, E')
    IN THAT ORDER
    *********/

    double m2 = E*E - p[0]*p[0] - p[1]*p[1] - p[2]*p[2];
    double m_prime2 = p_prime[3]*p_prime[3] - p_prime[0]*p_prime[0] - p_prime[1]*p_prime[1] - p_prime[2]*p_prime[2];
    return p_prime;
}

std::vector<double> HNLDecayEnergy::RotateVector(std::vector<double> p, double alpha, int dir){
	std::vector<double> p_rot;

	if(dir==0){
		p_rot.push_back(p.at(0));
		p_rot.push_back(cos(alpha)*p.at(1) - sin(alpha)*p.at(2));
		p_rot.push_back(sin(alpha)*p.at(1) + cos(alpha)*p.at(2));
	}
	else if(dir==2){
		p_rot.push_back(cos(alpha)*p.at(0) - sin(alpha)*p.at(1));
		p_rot.push_back(sin(alpha)*p.at(0) + cos(alpha)*p.at(1));
		p_rot.push_back(p.at(2));
	}

	return p_rot;
}

// Define all explicit decay modes, so they can be used as (external) functions (through pybingings)

// nu_nu_nu
double HNLDecayEnergy::gamma_nu_nu_nu_overload(double m_HNL){
	return gamma_nu_nu_nu(G_F, m_HNL);
	}

// nu_e_e
double HNLDecayEnergy::gamma_e_e(double m_HNL){
	return gamma_nu_ll(G_F, sin2thetaW, m_electron, m_HNL);
	}

// pi0_nu
double HNLDecayEnergy::gamma_p0_nu(double m_HNL){
	return gamma_meson_nu(G_F, f_pi, m_pi0, m_HNL);
	}

// nu_mu_mu
double HNLDecayEnergy::gamma_nu_mu_mu(double m_HNL){
	return gamma_nu_ll(G_F, sin2thetaW, m_muon, m_HNL);
	}

// eta_nu
double HNLDecayEnergy::gamma_eta_nu(double m_HNL){
	return gamma_meson_nu(G_F, f_eta, m_eta, m_HNL);
	}

// rho0_nu
double HNLDecayEnergy::gamma_rho0_nu(double m_HNL){
	double g_v = 1-2*sin2thetaW;
	return gamma_V_nu(f_rho, g_v, G_F, m_rho, m_HNL);
	}

// omega_nu
double HNLDecayEnergy::gamma_omega_nu(double m_HNL){
	double g_v = (-2./3.)*sin2thetaW;
	return gamma_V_nu(f_omega, g_v, G_F, m_omega, m_HNL);
	}

// K*0_nu and K*0bar_nu are not accessible at tree level, apparently.

// etaprime_nu
double HNLDecayEnergy::gamma_etaprime_nu(double m_HNL){
	return gamma_meson_nu(G_F, f_etaprime, m_eta_prime, m_HNL);
	}

// phi_nu
double HNLDecayEnergy::gamma_phi_nu(double m_HNL){
	double g_v = (-1.)*sqrt(2)*(0.5-(2./3.)*sin2thetaW);
	return gamma_V_nu(f_phi, g_v, G_F, m_phi, m_HNL);
	}

// nue_e_tau
double HNLDecayEnergy::gamma_nue_e_tau(double m_HNL){
	return gamma_nu1_l1_l2(G_F, m_electron, m_tau, m_HNL);
	}

// numu_mu_tau
double HNLDecayEnergy::gamma_numu_mu_tau(double m_HNL){
	return gamma_nu1_l1_l2(G_F, m_muon, m_tau, m_HNL);
	}

// tau_pi
double HNLDecayEnergy::gamma_tau_pi(double m_HNL){
	return gamma_H_l(G_F, V_ud, f_pi, m_pi_plus, m_tau, m_HNL);
	}

// tau_K
double HNLDecayEnergy::gamma_tau_K(double m_HNL){
	return gamma_H_l(G_F, V_us, f_kaon, m_k_plus, m_tau, m_HNL);
	}

// tau_rho
double HNLDecayEnergy::gamma_tau_rho(double m_HNL){
	return gamma_rho_l(f_rho, G_F, V_ud, m_rho, m_tau, m_HNL);
	}


double HNLDecayEnergy::FullWidth(double m_HNL){
    return (
		gamma_nu_nu_nu_overload(m_HNL) +		// nu_nu_nu
		gamma_e_e(m_HNL) +            		 	// nu_e_e
		gamma_p0_nu(m_HNL) +                    // pi0_nu
		gamma_nu_mu_mu(m_HNL) +                 // nu_mu_mu
		gamma_eta_nu(m_HNL) +                   // eta_nu
		gamma_rho0_nu(m_HNL) +         			// rho0_nu
		gamma_omega_nu(m_HNL) +					// omega_nu
		gamma_etaprime_nu(m_HNL) +       		// etaprime_nu
		gamma_phi_nu(m_HNL) +					// phi_nu
		gamma_nue_e_tau(m_HNL) +                // nue_e_tau
		gamma_numu_mu_tau(m_HNL) +              // numu_mu_tau
		gamma_tau_pi(m_HNL) +      		     	// tau_pi
		gamma_tau_K(m_HNL) +       			 	// tau_K
		gamma_tau_rho(m_HNL)       		     	// tau_rho
	);
}

int HNLDecayEnergy::DecayChannel(double rand, double m_HNL){
	double fullwidth = FullWidth(m_HNL);

//     double x[13];
// 	x[0]   = gamma_nu_nu_nu_overload(m_HNL)/fullwidth;
// 	x[1]   = x[0]+gamma_e_e(m_HNL)/fullwidth;
// 	x[2]   = x[1]+gamma_p0_nu(m_HNL)/fullwidth;
// 	x[3]   = x[2]+gamma_nu_mu_mu(m_HNL)/fullwidth;
// 	x[4]   = x[3]+gamma_eta_nu(m_HNL)/fullwidth;
// 	x[5]   = x[4]+gamma_rho0_nu(m_HNL)/fullwidth;
// 	x[6]   = x[5]+gamma_omega_nu(m_HNL)/fullwidth;
// 	x[7]   = x[6]+gamma_etaprime_nu(m_HNL)/fullwidth;
// 	x[8]   = x[7]+gamma_phi_nu(m_HNL)/fullwidth;
// 	x[9]   = x[8]+gamma_nue_e_tau(m_HNL)/fullwidth;
// 	x[10]  = x[9]+gamma_numu_mu_tau(m_HNL)/fullwidth;
// 	x[11]  = x[10]+gamma_tau_pi(m_HNL)/fullwidth;
// 	x[12]  = x[11]+gamma_tau_K(m_HNL)/fullwidth;
// 	//double x13 = gamma_tau_rho(m_HNL)/fullwidth;  // captured by last return statement
    
    std::vector<double> x;
    x.push_back(gamma_nu_nu_nu_overload(m_HNL)/fullwidth);
	x.push_back(x[0]+gamma_e_e(m_HNL)/fullwidth);
	x.push_back(x[1]+gamma_p0_nu(m_HNL)/fullwidth);
	x.push_back(x[2]+gamma_nu_mu_mu(m_HNL)/fullwidth);
	x.push_back(x[3]+gamma_eta_nu(m_HNL)/fullwidth);
	x.push_back(x[4]+gamma_rho0_nu(m_HNL)/fullwidth);
	x.push_back(x[5]+gamma_omega_nu(m_HNL)/fullwidth);
	x.push_back(x[6]+gamma_etaprime_nu(m_HNL)/fullwidth);
	x.push_back(x[7]+gamma_phi_nu(m_HNL)/fullwidth);
	x.push_back(x[8]+gamma_nue_e_tau(m_HNL)/fullwidth);
	x.push_back(x[9]+gamma_numu_mu_tau(m_HNL)/fullwidth);
	x.push_back(x[10]+gamma_tau_pi(m_HNL)/fullwidth);
	x.push_back(x[11]+gamma_tau_K(m_HNL)/fullwidth);

// 	double x[13] = {
//                     x0,
//                     x0+x1,
//                     x0+x1+x2,
//                     x0+x1+x2+x3,
//                     x0+x1+x2+x3+x4,
//                     x0+x1+x2+x3+x4+x5,
//                     x0+x1+x2+x3+x4+x5+x6,
//                     x0+x1+x2+x3+x4+x5+x6+x7,
//                     x0+x1+x2+x3+x4+x5+x6+x7+x8,
//                     x0+x1+x2+x3+x4+x5+x6+x7+x8+x9,
//                     x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
//                     x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11,
//                     x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12,
//                     };///Added/removed print statement of x bins here
    
//     std::cout<<"Bin edges: "<<std::endl;
//     for (int i=0; i<13; i++){
//         std::cout<<x[i]<< ", ";
//     }
//     std::cout<<std::endl;
//     std::cout<<"rand: "<< rand <<std::endl;
    
//         std::cout << "Some generic print statement" << std::endl;
	for(int i=0; i<13; i++){
		if(rand<=x[i]) return i+1;
//         if(rand == 1) return i+1;
	}
	return 14;
}

// Pushes back the following elements to the vector: the pdgid, the energy and the direction of the daughter produced in the HNL 2-body decay 
void HNLDecayEnergy::TwoBodyDecay(std::vector<double> &decay_products, bool isParticle, double pdgId_daughter, double m_daughter1, double m_daughter2, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double coszenith_daughter, double azimuth_daughter){
    // Compute boost
    double speed_HNL = c*Beta(E_HNL, m_HNL);
    double v_HNL[3] = {speed_HNL*sin(zenith_HNL)*cos(azimuth_HNL),speed_HNL*sin(zenith_HNL)*sin(azimuth_HNL),speed_HNL*cos(zenith_HNL)};
    //Incorrect trig in velocity calculation? Nope, was correct.

    // Compute rest frame variables
	double E_daughter = TwoBodyRestE1(m_HNL, m_daughter1, m_daughter2);
	double beta_daughter = Beta(E_daughter, m_daughter1);
    
//     double p_daughter_x = beta_daughter*c*m_daughter1*sin(acos(coszenith_daughter))*cos(azimuth_daughter);
//     double p_daughter_y = beta_daughter*c*m_daughter1*sin(acos(coszenith_daughter))*sin(azimuth_daughter);
//     double p_daughter_z = beta_daughter*c*m_daughter1*coszenith_daughter;
//     //Incorrect trig? Nope, was correct.
    
    double p_magnitude = (pow(m_HNL, 2) - pow(m_daughter1, 2))/(2*m_HNL);
    double p_daughter_x = p_magnitude*sin(acos(coszenith_daughter))*cos(azimuth_daughter);
    double p_daughter_y = p_magnitude*sin(acos(coszenith_daughter))*sin(azimuth_daughter);
    double p_daughter_z = p_magnitude*coszenith_daughter;

    
	std::vector<double> p_daughter{p_daughter_x,p_daughter_y,p_daughter_z};

    // Lorentz-transform the momentum components
	std::vector<double> p_daughter_prime;
	p_daughter_prime = Lorentz(E_daughter,p_daughter,v_HNL,GammaFromE(E_HNL,m_HNL));

	// Compute daughter energy and angles
	double E_daughter_prime = p_daughter_prime.at(3);
	double zenith_daughter_prime = acos(p_daughter_prime.at(2)/sqrt(pow(p_daughter_prime.at(0),2)+pow(p_daughter_prime.at(1),2)+pow(p_daughter_prime.at(2),2)));
    //fix zenith so cos(zenith)=1 is upgoing?
    //zenith_daughter_prime += M_PI;
	double azimuth_daughter_prime = atan(p_daughter_prime.at(1)/p_daughter_prime.at(0));
	if(p_daughter_prime.at(0)>0 && p_daughter_prime.at(1)<0) azimuth_daughter_prime = 2*M_PI + azimuth_daughter_prime;
	else if(p_daughter_prime.at(0)<0 && p_daughter_prime.at(1)>0) azimuth_daughter_prime = M_PI + azimuth_daughter_prime;
	else if(p_daughter_prime.at(0)<0 && p_daughter_prime.at(1)<0) azimuth_daughter_prime = M_PI + azimuth_daughter_prime;
    
    double E_neutrino = m_HNL-E_daughter;

	double pdgId_neutrino = 16;
	if(!isParticle) pdgId_neutrino = -16;
    std::vector<double> p_neutrino{-p_daughter_x, -p_daughter_y, -p_daughter_z};
    std::vector<double> p_neutrino_prime;
	p_neutrino_prime = Lorentz(E_neutrino, p_neutrino, v_HNL, GammaFromE(E_HNL,m_HNL));
    
    // Compute neutrino energy and angles
	double E_neutrino_prime = p_neutrino_prime.at(3);
	double zenith_neutrino_prime = acos(p_neutrino_prime.at(2)/sqrt(pow(p_neutrino_prime.at(0),2) + pow(p_neutrino_prime.at(1),2) + pow(p_neutrino_prime.at(2),2)));
    //fix zenith so cos(zenith)=1 is upgoing?
    //zenith_neutrino_prime += M_PI;
	double azimuth_neutrino_prime = atan(p_neutrino_prime.at(1)/p_neutrino_prime.at(0));
	if(p_neutrino_prime.at(0)>0 && p_neutrino_prime.at(1)<0) azimuth_neutrino_prime = 2*M_PI + azimuth_neutrino_prime;
	else if(p_neutrino_prime.at(0)<0 && p_neutrino_prime.at(1)>0) azimuth_neutrino_prime = M_PI + azimuth_neutrino_prime;
	else if(p_neutrino_prime.at(0)<0 && p_neutrino_prime.at(1)<0) azimuth_neutrino_prime = M_PI + azimuth_neutrino_prime;
    
    // Fill vector
	decay_products.push_back(pdgId_daughter);
	decay_products.push_back(E_daughter_prime);
	decay_products.push_back(zenith_daughter_prime);
	decay_products.push_back(azimuth_daughter_prime);
    
    decay_products.push_back(pdgId_neutrino);  // PDG ID for tau neutrino/antineutrino
	decay_products.push_back(E_neutrino_prime);
	decay_products.push_back(zenith_neutrino_prime);
	decay_products.push_back(azimuth_neutrino_prime);

    if(E_daughter_prime <0 || E_neutrino_prime <0){
        std::cout<<" Energies: " << E_daughter_prime << ", " << E_neutrino_prime <<std::endl;
    }
    
    //////// Write out extra information for debugging purposes
    
    std::ofstream outputFile; // Declare an output file stream object
    // Open the file in append mode
    outputFile.open("/data/user/jbook/I3_HNL_Search/validation/new_daughter_energies1.txt", std::ios::app); 
    outputFile << E_daughter_prime << ", " << E_neutrino_prime << ", "<< E_daughter << ", " << E_neutrino<< ", " << E_HNL <<std::endl; // Write the new data to the file
    outputFile.close(); // Close the file
    
}

// Pushes back the following elements to the vector: the pdgid, the energy and the direction of the two visible daughters produced in the HNL 3-body decay 
void HNLDecayEnergy::ThreeBodyDecay(std::vector<double> &decay_products, bool isParticle, std::string file, double m_1, double m_2, double m_3, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, int rand_event, int channel){

    // Compute boost
    double v_HNL[3] = {c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*cos(azimuth_HNL),c*Beta(E_HNL,m_HNL)*sin(zenith_HNL)*sin(azimuth_HNL),c*Beta(E_HNL,m_HNL)*cos(zenith_HNL)};
    
    double norm_sqr_v = v_HNL[0]*v_HNL[0] + v_HNL[1]*v_HNL[1] + v_HNL[2]*v_HNL[2];
    // Open file and get variables
	int pdgId_1, pdgId_2, pdgId_3;
	double E_1, E_2, E_3;
    double px_1, px_2, px_3;
    double py_1, py_2, py_3;
    double pz_1, pz_2, pz_3;
    
	std::ifstream in(file);
	std::string line;
	std::vector<std::string> splitted_line;
	int i = 1;
    
	while(std::getline(in,line)){
        if(i==rand_event){
			std::stringstream ss(line);
			std::vector<std::string> splitted_line(std::istream_iterator<std::string>{ss},std::istream_iterator<std::string>());
			
            pdgId_1 = std::stoi(splitted_line[0]);
			pdgId_2 = std::stoi(splitted_line[4]);
            pdgId_3 = std::stoi(splitted_line[8]);
            
            px_1 = std::stod(splitted_line[1]);
            px_2 = std::stod(splitted_line[5]);
            px_3 = std::stod(splitted_line[9]);
            
            double px_tot = px_1 + px_2 + px_3;
            
            py_1 = std::stod(splitted_line[2]);
			py_2 = std::stod(splitted_line[6]);
            py_3 = std::stod(splitted_line[10]);
            
            double py_tot = py_1 + py_2 + py_3;
            
			pz_1 = std::stod(splitted_line[3]);
			pz_2 = std::stod(splitted_line[7]);
            pz_3 = std::stod(splitted_line[11]);
            
            double pz_tot = pz_1 + pz_2 + pz_3;

            if(sqrt(px_1*px_1 + py_1*py_1 + pz_1*pz_1) + sqrt(px_2*px_2 + py_2*py_2 + pz_2*pz_2) + sqrt(px_3*px_3 + py_3*py_3 + pz_3*pz_3) <= m_HNL + 1.e-4){
                break;
            }
            std::cout << "****** Rejected ********" << std::endl;
		}
		i++;
	} 

	std::vector<double> p_1 = {px_1, py_1, pz_1};
	std::vector<double> p_2 = {px_2, py_2, pz_2};
    std::vector<double> p_3 = {px_3, py_3, pz_3};

    E_1 = c * sqrt( (px_1*px_1 + py_1*py_1 + pz_1*pz_1) + pow(m_1*c, 2));
    E_2 = c * sqrt( (px_2*px_2 + py_2*py_2 + pz_2*pz_2) + pow(m_2*c, 2));
    E_3 = c * sqrt( (px_3*px_3 + py_3*py_3 + pz_3*pz_3) + pow(m_3*c, 2));

    // Lorentz-transform the momentum components
	std::vector<double> p_1_prime, p_2_prime, p_3_prime;
	p_1_prime = Lorentz(E_1,p_1,v_HNL,GammaFromE(E_HNL,m_HNL));
	p_2_prime = Lorentz(E_2,p_2,v_HNL,GammaFromE(E_HNL,m_HNL));
    p_3_prime = Lorentz(E_3,p_3,v_HNL,GammaFromE(E_HNL,m_HNL));
	
    // Compute daughter_1 energy and angles
	double E_1_prime = p_1_prime.at(3);
	double zenith_1_prime = acos(p_1_prime.at(2)/sqrt(pow(p_1_prime.at(0),2)+pow(p_1_prime.at(1),2)+pow(p_1_prime.at(2),2)));
	double azimuth_1_prime = atan(p_1_prime.at(1)/p_1_prime.at(0));
    
	if(p_1_prime.at(0)>0 && p_1_prime.at(1)<0){
        azimuth_1_prime = 2*M_PI + azimuth_1_prime;
    }
	else if(p_1_prime.at(0)<0 && p_1_prime.at(1)>0){
        azimuth_1_prime = M_PI + azimuth_1_prime;
    }
	else if(p_1_prime.at(0)<0 && p_1_prime.at(1)<0){
        azimuth_1_prime = M_PI + azimuth_1_prime;
    }
	
    // Compute daughter_2 energy and angles
	double E_2_prime = p_2_prime.at(3);
	double zenith_2_prime = acos(p_2_prime.at(2)/sqrt(pow(p_2_prime.at(0),2)+pow(p_2_prime.at(1),2)+pow(p_2_prime.at(2),2)));
	double azimuth_2_prime = atan(p_2_prime.at(1)/p_2_prime.at(0));
	
    if(p_2_prime.at(0)>0 && p_2_prime.at(1)<0){
        azimuth_2_prime = 2*M_PI + azimuth_2_prime;
    }
	else if(p_2_prime.at(0)<0 && p_2_prime.at(1)>0){
        azimuth_2_prime = M_PI + azimuth_2_prime;
    }
	else if(p_2_prime.at(0)<0 && p_2_prime.at(1)<0){
        azimuth_2_prime = M_PI + azimuth_2_prime;
    }
    
    // Compute daughter_3 energy and angles
	double E_3_prime = p_3_prime.at(3);
	double zenith_3_prime = acos(p_3_prime.at(2)/sqrt(pow(p_3_prime.at(0),2)+pow(p_3_prime.at(1),2)+pow(p_3_prime.at(2),2)));
	double azimuth_3_prime = atan(p_3_prime.at(1)/p_3_prime.at(0));
	
    if(p_3_prime.at(0)>0 && p_3_prime.at(1)<0){
        azimuth_3_prime = 2*M_PI + azimuth_3_prime;
    }
	else if(p_3_prime.at(0)<0 && p_3_prime.at(1)>0){
        azimuth_3_prime = M_PI + azimuth_3_prime;
    }
	else if(p_3_prime.at(0)<0 && p_3_prime.at(1)<0){
        azimuth_3_prime = M_PI + azimuth_3_prime;
    }

	if(!isParticle){
		pdgId_1 *= -1;
		pdgId_2 *= -1;
		pdgId_3 *= -1;
	}
    
    /////////////////////////////

    // Fill vector
	decay_products.push_back(pdgId_1);
	decay_products.push_back(E_1_prime);
	decay_products.push_back(zenith_1_prime);
	decay_products.push_back(azimuth_1_prime);
    
	decay_products.push_back(pdgId_2);
	decay_products.push_back(E_2_prime);
	decay_products.push_back(zenith_2_prime);
	decay_products.push_back(azimuth_2_prime);
    
	decay_products.push_back(pdgId_3);
	decay_products.push_back(E_3_prime);
	decay_products.push_back(zenith_3_prime);
	decay_products.push_back(azimuth_3_prime);
}

std::vector<double> HNLDecayEnergy::GetDaughterVector(bool isParticle, double m_HNL, double E_HNL, double zenith_HNL, double azimuth_HNL, double rand_channel, double rand_coszenith, double rand_azimuth, int rand_event){

    bool majorana = false; 
	
    // Declare vector
	std::vector<double> decay_products;

	// Determine the decay channel
	int channel = DecayChannel(rand_channel,m_HNL);
//     std::cout<<"rand channel = " <<rand_channel<<std::endl;
//     std::cout<<"channel = " <<channel<<std::endl;
    
    std::string i3_src = std::getenv("I3_SRC");
	std::string file_prefix;
	std::string mass_string;
	
	if (majorana){
   		file_prefix = i3_src+"/LeptonInjector/resources/madgraph/DatosMajorana/";
    	if(m_HNL == 0.1){
            mass_string = "HNL_Decay_Majorana_0.1GeV";
        }
        else if (m_HNL == 0.3){
            mass_string = "HNL_Decay_Majorana_0.3GeV";
        }
        else if(m_HNL == 0.6){
            mass_string = "HNL_Decay_Majorana_0.6GeV";
        }
        else if(m_HNL == 1.0){
            mass_string = "HNL_Decay_Majorana_1GeV";
        }
        else{
            std::cout << "HNL mass = " << m_HNL << std::endl;
            std::cout << "Unsupported mass choice - three body decays will break" << std::endl;
        }
    }
    else{
       file_prefix = i3_src+"/LeptonInjector/resources/madgraph/DatosDirac/";
        if(m_HNL == 0.1){
            mass_string = "HNL_Decay_Dirac_0.1GeV";
        }
        else if (m_HNL == 0.3){
            mass_string = "HNL_Decay_Dirac_0.3GeV";
        }
        else if(m_HNL == 0.6){
            mass_string = "HNL_Decay_Dirac_0.6GeV";
        }
        else if(m_HNL == 1.0){
            mass_string = "HNL_Decay_Dirac_1GeV";
        }
        else{
            std::cout << "HNL mass = " << m_HNL << std::endl;
            std::cout << "Unsupported mass choice - three body decays will break" << std::endl;
        }
        
    }

	if(channel==1){
        std::string madgraph_file = file_prefix+mass_string+"_vtveve.dat";
        ThreeBodyDecay(decay_products, isParticle, madgraph_file, m_neutrino, m_neutrino, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_event, channel);
    }
	if(channel==2){
        std::string madgraph_file = file_prefix+mass_string+"_vee.dat";
        ThreeBodyDecay(decay_products, isParticle, madgraph_file, m_neutrino, m_electron, m_electron, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_event, channel);
    }
	else if(channel==3) TwoBodyDecay(decay_products, isParticle, 11, m_pi0, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_coszenith, rand_azimuth);
	else if(channel==4){
		std::string madgraph_file = file_prefix+mass_string+"_vmm.dat";
		ThreeBodyDecay(decay_products, isParticle, madgraph_file, m_neutrino, m_muon, m_muon, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_event, channel);
	}
	else if(channel==5) TwoBodyDecay(decay_products, isParticle, 221, m_eta, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_coszenith, rand_azimuth);
	else if(channel==6) TwoBodyDecay(decay_products, isParticle, 113, m_rho, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_coszenith, rand_azimuth);
	else if(channel==7) TwoBodyDecay(decay_products, isParticle, 223, m_omega, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_coszenith, rand_azimuth);
	else if(channel==8) TwoBodyDecay(decay_products, isParticle, 331, m_eta_prime, m_neutrino, m_HNL, E_HNL, zenith_HNL, azimuth_HNL, rand_coszenith, rand_azimuth);

    // Save decay channel
    decay_products.push_back(channel);
    
    //Return the vector
	return decay_products;
}
