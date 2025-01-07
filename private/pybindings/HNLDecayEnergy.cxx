#include <icetray/load_project.h>
#include "LeptonInjector_HNL/HNLDecayEnergy.h"

using namespace boost::python;
namespace bp = boost::python;

void register_HNLDecayEnergy()
{
	// map the HNLDecayEnergy namespace to a sub-module  
	bp::object HNLDecayEnergyModule(bp::handle<>(bp::borrowed(PyImport_AddModule("LeptonInjector_HNL.HNLDecayEnergy"))));  
	// make "from icecube.LeptonInjector_HNL import HNLDecayEnergy" work  
	bp::scope().attr("HNLDecayEnergy") = HNLDecayEnergyModule;  
	// set the current scope to the new sub-module  
	bp::scope HNLDecayEnergy_scope = HNLDecayEnergyModule;  
    
	// export FullWidth in the HNLDecayEnergy namespace
	def("FullWidth", HNLDecayEnergy::FullWidth,
	  "Calculate total decay width of HNL particle.",
	  (bp::arg("m_HNL")));

	// export individual widths in the HNLDecayEnergy namespace
	def("gamma_nu_nu_nu_overload", HNLDecayEnergy::gamma_nu_nu_nu_overload,
	  "Calculate decay width of HNL to 3 neutrinos.",
	  (bp::arg("m_HNL")));

	def("gamma_e_e", HNLDecayEnergy::gamma_e_e,
	  "Calculate decay width of HNL to neutrino+electron+positron.",
	  (bp::arg("m_HNL")));

	def("gamma_p0_nu", HNLDecayEnergy::gamma_p0_nu,
	  "Calculate decay width of HNL to neutrino+pi0.",
	  (bp::arg("m_HNL")));

	def("gamma_nu_mu_mu", HNLDecayEnergy::gamma_nu_mu_mu,
	  "Calculate decay width of HNL to neutrino+muon+antimuon.",
	  (bp::arg("m_HNL")));

	def("gamma_eta_nu", HNLDecayEnergy::gamma_eta_nu,
	  "Calculate decay width of HNL to neutrino+eta.",
	  (bp::arg("m_HNL")));

	def("gamma_rho0_nu", HNLDecayEnergy::gamma_rho0_nu,
	  "Calculate decay width of HNL to neutrino+rho0.",
	  (bp::arg("m_HNL")));

	def("gamma_omega_nu", HNLDecayEnergy::gamma_omega_nu,
	  "Calculate decay width of HNL to neutrino+omega.",
	  (bp::arg("m_HNL")));

	def("gamma_etaprime_nu", HNLDecayEnergy::gamma_etaprime_nu,
	  "Calculate decay width of HNL to neutrino+etaprime.",
	  (bp::arg("m_HNL")));

	def("gamma_phi_nu", HNLDecayEnergy::gamma_phi_nu,
	  "Calculate decay width of HNL to neutrino+phi.",
	  (bp::arg("m_HNL")));

	def("gamma_nue_e_tau", HNLDecayEnergy::gamma_nue_e_tau,
	  "Calculate decay width of HNL to tau+nue+e.",
	  (bp::arg("m_HNL")));

	def("gamma_numu_mu_tau", HNLDecayEnergy::gamma_numu_mu_tau,
	  "Calculate decay width of HNL to tau+numu+mu.",
	  (bp::arg("m_HNL")));

	def("gamma_tau_pi", HNLDecayEnergy::gamma_tau_pi,
	  "Calculate decay width of HNL to tau+pi.",
	  (bp::arg("m_HNL")));

	def("gamma_tau_K", HNLDecayEnergy::gamma_tau_K,
	  "Calculate decay width of HNL to tau+K.",
	  (bp::arg("m_HNL")));

	def("gamma_tau_rho", HNLDecayEnergy::gamma_tau_rho,
	  "Calculate decay width of HNL to tau+rho.",
	  (bp::arg("m_HNL")));
}
