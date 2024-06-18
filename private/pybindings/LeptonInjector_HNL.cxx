#include <icetray/load_project.h>
#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/stream_to_string.hpp>
#include <dataclasses/ostream_overloads.hpp>
#include <tableio/converter/pybindings.h>
#include <LeptonInjector_HNL/HNLDecayEnergy.h>
#include <LeptonInjector_HNL/LeptonInjector.h>
#include <LeptonInjector_HNL/converter/LeptonInjectionConfigurationConverter.h>

//void register_HNLDecayEnergy();
using namespace boost::python;
namespace bp = boost::python;

void register_HNLDecayEnergy()
{
	// map the HNLDecayEnergy namespace to a sub-module  
	bp::object HNLDecayEnergyModule(bp::handle<>(bp::borrowed(PyImport_AddModule("LeptonInjector_HNL.HNLDecayEnergy"))));  
	// make "from icecube.LeptonInjector import HNLDecayEnergy" work  
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

using namespace boost::python;

I3_PYTHON_MODULE(LeptonInjector_HNL){
	using namespace LeptonInjector_HNL;
	load_project("libLeptonInjector_HNL", false);
	load_project("libtableio_HNL", false);
	
	class_<BasicInjectionConfiguration,bases<I3FrameObject> >("BasicInjectionConfiguration")
	.def_readonly("events",&BasicInjectionConfiguration::events)
	.def_readonly("energyMinimum",&BasicInjectionConfiguration::energyMinimum)
	.def_readonly("energyMaximum",&BasicInjectionConfiguration::energyMaximum)
	.def_readonly("powerlawIndex",&BasicInjectionConfiguration::powerlawIndex)
	.def_readonly("azimuthMinimum",&BasicInjectionConfiguration::azimuthMinimum)
	.def_readonly("azimuthMaximum",&BasicInjectionConfiguration::azimuthMaximum)
	.def_readonly("zenithMinimum",&BasicInjectionConfiguration::zenithMinimum)
	.def_readonly("zenithMaximum",&BasicInjectionConfiguration::zenithMaximum)
	.def_readonly("finalType1",&BasicInjectionConfiguration::finalType1)
	.def_readonly("finalType2",&BasicInjectionConfiguration::finalType2)
	.def_readonly("HNL_mass",&BasicInjectionConfiguration::HNL_mass)
	;
	
	class_<RangedInjectionConfiguration,bases<BasicInjectionConfiguration,I3FrameObject> >("RangedInjectionConfiguration")
	.def_readonly("injectionRadius",&RangedInjectionConfiguration::injectionRadius)
	.def_readonly("endcapLength",&RangedInjectionConfiguration::endcapLength)
	;
	
	class_<VolumeInjectionConfiguration,bases<BasicInjectionConfiguration,I3FrameObject> >("VolumeInjectionConfiguration")
	.def_readonly("cylinderRadius",&VolumeInjectionConfiguration::cylinderRadius)
	.def_readonly("cylinderHeight",&VolumeInjectionConfiguration::cylinderHeight)
	;
	
	class_<MinimalInjectionConfiguration, boost::shared_ptr<MinimalInjectionConfiguration> >("injector",
	  init<unsigned int,I3Particle::ParticleType,I3Particle::ParticleType,std::string,std::string,bool>(
	    (args("NEvents"),args("FinalType1"),args("FinalType2"),args("DoublyDifferentialCrossSectionFile"),args("TotalCrossSectionFile"),args("Ranged"))
	  )
	)
	.def_readwrite("events",&MinimalInjectionConfiguration::events)
	.def_readwrite("finalType1",&MinimalInjectionConfiguration::finalType1)
	.def_readwrite("finalType2",&MinimalInjectionConfiguration::finalType2)
	.def_readwrite("crossSectionPath",&MinimalInjectionConfiguration::crossSectionPath)
	.def_readwrite("totalCrossSectionPath",&MinimalInjectionConfiguration::totalCrossSectionPath)
	.def_readwrite("ranged",&MinimalInjectionConfiguration::ranged)
	;
	
	class_<std::vector<MinimalInjectionConfiguration> >("MultiConfigList")
	.def(vector_indexing_suite<std::vector<MinimalInjectionConfiguration> >())
	;
	from_python_sequence<std::vector<MinimalInjectionConfiguration>, variable_capacity_policy>();
	
	class_<BasicEventProperties,bases<I3FrameObject> >("BasicEventProperties")
	.def_readonly("totalEnergy",&BasicEventProperties::totalEnergy)
	.def_readonly("zenith",&BasicEventProperties::zenith)
	.def_readonly("azimuth",&BasicEventProperties::azimuth)
	.def_readonly("finalStateX",&BasicEventProperties::finalStateX)
	.def_readonly("finalStateY",&BasicEventProperties::finalStateY)
	.def_readonly("finalType1",&BasicEventProperties::finalType1)
	.def_readonly("finalType2",&BasicEventProperties::finalType2)
	.def_readonly("initialType",&BasicEventProperties::initialType)
	.def_readonly("mHNL",&BasicEventProperties::mHNL)
	.def_readonly("cylinderCenter",&BasicEventProperties::cylinderCenter)
	.def_readonly("distance",&BasicEventProperties::distance)
	.def_readonly("distanceMin",&BasicEventProperties::distanceMin)
	.def_readonly("distanceMax",&BasicEventProperties::distanceMax)
	.def_readonly("lifetime",&BasicEventProperties::lifetime)
    .def_readonly("decay_channel",&BasicEventProperties::decay_channel)
    .def_readonly("outgoingNeutrinoEnergy",&BasicEventProperties::outgoingNeutrinoEnergy)
    .def_readonly("physical",&BasicEventProperties::physical)
	;
	
	class_<RangedEventProperties,bases<BasicEventProperties,I3FrameObject> >("RangedEventProperties")
	.def_readonly("impactParameter",&RangedEventProperties::impactParameter)
	.def_readonly("totalColumnDepth",&RangedEventProperties::totalColumnDepth)
	;
	
	class_<VolumeEventProperties,bases<BasicEventProperties,I3FrameObject> >("VolumeEventProperties")
	.def_readonly("radius",&VolumeEventProperties::radius)
	.def_readonly("z",&VolumeEventProperties::z)
	.def_readonly("totalColumnDepth",&VolumeEventProperties::totalColumnDepth)
	;


	register_HNLDecayEnergy();

	{
	I3CONVERTER_NAMESPACE(LeptonInjector_HNL);
	I3CONVERTER_EXPORT_DEFAULT(EventPropertiesConverter,"Converts an EventProperties");
	}
}
