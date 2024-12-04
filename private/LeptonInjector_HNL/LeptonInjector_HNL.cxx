#include <LeptonInjector_HNL/LeptonInjector_HNL.h>

#include <cassert>
#include <fstream>
#include <tuple>
#include <utility> // std::move

#include <boost/math/constants/constants.hpp>
#include <boost/make_shared.hpp>

#include <icetray/open.h>
#include <icetray/I3Units.h>
#include <dataclasses/physics/I3MCTree.h>
#include <LeptonInjector_HNL/HNLDecayEnergy.h>
#include <LeptonInjector_HNL/HNLCrossSection.h>


namespace constants = boost::math::constants;

namespace LeptonInjector_HNL{
	
	//------------------------------------------------
	//Helper functions for dealing with particle types
	
	bool isLepton(I3Particle::ParticleType p){
		return(p==I3Particle::EMinus   || p==I3Particle::EPlus ||
			   p==I3Particle::MuMinus  || p==I3Particle::MuPlus ||
			   p==I3Particle::TauMinus || p==I3Particle::TauPlus ||
			   p==I3Particle::NuE      || p==I3Particle::NuEBar ||
			   p==I3Particle::NuMu     || p==I3Particle::NuMuBar ||
			   p==I3Particle::NuTau    || p==I3Particle::NuTauBar ||
			   p==I3Particle::HNL	   || p==I3Particle::HNLBar);
	}
	
	bool isCharged(I3Particle::ParticleType p){
		assert(isLepton(p) || p==I3Particle::Hadrons);
		return(p==I3Particle::EMinus   || p==I3Particle::EPlus ||
			   p==I3Particle::MuMinus  || p==I3Particle::MuPlus ||
			   p==I3Particle::TauMinus || p==I3Particle::TauPlus ||
			   p==I3Particle::Hadrons);
	}
	
	std::string particleName(I3Particle::ParticleType p){
		return(I3Particle(I3Particle::Null,p).GetTypeString());
	}
	
	double particleMass(I3Particle::ParticleType type){
		I3Particle p(I3Particle::Null,type);
		if(!p.HasMass()){
			log_debug_stream("Treating particle of type " << p.GetTypeString()
							 << " with unknown mass as massless");
			return(0);
		}
		return(p.GetMass());
	}

	std::vector<char> dumpSpline(const photospline::splinetable<>& spline){
		auto buf = spline.write_fits_mem();
		std::vector<char> blob(buf.second);
		std::copy((char*)buf.first,(char*)buf.first+buf.second,blob.data());
		free(buf.first);
		return blob;
	}
	
	double kineticEnergy(I3Particle::ParticleType type, double totalEnergy){
		double mass=particleMass(type);
		if(totalEnergy<mass){
			log_warn_stream("Negative kinetic energy (particle type = " << particleName(type) << ", mass = " << mass << ", total energy = " << totalEnergy << ')');
			return(0.);
		}
		return(sqrt(totalEnergy*totalEnergy-mass*mass));
	}
	
	double kineticEnergy(double mass, double totalEnergy){
		if(totalEnergy<mass){
			log_warn_stream("Negative kinetic energy, mass = " << mass << ", total energy = " << totalEnergy << ')');
			return(0.);
		}
		return(sqrt(totalEnergy*totalEnergy-mass*mass));
	}

	double particleSpeed(I3Particle::ParticleType type, double kineticEnergy){
		I3Particle p=I3Particle(I3Particle::Null,type);
		if(!p.HasMass()){
			log_debug_stream("Treating particle of type " << p.GetTypeString()
							 << " with unknown mass as massless");
			return(I3Constants::c);
		}
		double mass=p.GetMass();
		if(kineticEnergy<0){
			log_warn("Negative kinetic energy");
			return(0.);
		}
		double r=(mass>0 ? mass/(kineticEnergy+mass) : 0.);
		return(I3Constants::c*sqrt(1-r*r));
	}
	
	double particleSpeed(double mass, double kineticEnergy){
		if(kineticEnergy<0){
			log_warn("Negative kinetic energy");
			return(0.);
		}
		double r=(mass>0 ? mass/(sqrt(kineticEnergy*kineticEnergy+mass*mass)) : 0.);
		return(I3Constants::c*sqrt(1-r*r));
	}

	I3Particle::ParticleShape decideShape(I3Particle::ParticleType t){
		switch(t){
			case I3Particle::MuMinus:  case I3Particle::MuPlus:
			case I3Particle::TauMinus: case I3Particle::TauPlus:
			case I3Particle::NuE:      case I3Particle::NuEBar:
			case I3Particle::NuMu:     case I3Particle::NuMuBar:
			case I3Particle::NuTau:    case I3Particle::NuTauBar:
				return(I3Particle::MCTrack);
			case I3Particle::EMinus: case I3Particle::EPlus:
			case I3Particle::Hadrons:
				return(I3Particle::Cascade);
			case I3Particle::HNL: case I3Particle::HNLBar:
				return(I3Particle::Dark);
			default:
				log_fatal_stream("Unable to decide shape for unexpected particle type: " << particleName(t));
		}
	}
	
	I3Particle::ParticleType deduceInitialType(I3Particle::ParticleType pType1, I3Particle::ParticleType pType2){
		
		//Assume all HNL events are produced by tau neutrinos
		if(pType1==I3Particle::HNL && pType2==I3Particle::Hadrons){
			log_info_stream("NuTau -> HNL event");
			return I3Particle::NuTau;
		}
		else if(pType1==I3Particle::HNLBar && pType2==I3Particle::Hadrons){
			log_info_stream("NuTauBar -> HNLBar event");
			return I3Particle::NuTauBar;
        	}

		//only accept certain particle types in general
		if(!isLepton(pType1) && pType1!=I3Particle::Hadrons)
			log_fatal_stream("Unexpected particle type: "
							 << particleName(pType1)
							 << ";\nonly leptons and 'Hadrons' are supported");
		if(!isLepton(pType2) && pType2!=I3Particle::Hadrons)
			log_fatal_stream("Unexpected particle type: "
							 << particleName(pType2)
							 << ";\nonly leptons and 'Hadrons' are supported");
		
		bool c1=isCharged(pType1);
		bool c2=isCharged(pType2);
		bool l1=isLepton(pType1);
		bool l2=isLepton(pType2);
		
		//at least one particle should be charged
		if(!c1 && !c2)
			log_fatal_stream("Final state must contain at least one charged particle\n"
							 << "specified particles were " << particleName(pType1)
							 << " and " << particleName(pType2));
		
		//first particle is charged, second is not
		if(c1 && !c2){
			//valid cases are charged lepton + matching antineutrino for GR
			if(l1){
				//!c2 => pType2 is a neutrino
				if(!((pType1==I3Particle::EMinus   && pType2==I3Particle::NuEBar) ||
					 (pType1==I3Particle::EPlus    && pType2==I3Particle::NuE) ||
					 (pType1==I3Particle::MuMinus  && pType2==I3Particle::NuMuBar) ||
					 (pType1==I3Particle::MuPlus   && pType2==I3Particle::NuMu) ||
					 (pType1==I3Particle::TauMinus && pType2==I3Particle::NuTauBar) ||
					 (pType1==I3Particle::TauPlus  && pType2==I3Particle::NuTau)))
					log_fatal_stream("Final states with a charged lepton must have an anti-matching neutrino.\n"
									 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (leptonic)");
				return(I3Particle::NuEBar);
			}
			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//first particle is neutral, second is charged
		if(!c1 && c2){
			if(l1 && pType2==I3Particle::Hadrons){
				//particle 1 is a neutral lepton, so it must be a neutrino
				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Neutral Current");
				return(pType1); //the incoming neutrino type is the same as the outgoing
			}
			log_fatal_stream("Unrecognized final state type: " << particleName(pType1) << " and " << particleName(pType2));
		}
		
		//have two charged particles
		if(c1 && c2){
			//no two charged lepton states
			if(l1 && l2)
				log_fatal_stream("Two charged lepton final states are not allowed.\n"
								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			//lepton should be given first
			if(!l1 && l2)
				log_fatal_stream("Final states should specify charged leptons before 'Hadrons'.\n"
								 << "Specified particles were " << particleName(pType1) << " and " << particleName(pType2));
			
			if(l1 && !l2){ //valid: charged lepton + Hadrons for CC
				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Charged Current");
				switch(pType1){
					case I3Particle::EMinus: return(I3Particle::NuE);
					case I3Particle::EPlus: return(I3Particle::NuEBar);
					case I3Particle::MuMinus: return(I3Particle::NuMu);
					case I3Particle::MuPlus: return(I3Particle::NuMuBar);
					case I3Particle::TauMinus: return(I3Particle::NuTau);
					case I3Particle::TauPlus: return(I3Particle::NuTauBar);
					default: assert(false && "This point should be unreachable");
				}
			}
			if(!l1 && !l2){ //valid: two hadrons (for GR)
				log_info_stream(particleName(pType1) << ", " << particleName(pType2) << " identified as Glashow Resonance (hadronic)");
				return(I3Particle::NuEBar);
			}
		}
		log_fatal("Logic error; this point should be unreachable");
	}
	
	//--------------
	//Config objects
	
	BasicInjectionConfiguration::BasicInjectionConfiguration():
	events(1),
	energyMinimum(10*I3Units::GeV),
	energyMaximum(1e9*I3Units::GeV),
	powerlawIndex(1.0),
	azimuthMinimum(0),
	azimuthMaximum(2*constants::pi<double>()),
	zenithMinimum(0),
	zenithMaximum(constants::pi<double>()),
	finalType1(I3Particle::MuMinus),
	finalType2(I3Particle::Hadrons)
	{}
	
	BasicInjectionConfiguration::~BasicInjectionConfiguration(){}
	
	RangedInjectionConfiguration::RangedInjectionConfiguration():
	injectionRadius(1200*I3Units::meter),
	endcapLength(1200*I3Units::meter)
	{}
	
	RangedInjectionConfiguration::~RangedInjectionConfiguration(){}
	
	VolumeInjectionConfiguration::VolumeInjectionConfiguration():
	cylinderRadius(1200*I3Units::meter),
	cylinderHeight(1200*I3Units::meter)
	{}
	
	VolumeInjectionConfiguration::~VolumeInjectionConfiguration(){}
	
	void BasicInjectionConfiguration::setCrossSection(const 		photospline::splinetable<>& crossSection, const photospline::splinetable<>& totalCrossSection){
		crossSectionBlob = dumpSpline(crossSection);
		totalCrossSectionBlob = dumpSpline(totalCrossSection);
	}
	
	//---------------------
	//Event propery objects
	
	BasicEventProperties::~BasicEventProperties(){}
	
	RangedEventProperties::~RangedEventProperties(){}
	
	VolumeEventProperties::~VolumeEventProperties(){}
	
	//-----------
	//Module base
	
	LeptonInjectorBase::LeptonInjectorBase(const I3Context& context,BasicInjectionConfiguration& config):
	I3ConditionalModule(context),
	config(config),
	eventsGenerated(0),
	wroteConfigFrame(false),
	suspendOnCompletion(true){
		//do NOTHING with config in this constructor, as it is not yet fully constructed
		AddOutBox("OutBox");
	}
	
	LeptonInjectorBase::~LeptonInjectorBase(){
	}
	
	void LeptonInjectorBase::AddBaseParameters(){
		AddParameter("NEvents",
					 "Number of events to generate",
					 config.events);
		AddParameter("MinimumEnergy",
					 "Minimum total event energy to inject",
					 config.energyMinimum);
		AddParameter("MaximumEnergy",
					 "Maximum total event energy to inject",
					 config.energyMaximum);
		AddParameter("PowerlawIndex",
					 "Powerlaw index of the energy spectrum to inject "
					 "(should be positive)",
					 config.powerlawIndex);
		AddParameter("MinimumAzimuth",
					 "Minimum azimuth angle for injected events",
					 config.azimuthMinimum);
		AddParameter("MaximumAzimuth",
					 "Maximum azimuth angle for injected events",
					 config.azimuthMaximum);
		AddParameter("MinimumZenith",
					 "Minimum zenith angle for injected events",
					 config.zenithMinimum);
		AddParameter("MaximumZenith",
					 "Maximum zenith angle for injected events",
					 config.zenithMaximum);
		AddParameter("FinalType1",
					 "The first particle type in the final state",
					 config.finalType1);
		AddParameter("FinalType2",
					 "The seocnd particle type in the final state",
					 config.finalType2);
		AddParameter("RandomService",
					 "Name of the random service to use",
					 "I3RandomService");
		AddParameter("DoublyDifferentialCrossSectionFile",
					 "Path to the spline FITS file representing the doubly-differential cross section",
					 "");
		AddParameter("TotalCrossSectionFile",
					 "Path to the spline FITS file representing the total cross section as a function of energy"
					 " (same as DoublyDifferentialCrossSectionFile but integrated over x and y)",
					 "");
		AddParameter("SuspendOnCompletion",
					 "Suspend the tray after all events have been generated",
					 suspendOnCompletion);
		AddParameter("HNL_mass",
					 "Specify a mass for the simulated HNLs",
					 config.HNL_mass);
        	AddParameter("CylinderCenter",
					 "Set the custom center of the cylinder",
					 config.cylinderCenter);
	}
	
	void LeptonInjectorBase::BaseConfigure(){
		std::string randomServiceName;
		std::string dd_crossSectionFile;
		std::string total_crossSectionFile;
		
		GetParameter("NEvents",config.events);
		GetParameter("MinimumEnergy",config.energyMinimum);
		GetParameter("MaximumEnergy",config.energyMaximum);
		GetParameter("PowerlawIndex",config.powerlawIndex);
		GetParameter("MinimumAzimuth",config.azimuthMinimum);
		GetParameter("MaximumAzimuth",config.azimuthMaximum);
		GetParameter("MinimumZenith",config.zenithMinimum);
		GetParameter("MaximumZenith",config.zenithMaximum);
		GetParameter("FinalType1",config.finalType1);
		GetParameter("FinalType2",config.finalType2);
		GetParameter("RandomService",randomServiceName);
		GetParameter("DoublyDifferentialCrossSectionFile",dd_crossSectionFile);
		GetParameter("TotalCrossSectionFile",total_crossSectionFile);
		GetParameter("SuspendOnCompletion",suspendOnCompletion);
		
		GetParameter("HNL_mass", config.HNL_mass); //Added to support single mass samples
		GetParameter("CylinderCenter", config.cylinderCenter);
				
		if(config.events==0)
			log_fatal_stream(GetName() << ": there's no point in running this if you don't generate at least one event");
		if(config.energyMinimum<=0)
			log_fatal_stream(GetName() << ": minimum energy must be positive");
		if(config.energyMaximum<=0)
			log_fatal_stream(GetName() << ": maximum energy must be positive");
		if(config.energyMaximum<config.energyMinimum)
			log_fatal_stream(GetName() << ": maximum energy must be greater than or equal to minimum energy");
		if(config.azimuthMinimum<0.0)
			log_fatal_stream(GetName() << ": minimum azimuth angle must be greater than or equal to zero");
		if(config.azimuthMaximum>2*constants::pi<double>())
			log_fatal_stream(GetName() << ": maximum azimuth angle must be less than or equal to 2 pi");
		if(config.azimuthMinimum>config.azimuthMaximum)
			log_fatal_stream(GetName() << ": minimum azimuth angle must be less than or equal to maximum azimuth angle");
		if(config.zenithMinimum<0.0)
			log_fatal_stream(GetName() << ": minimum zenith angle must be greater than or equal to zero");
		if(config.zenithMaximum>constants::pi<double>())
			log_fatal_stream(GetName() << ": maximum zenith angle must be less than or equal to pi");
		if(config.zenithMinimum>config.zenithMaximum)
			log_fatal_stream(GetName() << ": minimum zenith angle must be less than or equal to maximum zenith angle");
	
		//Added to support single mass samples
		if(config.HNL_mass > config.energyMaximum)
			log_fatal_stream(GetName() << ": HNL mass must be less than the maximum energy");
		try{
			initialType=deduceInitialType(config.finalType1,config.finalType2);
		}catch(std::runtime_error& re){
			log_error_stream("While configuring " << GetName());
			throw;
		}
		random = context_.Get<boost::shared_ptr<I3RandomService> >(randomServiceName);
		if(!random)
			log_fatal_stream(GetName() << ": A random service is required");
		if(dd_crossSectionFile.empty())
			log_fatal_stream(GetName() << ": DoublyDifferentialCrossSectionFile must be specified");
		else if(total_crossSectionFile.empty())
			log_fatal_stream(GetName() << ": TotalCrossSectionFile must be specified");
		else
			crossSection.load(dd_crossSectionFile,total_crossSectionFile);
	}
	
	void LeptonInjectorBase::Finish(){
		if(eventsGenerated!=config.events)
			log_error_stream(GetName() << ": Only " << eventsGenerated <<
							 " event have been output out of a requested total of " << config.events);
	}
	
	I3Position LeptonInjectorBase::SampleFromDisk(double radius, double zenith, double azimuth){
		//choose a random point on a disk laying in the xy plane
		double t=random->Uniform(0,2*constants::pi<double>());
		double u=random->Uniform()+random->Uniform();
		double r=(u>1.?2.-u:u)*radius;
		I3Position pos(r*cos(t),r*sin(t),0.0);
		//now rotate to make the disc perpendicular to the requested normal vector
		pos.RotateY(zenith);
		pos.RotateZ(azimuth);
		return(pos);
	}
	
	double LeptonInjectorBase::SampleEnergy(){
		if(config.energyMinimum==config.energyMaximum)
			return(config.energyMinimum); //return the only allowed energy
			
		double minE = std::max(config.energyMinimum, config.HNL_mass); //prevent chosen energy from being less than the mass
		if(config.powerlawIndex==1.0) //sample uniformly in log space
            		return(pow(10.0,random->Uniform(log10(minE),log10(config.energyMaximum))));
		else{
			double u=random->Uniform();
            		double energyP=(1-u)*pow(minE,1-config.powerlawIndex) + u*pow(config.energyMaximum,1-config.powerlawIndex);
			return(pow(energyP,1/(1-config.powerlawIndex)));
		}
	}
	
	std::pair<double,double> LeptonInjectorBase::computeFinalStateAngles(double E_total, double x, double y){
		const double M_N = crossSection.getTargetMass();
		double theta1=0, theta2=0;
		
		//first particle is a lepton, which covers CC, NC, and leptonic GR
		if(isLepton(config.finalType1)){
			double m1=I3Particle::GetMassForType(config.finalType1);
			double E1 = (1 - y) * E_total;
			double cos_theta1, kE1;

			if (E1<0){
             log_warn_stream("negative energy assigned: y = " << y << ", E = " << E_total);   
            }//I think we don't need this anymore but it's worth checking
			
			if(!isLepton(config.finalType2)){ //CC and NC have Hadrons as particle 2
				//squared kinetic energy of final state particle 1:
				double kE1sq=E1*E1 - m1*m1;
				if(kE1sq<=0)
					log_fatal_stream("Underflow: Squared total lepton energy smaller than squared rest mass"
						"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y << " m=" << m1 << ")");
				cos_theta1=(E1 - x*y*M_N - m1*m1/(2*E_total))/sqrt(kE1sq);
				kE1=sqrt(kE1sq);
			}
			else{ //leptonic GR
				double m_e = I3Particle::GetMassForType(I3Particle::EMinus);
				
				if(E1<=0)
					log_fatal("Underflow: Final state neutrino energy is unphysical");
				
				cos_theta1=1 - (m_e*m_e + 2*m_e*E_total - m1*m1)/(2*E_total*E1);
				kE1=E1;
			}
			
			if(cos_theta1<-1){
				log_warn_stream("cos(theta) underflow (" << cos_theta1 << "); rounding up to -1"
					"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y << ")");
				cos_theta1=-1;
			}
			else if(cos_theta1>1){
				//tell the user if the difference was large enough to plausibly not be just round-off
				if((cos_theta1-1)>1e-3)
					log_warn_stream("cos(theta) overflow (" << cos_theta1 << "); rounding down to 1"
						"\n(E_total=" << E_total/I3Units::GeV << " x=" << x << " y=" << y <<")");
				cos_theta1=1;
			}
			
			theta1=acos(cos_theta1);
			
			//longitudinal component of final state particle 2 momentum:
			double p_long=E_total-kE1*cos_theta1;
			//transverse component of final state particle 2 momentum:
			double p_trans=kE1*sin(theta1);
			theta2=atan(p_trans/p_long);
		}
		//otherwise we have hadronic GR, so both final state masses are unknown
		//and there isn't much we can do, so leave everything colinear
		return(std::make_pair(theta1,theta2));
	}
	
	I3Direction rotateRelative(I3Direction base, double zenith, double azimuth){
		I3Direction result(zenith,azimuth);
		result.RotateY(base.GetZenith());
		result.RotateZ(base.GetAzimuth());
		return(result);
	}
	
	boost::shared_ptr<I3MCTree> LeptonInjectorBase::FillTree(I3Position vertex, I3Direction dir, double energy, BasicEventProperties& properties){
        
		double mHNL = config.HNL_mass;
        
        
		//Implement kinematic bounds on the cross sections
        //following HEP PH 0407371 Eq. 7
        
// 		double proton_mass = 0.938272;
//         double neutron_mass = 0.939565;
//         double mI = (proton_mass + neutron_mass)/2;
//         double d_lepton = mHNL/(2*mI*energy);
        
//         bool physical = false;
//         int iter = 0;
//         const I3CrossSection::finalStateRecord& fs;
        
//         while(physical == false){
//             const I3CrossSection::finalStateRecord& temp_fs = crossSection.sampleFinalState(energy,config.finalType1,random);
//             double h = temp_fs.x*temp_fs.y*d_lepton;
//             if( ((1-temp_fs.y)*energy>mHNL) && (pow(h, 2)*(1 + temp_fs.x*d_lepton) - h*temp_fs.x*d_lepton + temp_fs.x*d_lepton <= 0 ) ){
//                 fs = temp_fs;
//                 physical = true;
//             }
//             if (iter>50){
//                 break;
//             }
//             iter += 1;
                     
//         }
// //         const I3CrossSection::finalStateRecord& fs = temp_fs;
        bool physical = true;
        const HNLCrossSection::finalStateRecord& fs = crossSection.sampleFinalState(mHNL, energy,config.finalType1,random);

		std::pair<double,double> relativeZeniths=computeFinalStateAngles(energy,fs.x,fs.y);
		double azimuth1=random->Uniform(0,2*constants::pi<double>());
		double azimuth2=azimuth1+(azimuth1<constants::pi<double>() ? 1 : -1)*constants::pi<double>();
		//Why are we using constants::pi here and M_PI elsewhere?

		//Confirm HNL mass
		//double mHNL = config.HNL_mass;
		
		//TODO: re-enable continuous mass sampling and interpolation from madgraph files for other masses
        	//Set a random HNL mass from the allowed masses if none is chosen
                if (mHNL == 0.0){
					int random_mass = random -> Integer(4);
					double masses[4] = {0.1, 0.3, 0.6, 1.0};
                    mHNL = masses[random_mass];
                }
                
		// Check that one of the supported masses has been chosen
                if (!(mHNL == 0.1 || mHNL == 0.3 || mHNL == 0.6 || mHNL == 1.0)){
					log_fatal_stream("Unsupported particle mass: " << mHNL << ". Please choose 0.1, 0.3, 0.6, or 1.0 Gev");                                                                           }


		
		//Make the first final state particle
		I3Particle p1(decideShape(config.finalType1),config.finalType1);
		p1.SetLocationType(I3Particle::InIce);
		p1.SetPos(vertex);
		p1.SetDir(rotateRelative(dir,relativeZeniths.first,azimuth1));
		p1.SetEnergy(kineticEnergy(mHNL, (1-fs.y)*energy));
		p1.SetSpeed(particleSpeed(mHNL, p1.GetEnergy()));
		p1.SetTime(0.0);
		
		
		//Make the HNL

		//Choose a length range
		//TODO: Implement length range as a user parameter

		//larger length range, aimed at double cascade detection
		double distanceMin = 1.*I3Units::m;
		double distanceMax = 1000.*I3Units::m;

		//smaller length range to cover single cascade region
		// double distanceMin = 0.0001*I3Units::m;
		// double distanceMax = 1.*I3Units::m;
		
		//Find a lifetime from the distance range

		double gamma_1 = ((1-fs.y)*energy)/mHNL;
		double speed_1 = sqrt(1-pow(1./gamma_1,2))*I3Constants::c;
		double lifetime_min = distanceMin/(gamma_1*speed_1);
		double lifetime_max = distanceMax/(gamma_1*speed_1); 												       
		double x = random->Uniform(0,log(lifetime_max)-log(lifetime_min));
		double lifetime_rest = exp(log(lifetime_min) + x);  // bracket was wrong!
		double lifetime_boosted = gamma_1*lifetime_rest;
		I3Position vertex_daughters = vertex + p1.GetDir()*speed_1*lifetime_boosted;

		//Set random number to choose decay channel
		double rand_channel = random->Uniform(0,1);

		//Random numbers for decay daughters
		double rand_coszenith = random->Uniform(-1,1);
		double rand_azimuth = random->Uniform(0,2*M_PI);
        int rand_event = random->Integer(999999)+1;
		
        bool isHNL = (config.finalType1==I3Particle::HNL);
		
		//Create the daughter(s) vector, get from HNLDecayEnergy routine		
		std::vector<double> daughters_vector = HNLDecayEnergy::GetDaughterVector(isHNL,mHNL,(1-fs.y)*energy,M_PI-p1.GetZenith(),M_PI+p1.GetAzimuth(),rand_channel,rand_coszenith,rand_azimuth,rand_event);
		int n_daughters = (daughters_vector.size()-1)/4;  // Each daughter has 4 elements. Exclude last value, which is the decay channel
            
		//Save decay channel
		double decay_channel = daughters_vector[daughters_vector.size()-1];
		I3Particle p_daughter_1(I3Particle::Cascade,I3Particle::unknown);

		//Initialize daughters 2 and 3 (necessary even if not used)
		I3Particle p_daughter_2(I3Particle::Cascade,I3Particle::unknown);
		I3Particle p_daughter_3(I3Particle::Cascade,I3Particle::unknown);

		//Make daughter particle(s)
		p_daughter_1.SetLocationType(I3Particle::InIce);
        p_daughter_1.SetType(I3Particle::ParticleType(daughters_vector.at(0)));
		p_daughter_1.SetPos(vertex_daughters);
		double azimuth_daughter_1 = ((daughters_vector.at(3)-M_PI)<0) ? (M_PI+daughters_vector.at(3)) : (daughters_vector.at(3)-M_PI);
		p_daughter_1.SetDir(M_PI-daughters_vector.at(2),azimuth_daughter_1);
		p_daughter_1.SetEnergy(kineticEnergy(p_daughter_1.GetMass(),daughters_vector.at(1)));  // is daughters_vector.at(1) the total energy? Yes, boosted to the lab frame
		p_daughter_1.SetSpeed(particleSpeed(p_daughter_1.GetType(),p_daughter_1.GetEnergy()));
		p_daughter_1.SetTime(lifetime_boosted);
        
        if(n_daughters>=2){
			p_daughter_2.SetLocationType(I3Particle::InIce);
			p_daughter_2.SetType(I3Particle::ParticleType(daughters_vector.at(4)));
			p_daughter_2.SetPos(vertex_daughters);
			double azimuth_daughter_2 = ((daughters_vector.at(7)-M_PI)<0) ? (M_PI+daughters_vector.at(7)) : (daughters_vector.at(7)-M_PI);
			p_daughter_2.SetDir(M_PI-daughters_vector.at(6),azimuth_daughter_2);
			p_daughter_2.SetEnergy(kineticEnergy(p_daughter_2.GetMass(),daughters_vector.at(5)));  // is daughters_vector.at(5) the total energy?
			p_daughter_2.SetSpeed(particleSpeed(p_daughter_2.GetType(),p_daughter_2.GetEnergy()));
			p_daughter_2.SetTime(lifetime_boosted);
		}
		if(n_daughters==3){
			p_daughter_3.SetLocationType(I3Particle::InIce);
			p_daughter_3.SetType(I3Particle::ParticleType(daughters_vector.at(8)));
			p_daughter_3.SetPos(vertex_daughters);
			double azimuth_daughter_3 = ((daughters_vector.at(11)-M_PI)<0) ? (M_PI+daughters_vector.at(11)) : (daughters_vector.at(11)-M_PI);
			p_daughter_3.SetDir(M_PI-daughters_vector.at(10),azimuth_daughter_3);
			p_daughter_3.SetEnergy(kineticEnergy(p_daughter_3.GetMass(),daughters_vector.at(9)));  // is daughters_vector.at(9) the total energy?
			p_daughter_3.SetSpeed(particleSpeed(p_daughter_3.GetType(),p_daughter_3.GetEnergy()));
			p_daughter_3.SetTime(lifetime_boosted);
		}
		//Making the "Hadrons" object

		//Make the second final state particle
		I3Particle p2(decideShape(config.finalType2),config.finalType2);
		p2.SetLocationType(I3Particle::InIce);
		p2.SetPos(vertex);
		p2.SetDir(rotateRelative(dir,relativeZeniths.second,azimuth2));
		p2.SetEnergy(kineticEnergy(p2.GetType(),fs.y*energy));
		p2.SetSpeed(particleSpeed(p2.GetType(),p2.GetEnergy()));
		p2.SetTime(0.0);
        
		boost::shared_ptr<I3MCTree> mctree(new I3MCTree);

		//Make a dummy primary
		I3Particle primary(I3Particle::Primary,initialType);
		primary.SetEnergy(kineticEnergy(primary.GetType(),energy));
		primary.SetSpeed(particleSpeed(primary.GetType(),primary.GetEnergy()));
		primary.SetDir(dir);
		primary.SetPos(vertex);
		primary.SetTime(0.0);
		primary.SetLength(0.0);
		I3MCTree::iterator primaryIt=mctree->insert(mctree->end(), primary);
		I3MCTree::iterator p1It=mctree->append_child(primaryIt, p1);
		I3MCTree::iterator p_daughter_1It = mctree->append_child(p1It, p_daughter_1);
		if(n_daughters==2) mctree->insert_after(p_daughter_1It, p_daughter_2);
		if(n_daughters==3) {
			mctree->insert_after(p_daughter_1It, p_daughter_3);
			mctree->insert_after(p_daughter_1It, p_daughter_2);
		}		
		mctree->insert_after(p1It, p2);
		properties.totalEnergy=energy;
		properties.zenith=dir.GetZenith();
		properties.azimuth=dir.GetAzimuth();
		properties.finalStateX=fs.x;
		properties.finalStateY=fs.y;
		properties.finalType1=p1.GetType(); //Should give 19 (HNL)
		properties.finalType2=p2.GetType(); //Should give -2e9 (Hadrons object)
		properties.initialType=primary.GetType();
	
		// Keep track of which particle is the outgoing neutrino
		// I'm fairly certain this does not behave as expected for the 3-body decay case.
		int neutrino_index = -1;
		for (int i=0; i<n_daughters; i++){
			if(std::abs(daughters_vector.at(i*4)) == 12 || std::abs(daughters_vector.at(i*4))==14 || std::abs(daughters_vector.at(i*4))==16){
				neutrino_index = i*4;
				break;
			}
		}

		properties.outgoingNeutrinoEnergy = daughters_vector.at(neutrino_index+1);
                    
		properties.mHNL=mHNL;
		properties.cylinderCenter=config.cylinderCenter;
		properties.distance=speed_1*lifetime_boosted;
		properties.distanceMin=distanceMin;
		properties.distanceMax=distanceMax;
		properties.lifetime=lifetime_rest;
        properties.decay_channel=decay_channel;
		properties.physical=physical;
		
		return(mctree);
	}
	
	//-----------------------
	//Ranged injection module
	
	RangedLeptonInjector::RangedLeptonInjector(const I3Context& context):
	LeptonInjectorBase(context,config){
		init();
	}
	
	RangedLeptonInjector::RangedLeptonInjector(const I3Context& context, RangedInjectionConfiguration config_):
	LeptonInjectorBase(context,config),config(config_){
		init();
	}
	
	void RangedLeptonInjector::init(){
		AddBaseParameters();
		AddParameter("InjectionRadius",
					 "Radius around the origin within which to target events",
					 config.injectionRadius);
		AddParameter("EndcapLength",
					 "Length of the fixed endcaps add to each end of the distance "
					 "along which to sample interactions",
					 config.endcapLength);
		AddParameter("EarthModel",
					 "Name of the Earth model service to use",
					 "");
		AddOutBox("OutBox");
	}
	
	void RangedLeptonInjector::Configure(){
		BaseConfigure();
		GetParameter("InjectionRadius",config.injectionRadius);
		GetParameter("EndcapLength",config.endcapLength);
		std::string earthModelName;
		GetParameter("EarthModel",earthModelName);
		
		if(config.injectionRadius<0)
			log_fatal_stream(GetName() << ": InjectionRadius must be non-negative");
		if(config.endcapLength<0)
			log_fatal_stream(GetName() << ": EndcapLength must be non-negative");
		earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
		if(!earthModel)
			log_fatal_stream(GetName() << ": an Earth model service is required");
	}
	
	void RangedLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		//first, make sure configuration gets written once
		if(!wroteConfigFrame){
			boost::shared_ptr<I3Frame> sframe(new I3Frame('S'));
            boost::shared_ptr<RangedInjectionConfiguration> sconfig(new RangedInjectionConfiguration(config));
			sconfig->setCrossSection(getCrossSection(),getTotalCrossSection());
			sframe->Put("LeptonInjectorProperties",sconfig);
			PushFrame(sframe);
			wroteConfigFrame=true;
		}
		if(DoneGenerating()){
			PushFrame(frame);
            std::cout<<"Frame pushed"<<std::endl;
			return;
		}
		
		//Choose an energy
		double energy=SampleEnergy();
		
		//Pick a direction on the sphere
		I3Direction dir(acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
						random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
		log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
		
		//decide the point of closest approach
		I3Position pca=SampleFromDisk(config.injectionRadius,dir.GetZenith(),dir.GetAzimuth());
		log_trace_stream("pca=(" << pca.GetX() << ',' << pca.GetY() << ',' << pca.GetZ() << ')');
		
		//Figure out where we want the vertex
		//Add up the column depth for the range of a muon at this energy with the
		//column depth for the fixed endcaps to ensure that the whole detector is
		//covered
		using namespace earthmodel::EarthModelCalculator;
		
		bool is_tau = (config.finalType1 == I3Particle::TauMinus) or (config.finalType1 == I3Particle::TauPlus);
		double lepton_range = GetLeptonRange(energy, is_tau);
		
		
		double totalColumnDepth=MWEtoColumnDepthCGS(lepton_range)
		+earthModel->GetColumnDepthInCGS(pca-config.endcapLength*dir,pca+config.endcapLength*dir);
		//See whether that much column depth actually exists along the chosen path
		{
			double maxDist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth)-config.endcapLength;
			double actualColumnDepth=earthModel->GetColumnDepthInCGS(pca+config.endcapLength*dir,pca-maxDist*dir);
			if(actualColumnDepth<(totalColumnDepth-1)){ //if actually smaller, clip as needed, but for tiny differences we don't care
				log_debug_stream("Wanted column depth of " << totalColumnDepth << " but found only " << actualColumnDepth << " g/cm^2");
				totalColumnDepth=actualColumnDepth;
			}
		}
		//Choose how much of the total column depth this event should have to traverse
		double traversedColumnDepth=totalColumnDepth*random->Uniform();
		//endcapLength is subtracted so that dist==0 corresponds to pca
		double dist=earthModel->DistanceForColumnDepthToPoint(pca+config.endcapLength*dir,dir,totalColumnDepth-traversedColumnDepth)-config.endcapLength;
		
		{ //ensure that the point we picked is inside the atmosphere
			I3Position atmoEntry, atmoExit;
			int isect=GetIntersectionsWithSphere(earthModel->GetEarthCoordPosFromDetCoordPos(pca),
												 earthModel->GetEarthCoordDirFromDetCoordDir(dir),
												 earthModel->GetAtmoRadius(),atmoEntry,atmoExit);
			if(isect<2)
				log_fatal_stream("PCA not inside atmosphere: " << pca << " (" << earthModel->GetEarthCoordPosFromDetCoordPos(pca) << ')');
			atmoEntry=earthModel->GetDetCoordPosFromEarthCoordPos(atmoEntry);
			double atmoDist=(pca-atmoEntry).Magnitude();
			if(std::abs(dist-atmoDist)<100.0)
				dist=std::min(dist,atmoDist);
		}
		I3Position vertex=pca-dist*dir;
		
		//assemble the MCTree
		boost::shared_ptr<RangedEventProperties> properties(new RangedEventProperties);
		boost::shared_ptr<I3MCTree> mctree=FillTree(vertex,dir,energy,*properties);
		
		//set subclass properties
		properties->impactParameter=(pca-I3Position(0,0,0)).Magnitude();
		properties->totalColumnDepth=totalColumnDepth;
		
		//package up output and send it
		frame->Put(mctree);
		frame->Put("EventProperties",properties);
		PushFrame(frame);
		
		//update event count and check for completion
		eventsGenerated++;
		if(eventsGenerated==config.events && suspendOnCompletion)
			RequestSuspension();
	}
	
	I3_MODULE(RangedLeptonInjector);
	
	//-----------------------
	//Volume injection module
	
	VolumeLeptonInjector::VolumeLeptonInjector(const I3Context& context):
	LeptonInjectorBase(context,config){
		init();
	}
	
	VolumeLeptonInjector::VolumeLeptonInjector(const I3Context& context, VolumeInjectionConfiguration config_):
	LeptonInjectorBase(context,config),config(config_){
		init();
	}
	
	void VolumeLeptonInjector::init(){
		AddBaseParameters();
		AddParameter("CylinderRadius",
					 "Radius of the vertical cylinder around the origin within "
					 "which to place events",
					 config.cylinderRadius);
		AddParameter("CylinderHeight",
					 "Height of the vertical cylinder around the origin within "
					 "which to place events",
					 config.cylinderHeight);
		AddParameter("EarthModel",
					 "Name of the Earth model service to use",
					 "");
		AddOutBox("OutBox");
	}
	
	void VolumeLeptonInjector::Configure(){
		BaseConfigure();
		GetParameter("CylinderRadius",config.cylinderRadius);
		GetParameter("CylinderHeight",config.cylinderHeight);
		std::string earthModelName;
		GetParameter("EarthModel",earthModelName);
		if(config.cylinderRadius<0)
			log_fatal_stream(GetName() << ": CylinderRadius must be non-negative");
		if(config.cylinderHeight<0)
			log_fatal_stream(GetName() << ": CylinderHeight must be non-negative");
		earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
		if(!earthModel)
			log_fatal_stream(GetName() << ": an Earth model service is required");
	}
	
	void VolumeLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		//first, make sure configuration gets written once
		
		if(!wroteConfigFrame){
			boost::shared_ptr<I3Frame> sframe(new I3Frame('S'));
			boost::shared_ptr<VolumeInjectionConfiguration> sconfig(new VolumeInjectionConfiguration(config));
			sconfig->setCrossSection(getCrossSection(),getTotalCrossSection());
			sframe->Put("LeptonInjectorProperties",sconfig);
			PushFrame(sframe);
			
			wroteConfigFrame=true;
		}
		if(DoneGenerating()){
			PushFrame(frame);
			
			return;
		}
		
		//Choose an energy
		double energy=SampleEnergy();
		
// 		//Pick a direction on the sphere
// 		I3Direction dir(acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum))),
// 						random->Uniform(config.azimuthMinimum,config.azimuthMaximum));
// 		log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
        double zenith = acos(random->Uniform(cos(config.zenithMaximum),cos(config.zenithMinimum)));
        double azimuth = random->Uniform(config.azimuthMinimum,config.azimuthMaximum);
                    
        I3Direction dir(zenith, azimuth);
        
        
		log_trace_stream("dir=(" << dir.GetX() << ',' << dir.GetY() << ',' << dir.GetZ() << ')');
        

		//Pick a position in the xy-plane
		I3Position vertex=SampleFromDisk(config.cylinderRadius);
		//Add on the vertical component
		vertex.SetZ(random->Uniform(config.cylinderCenter - config.cylinderHeight/2,config.cylinderCenter + config.cylinderHeight/2));
		log_trace_stream("vtx=(" << vertex.GetX() << ',' << vertex.GetY() << ',' << vertex.GetZ() << ')');
		
		//assemble the MCTree
		boost::shared_ptr<VolumeEventProperties> properties(new VolumeEventProperties);
		boost::shared_ptr<I3MCTree> mctree=FillTree(vertex,dir,energy,*properties);
		
		//set subclass properties
		properties->radius=vertex.GetRho();
		properties->z=vertex.GetZ();
        std::tuple<I3Position, I3Position> cylinder_intersections =
            computeCylinderIntersections(vertex, dir, config.cylinderRadius, config.cylinderCenter - config.cylinderHeight/2., config.cylinderCenter + config.cylinderHeight/2.);
        properties->totalColumnDepth =
            earthModel->GetColumnDepthInCGS(std::get<0>(cylinder_intersections), std::get<1>(cylinder_intersections));
			
        std::ofstream outputFile; // Declare an output file stream object
        // Open the file in append mode
        outputFile.open("/data/user/jbook/I3_HNL_Search/validation/initial_position.txt", std::ios::app); 
        outputFile << zenith << ", " << azimuth << ", " << dir.GetX() << ',' << dir.GetY() << ',' << vertex.GetZ()<<std::endl; // Write the new data to the file
        outputFile.close(); // Close the file
        
		//package up output and send it
		frame->Put(mctree);
		frame->Put("EventProperties",properties);
		PushFrame(frame);
		
		//update event count and check for completion
		eventsGenerated++;
		if(eventsGenerated==config.events && suspendOnCompletion)
			RequestSuspension();
	}
	
	I3_MODULE(VolumeLeptonInjector);
	
	bool operator==(const MinimalInjectionConfiguration& c1, const MinimalInjectionConfiguration& c2){
		return(c1.events==c2.events && c1.finalType1==c2.finalType1 && c1.finalType2==c2.finalType2
			   && c1.crossSectionPath==c2.crossSectionPath && c1.ranged==c2.ranged);
	}
	
	MultiLeptonInjector::MultiLeptonInjector(const I3Context& ctx):
	I3ConditionalModule(ctx),collector(new OutputCollector(innerContext)),
	results(collector->output){
		AddParameters();
	}
	
	MultiLeptonInjector::MultiLeptonInjector(const I3Context& ctx,
	                                         RangedInjectionConfiguration rconfig,
	                                         VolumeInjectionConfiguration vconfig):
	I3ConditionalModule(ctx),collector(new OutputCollector(innerContext)),
	results(collector->output),rangedConfig(rconfig),volumeConfig(vconfig){
		AddParameters();
	}
	
	void MultiLeptonInjector::AddParameters(){
		AddOutBox("OutBox");
		AddParameter("MinimumEnergy",
					 "Minimum total event energy to inject",
					 rangedConfig.energyMinimum);
		AddParameter("MaximumEnergy",
					 "Maximum total event energy to inject",
					 rangedConfig.energyMaximum);
		AddParameter("PowerlawIndex",
					 "Powerlaw index of the energy spectrum to inject "
					 "(should be positive)",
					 rangedConfig.powerlawIndex);
		AddParameter("MinimumAzimuth",
					 "Minimum azimuth angle for injected events",
					 rangedConfig.azimuthMinimum);
		AddParameter("MaximumAzimuth",
					 "Maximum azimuth angle for injected events",
					 rangedConfig.azimuthMaximum);
		AddParameter("MinimumZenith",
					 "Minimum zenith angle for injected events",
					 rangedConfig.zenithMinimum);
		AddParameter("MaximumZenith",
					 "Maximum zenith angle for injected events",
					 rangedConfig.zenithMaximum);
		AddParameter("RandomService",
					 "Name of the random service to use",
					 "I3RandomService");
		
		AddParameter("InjectionRadius",
					 "Radius around the origin within which to target events",
					 rangedConfig.injectionRadius);
		AddParameter("EndcapLength",
					 "Length of the fixed endcaps add to each end of the distance "
					 "along which to sample interactions",
					 rangedConfig.endcapLength);
		AddParameter("EarthModel",
					 "Name of the Earth model service to use",
					 "");
		
		AddParameter("CylinderRadius",
					 "Radius of the vertical cylinder around the origin within "
					 "which to place events",
					 volumeConfig.cylinderRadius);
		AddParameter("CylinderHeight",
					 "Height of the vertical cylinder around the origin within "
					 "which to place events",
					 volumeConfig.cylinderHeight);
		
		AddParameter("Generators","The collection of configurations to generate",generatorSettings);
		AddParameter("HNL_mass",
					 "Specify a mass for the simulated HNLs",
					 volumeConfig.HNL_mass);
		AddParameter("CylinderCenter",
					 "Custom vertical center of the cylinder",
					 volumeConfig.cylinderCenter);
	}
	
	void MultiLeptonInjector::Configure(){
		try{
			//get the set of generators to be run
			GetParameter("Generators",generatorSettings);
			
			if(generatorSettings.empty())
				log_fatal_stream(GetName()+": There is no point in running this module without specifying at least one generator");
			
			//figure out whether there are any ranged or volume injectors
			bool hasRanged=false, hasVolume=false;
			for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=generatorSettings.begin(), end=generatorSettings.end(); genSet!=end; genSet++){
				hasRanged |= genSet->ranged;
				hasVolume |= !genSet->ranged;
			}
			
			//get the properties shared by all generators
			std::string randomServiceName;
			boost::shared_ptr<I3RandomService> random;
			//fetch each parameter directly into one configuration object,
			//and clone it into the other
			GetParameter("MinimumEnergy",rangedConfig.energyMinimum);
			volumeConfig.energyMinimum=rangedConfig.energyMinimum;
			GetParameter("MaximumEnergy",rangedConfig.energyMaximum);
			volumeConfig.energyMaximum=rangedConfig.energyMaximum;
			GetParameter("PowerlawIndex",rangedConfig.powerlawIndex);
			volumeConfig.powerlawIndex=rangedConfig.powerlawIndex;
			GetParameter("MinimumAzimuth",rangedConfig.azimuthMinimum);
			volumeConfig.azimuthMinimum=rangedConfig.azimuthMinimum;
			GetParameter("MaximumAzimuth",rangedConfig.azimuthMaximum);
			volumeConfig.azimuthMaximum=rangedConfig.azimuthMaximum;
			GetParameter("MinimumZenith",rangedConfig.zenithMinimum);
			volumeConfig.zenithMinimum=rangedConfig.zenithMinimum;
			GetParameter("MaximumZenith",rangedConfig.zenithMaximum);
			volumeConfig.zenithMaximum=rangedConfig.zenithMaximum;
			GetParameter("RandomService",randomServiceName);
			GetParameter("HNL_mass", rangedConfig.HNL_mass);
            volumeConfig.HNL_mass = rangedConfig.HNL_mass;
            GetParameter("CylinderCenter", rangedConfig.cylinderCenter);
            volumeConfig.cylinderCenter = rangedConfig.cylinderCenter;
			
			if(rangedConfig.energyMinimum<=0)
				log_fatal_stream(GetName() << ": minimum energy must be positive");
			if(rangedConfig.energyMaximum<=0)
				log_fatal_stream(GetName() << ": maximum energy must be positive");
			if(rangedConfig.energyMaximum<rangedConfig.energyMinimum)
				log_fatal_stream(GetName() << ": maximum energy must be greater than or equal to minimum energy");
			if(rangedConfig.azimuthMinimum<0.0)
				log_fatal_stream(GetName() << ": minimum azimuth angle must be greater than or equal to zero");
			if(rangedConfig.azimuthMaximum>2*constants::pi<double>())
				log_fatal_stream(GetName() << ": maximum azimuth angle must be less than or equal to 2 pi");
			if(rangedConfig.azimuthMinimum>rangedConfig.azimuthMaximum)
				log_fatal_stream(GetName() << ": minimum azimuth angle must be less than or equal to maximum azimuth angle");
			if(rangedConfig.zenithMinimum<0.0)
				log_fatal_stream(GetName() << ": minimum zenith angle must be greater than or equal to zero");
			if(rangedConfig.zenithMaximum>constants::pi<double>())
				log_fatal_stream(GetName() << ": maximum zenith angle must be less than or equal to pi");
			if(rangedConfig.zenithMinimum>rangedConfig.zenithMaximum)
				log_fatal_stream(GetName() << ": minimum zenith angle must be less than or equal to maximum zenith angle");
			random = context_.Get<boost::shared_ptr<I3RandomService> >(randomServiceName);
			if(!random)
				log_fatal_stream(GetName() << ": A random service is required");
			
			innerContext.Put(random,randomServiceName);
			
			//get the properties for ranged injectors
			std::string earthModelName;
			boost::shared_ptr<earthmodel::EarthModelService> earthModel;
			if(hasRanged){
				GetParameter("InjectionRadius",rangedConfig.injectionRadius);
				GetParameter("EndcapLength",rangedConfig.endcapLength);
				GetParameter("EarthModel",earthModelName);
				
				if(rangedConfig.injectionRadius<0)
					log_fatal_stream(GetName() << ": InjectionRadius must be non-negative");
				if(rangedConfig.endcapLength<0)
					log_fatal_stream(GetName() << ": EndcapLength must be non-negative");
				earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
				if(!earthModel)
					log_fatal_stream(GetName() << ": an Earth model service is required");
				if(!innerContext.Has(earthModelName))
				    innerContext.Put(earthModel,earthModelName);
			}
			
			//get the properties for volume injectors
			if(hasVolume){
				GetParameter("CylinderRadius",volumeConfig.cylinderRadius);
				GetParameter("CylinderHeight",volumeConfig.cylinderHeight);
				GetParameter("EarthModel",earthModelName);
				
				if(volumeConfig.cylinderRadius<0)
					log_fatal_stream(GetName() << ": CylinderRadius must be non-negative");
				if(volumeConfig.cylinderHeight<0)
					log_fatal_stream(GetName() << ": CylinderHeight must be non-negative");
				earthModel = context_.Get<boost::shared_ptr<earthmodel::EarthModelService> >(earthModelName);
				if(!earthModel)
					log_fatal_stream(GetName() << ": an Earth model service is required");
				if(!innerContext.Has(earthModelName))
				    innerContext.Put(earthModel,earthModelName);
			}
// 			std::cout<<"Constructing Generators"<<std::endl;
			//construct all generators
			unsigned int i=0;
			for(std::vector<MinimalInjectionConfiguration>::const_iterator genSet=generatorSettings.begin(), end=generatorSettings.end(); genSet!=end; genSet++){
				log_debug_stream("Configuring injector " << i << ":");
                boost::movelib::unique_ptr<LeptonInjectorBase> generator=NULL;
				try{
					if(genSet->ranged){
						log_debug_stream(" this is a ranged injector");
						generator= boost::movelib::unique_ptr<RangedLeptonInjector>(new RangedLeptonInjector(innerContext,rangedConfig));
						generator->GetConfiguration().Set("EarthModel",boost::python::object(earthModelName));
					}
					else{ //volume
						log_debug_stream(" this is a volume injector");
						generator= boost::movelib::unique_ptr<VolumeLeptonInjector>(new VolumeLeptonInjector(innerContext,volumeConfig));
						generator->GetConfiguration().Set("EarthModel",boost::python::object(earthModelName));
					}
					
					//set properties not shared with other injectors, or which are not part of the config object
					generator->GetConfiguration().Set("NEvents",boost::python::object(genSet->events));
					generator->GetConfiguration().Set("FinalType1",boost::python::object(genSet->finalType1));
					generator->GetConfiguration().Set("FinalType2",boost::python::object(genSet->finalType2));
					generator->GetConfiguration().Set("RandomService",boost::python::object(randomServiceName));
					generator->GetConfiguration().Set("DoublyDifferentialCrossSectionFile",boost::python::object(genSet->crossSectionPath));
					generator->GetConfiguration().Set("TotalCrossSectionFile",boost::python::object(genSet->totalCrossSectionPath));
					generator->GetConfiguration().Set("SuspendOnCompletion",boost::python::object(false));
					
					generator->SetName(GetName()+"_Generator_"+boost::lexical_cast<std::string>(i++));
					generator->Configure();
				}catch(...){
					generator.reset(); //deletes the object contained in the unique ptr
					throw;
				}
                // need to explicitly /move/ the unique ptr to the generators deque
                //      so it doesn't try to copy anything! 
				generators.push_back(std::move(generator));
			}
			
			//bind the first generator to the collector
			generators.front()->ConnectOutBox("OutBox",collector);
		}catch(std::runtime_error& err){
			throw std::runtime_error("While configuring "+GetName()+":\n"+err.what());
		}
	}
	
	void MultiLeptonInjector::DAQ(boost::shared_ptr<I3Frame> frame){
		if(generators.empty()){
			RequestSuspension();
			return;
		}
		while(generators.front()->DoneGenerating()){
			generators.front().reset();
			generators.pop_front();
			//if there are no more generators, we are done
			if(generators.empty()){
				RequestSuspension();
				return;
			}
			//bind the next generator to the collector
			generators.front()->ConnectOutBox("OutBox",collector);
		}
//         std::cout<<"1145"<<std::endl;
		generators.front()->DAQ(frame);
		collector->Process();
        
//         std::cout<<"1149"<<std::endl;
        
		while(!collector->output.empty()){
//             std::cout<<"1152"<<std::endl;
			PushFrame(collector->output.front());
//             std::cout<<"1154"<<std::endl;
			collector->output.pop();
		}
	}

    std::tuple<I3Position, I3Position> computeCylinderIntersections(const I3Position& pos, const I3Direction& dir, double radius, double z_min, double z_max) {
        double x0 = pos.GetX();
        double y0 = pos.GetY();
        double z0 = pos.GetZ();
        double nx = dir.GetX();
        double ny = dir.GetY();
        double nz = dir.GetZ();
        double r = radius;
        double cz1 = z_min;
        double cz2 = z_max;
        // cylinder centered at 0,0,0
        // Assume that the position is within the cylinder
        //     i.e. sqrt(x0**2+y0**2)<=r and z<=cz2 and z>=cz1
        // Assume sqrt(nx**2+ny**2+nz**2)==1
        // x0,y0,z0: particle position relative to cylinder center
        // nx,ny,nz: particle direction vector
        // r: cylinder radius
        // cz1,cz2: cylinder z caps

        assert(z0 <= cz2);
        assert(z0 >= cz1);

        // Prep for the cylinder side intersection
        double nx2 = nx*nx;
        double ny2 = ny*ny;
        double nr2 = 2.*(nx2 + ny2);
        double n_sum = -2.*(nx*x0 + ny*y0);
        double r0_2 = x0*x0 + y0*y0;
        assert(r0_2 <= r*r);

        // Check if the direction is vertical
        if(nx == 0.0 and ny == 0.0) {
            assert(nz != 0.0);
            if(nz > 0.0)
                return std::tuple<I3Position, I3Position>(I3Position(x0, y0, cz1), I3Position(x0, y0, cz2));
            else
                return std::tuple<I3Position, I3Position>(I3Position(x0, y0, cz2), I3Position(x0, y0, cz1));
        }

        // Solving the quadratic
        double root = sqrt(n_sum*n_sum - 2. * nr2 * (r0_2 - r*r));

        // The two solutions
        double sol_1 = (n_sum - root) / nr2;
        double sol_2 = (n_sum + root) / nr2;

        // Corresponding positions
        double x1 = x0 + nx * sol_1;
        double y1 = y0 + ny * sol_1;
        double z1 = z0 + nz * sol_1;
        double x2 = x0 + nx * sol_2;
        double y2 = y0 + ny * sol_2;
        double z2 = z0 + nz * sol_2;

        // Check if the solutions are within the z boundaries
        bool b1_lower = z1 < cz1;
        bool b2_lower = z2 < cz1;

        bool b1_upper = z1 > cz2;
        bool b2_upper = z2 > cz2;

        bool bb_lower = b1_lower or b2_lower;
        bool bb_upper = b1_upper or b2_upper;
        bool bb = bb_lower or bb_upper;

        // Replace with endcap intersections otherwise
        if(bb) {
            if(bb_lower) {
                double t1 = (cz1 - z0)/nz;
                double xx = x0 + nx*t1;
                double yy = y0 + ny*t1;
                double zz = cz1;
                if(b1_lower) {
                    x1 = xx;
                    y1 = yy;
                    z1 = zz;
                }
                else {
                    x2 = xx;
                    y2 = yy;
                    z2 = zz;
                }
            }

            if(bb_upper) {
                double t2 = (cz2 - z0)/nz;
                double xx = x0 + nx*t2;
                double yy = y0 + ny*t2;
                double zz = cz2;
                if(b1_upper) {
                    x1 = xx;
                    y1 = yy;
                    z1 = zz;
                }
                else {
                    x2 = xx;
                    y2 = yy;
                    z2 = zz;
                }
            }
        }

        return std::tuple<I3Position, I3Position>(I3Position(x1, y1, z1), I3Position(x2, y2, z2));
    }
		
	I3_MODULE(MultiLeptonInjector);
	
} //namespace LeptonInjector
