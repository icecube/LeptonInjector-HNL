#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
#METAPROJECT /data/user/jbook/i3/build
import time
import os
import numpy as np

from icecube.icetray import I3Tray
from icecube import icetray, dataclasses, LeptonInjector_HNL, dataio

from icecube.icetray import I3Units
from icecube.hdfwriter import I3SimHDFWriter

keys_to_extract = [
    ###################
    # EventProperties #
    ###################
    "decay_channel", #integer values corresponsing to HNL decay mode
    "distance",      #decay length of the HNL
    "distanceMax",   #Maximum allowed "distance" value
    "distanceMin",   #Minimum allowed "distance" value
    "finalStateX",   #Bjorken X
    "finalStateY",   #Inelasticity
    "final_state_particle0", #PDG code for final state particles
    "final_state_particle1",
    "primary_type",  #PDG code for initial particles ((-)15 for tau (anti)neutrinos)
    "lifetime",      #HNL lifetime
    "mHNL",          #HNL mass - user parameter - 0.1, 0.3, 0.6, 1.0 GeV
    "outgoing_neutrino_energy",
    "totalEnergy",  
    "physical",      #Is this particle physically possible? Should always be 1
    "total_column_depth",
    
    ###################
    # I3MCTree        #
    ###################
    
    ### true HNL variables
    "HNL_true_x",
    "HNL_true_y",
    "HNL_true_z",
    "HNL_true_energy",
    "HNL_true_zenith",
    "HNL_true_azimuth",
    "HNL_true_time",  
    
    ### true primary variables
    "true_x",
    "true_y",
    "true_z",
    "true_energy",
    "true_zenith",
    "true_azimuth",
    "true_time",
    
     ### true first (DIS) cascade variables
     "casc0_true_x",
     "casc0_true_y",
     "casc0_true_z",
     "casc0_true_energy",
     "casc0_true_zenith",
     "casc0_true_azimuth",
     "casc0_true_time",  
    
     ### true second (HNL decay) cascade variables
     "casc1_true_x",
     "casc1_true_y",
     "casc1_true_z",
     "casc1_true_energy",
     "casc1_true_zenith",
     "casc1_true_azimuth",
     "casc1_true_time",    
     "nan_decay_energy",  #Was the neutrino simulated without a decay energy? Should be zero for all particles
]

# function to get HNL true variables
def store_mc_true_variables(frame):
    if not frame.Has("I3MCTree"):
        return False
    mctree = frame["I3MCTree"]
    p_true = mctree.primaries[0]

    p_daughters = mctree.get_daughters(p_true)
    assert len(p_daughters) == 2

    for p_daughter in p_daughters:
        if p_daughter.type == dataclasses.I3Particle.Hadrons:
            casc_0_true = p_daughter
        else:
            hnl_true = p_daughter

    hnl_daughters = mctree.get_daughters(hnl_true)
    assert len(hnl_daughters) > 0

    nan_energy = False

    for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):
        if not count_hnl_daughters:
            casc_1_true = hnl_daughter
        else:
            assert casc_1_true.pos == hnl_daughter.pos
            if np.isnan(hnl_daughter.energy):
                nan_energy = True
            casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy

    frame["true_x"] = dataclasses.I3Double(p_true.pos.x)
    frame["true_y"] = dataclasses.I3Double(p_true.pos.y)
    frame["true_z"] = dataclasses.I3Double(p_true.pos.z)
    frame["true_energy"] = dataclasses.I3Double(p_true.energy)
    frame["true_zenith"] = dataclasses.I3Double(p_true.dir.zenith)
    frame["true_azimuth"] = dataclasses.I3Double(p_true.dir.azimuth)
    frame["true_time"] = dataclasses.I3Double(p_true.time)

    frame["HNL_true_x"] = dataclasses.I3Double(hnl_true.pos.x)
    frame["HNL_true_y"] = dataclasses.I3Double(hnl_true.pos.y)
    frame["HNL_true_z"] = dataclasses.I3Double(hnl_true.pos.z)
    frame["HNL_true_energy"] = dataclasses.I3Double(hnl_true.energy)
    frame["HNL_true_zenith"] = dataclasses.I3Double(hnl_true.dir.zenith)
    frame["HNL_true_azimuth"] = dataclasses.I3Double(hnl_true.dir.azimuth)
    frame["HNL_true_time"] = dataclasses.I3Double(hnl_true.time)

    frame["casc0_true_x"] = dataclasses.I3Double(casc_0_true.pos.x)
    frame["casc0_true_y"] = dataclasses.I3Double(casc_0_true.pos.y)
    frame["casc0_true_z"] = dataclasses.I3Double(casc_0_true.pos.z)
    frame["casc0_true_energy"] = dataclasses.I3Double(casc_0_true.energy)
    frame["casc0_true_zenith"] = dataclasses.I3Double(casc_0_true.dir.zenith)
    frame["casc0_true_azimuth"] = dataclasses.I3Double(casc_0_true.dir.azimuth)
    frame["casc0_true_time"] = dataclasses.I3Double(casc_0_true.time)

    frame["casc1_true_x"] = dataclasses.I3Double(casc_1_true.pos.x)
    frame["casc1_true_y"] = dataclasses.I3Double(casc_1_true.pos.y)
    frame["casc1_true_z"] = dataclasses.I3Double(casc_1_true.pos.z)
    frame["casc1_true_energy"] = dataclasses.I3Double(casc_1_true.energy)
    frame["casc1_true_zenith"] = dataclasses.I3Double(casc_1_true.dir.zenith)
    frame["casc1_true_azimuth"] = dataclasses.I3Double(casc_1_true.dir.azimuth)
    frame["casc1_true_time"] = dataclasses.I3Double(casc_1_true.time)

    frame["nan_decay_energy"] = icetray.I3Bool(nan_energy)
    
    return True


# function to get EventProperties (LI)
def store_LI_event_properties(frame):
    if not frame.Has("EventProperties"):
        return False
    event_properties = frame["EventProperties"]

    frame["mHNL"] = dataclasses.I3Double(event_properties.mHNL)
    frame["distance"] = dataclasses.I3Double(event_properties.distance)
    frame["distanceMin"] = dataclasses.I3Double(event_properties.distanceMin)
    frame["distanceMax"] = dataclasses.I3Double(event_properties.distanceMax)
    frame["lifetime"] = dataclasses.I3Double(event_properties.lifetime)
    frame["totalEnergy"] = dataclasses.I3Double(event_properties.totalEnergy)
    frame["finalStateX"] = dataclasses.I3Double(event_properties.finalStateX)
    frame["finalStateY"] = dataclasses.I3Double(event_properties.finalStateY)
    frame['primary_type'] =  dataclasses.I3Double(event_properties.initialType)
    frame['final_state_particle0'] =  dataclasses.I3Double(event_properties.finalType1)
    frame['final_state_particle1'] =  dataclasses.I3Double(event_properties.finalType2)
    frame['total_column_depth'] =  dataclasses.I3Double(event_properties.totalColumnDepth)
    frame["decay_channel"] = dataclasses.I3Double(event_properties.decay_channel)
    frame["outgoing_neutrino_energy"] = dataclasses.I3Double(event_properties.outgoingNeutrinoEnergy)
    frame["physical"] = dataclasses.I3Double(event_properties.physical)

    return True

# function to skip nonphysical frames
def skip_nonphysical_frames(frame):
    '''
    Skip frames that are unphysical.
    Inputs:
        The frame
    Returns:
        True if the frame can be accessed (and therefore has contents)
        False otherwise
    '''
    if not frame.Has("EventProperties"):
        return False
    event_properties = frame["EventProperties"]
    return event_properties.physical


def main(mHNL):
    '''
    Main function: Load default inputs, set RNG and then call LeptonInjector. Skip empty frames and write hdf5.
    '''
    seed 		= int(np.random()*3e9)
    outfile 	= '$I3_SRC/LeptonInjector/resources/tests/testing_events{}.i3.zst'.format(mHNL)
    Emin 		= float(2)
    Emax 		= float(10000)
    index 		= float(2)
    nEvents 	= 5000
    Zmin 		= float(80)
    Zmax 		= float(180)
    radius 		= float(600)
    length		= float(600)
    HNL_mass    = mHNL
    write_hdf5  = True

    t0 = time.perf_counter() # time.clock()

    print('Starting ...')

    # The .lic file contains the simulation data used for weighting.
    generation_data_file_path = outfile.replace('.i3.zst', '.lic')

    xs_location = os.path.expandvars('$I3_SRC/LeptonInjector_HNL/resources/cross_sections/M_{:04.0f}MeV'.format(HNL_mass*1e3))

    tray = I3Tray()

    tray.context["I3FileStager"] = dataio.get_stagers()

    tray.AddService("I3GSLRandomServiceFactory")(("Seed",seed))
    tray.AddService("I3EarthModelServiceFactory","Earth")
    tray.AddModule("I3InfiniteSource")(("Stream",icetray.I3Frame.DAQ))

    # Set the cross section. Generate half of the events as neutrinos and the other half as antineutrinos.
    tray.AddModule("MultiLeptonInjector")(
            ("EarthModel","Earth"),
            ("Generators",[
                LeptonInjector.injector(
                    NEvents    = int(nEvents/2),
                    FinalType1 = dataclasses.I3Particle.ParticleType.HNL,
                    FinalType2 = dataclasses.I3Particle.ParticleType.Hadrons,
                    DoublyDifferentialCrossSectionFile 	= os.path.join(xs_location, "dsdxdy-nutau-N-nc-GRV98lo_patched_central.fits"),
                    TotalCrossSectionFile 				= os.path.join(xs_location, "sigma-nutau-N-nc-GRV98lo_patched_central.fits"),
                    Ranged     = False)
                ,
                LeptonInjector.injector(
                    NEvents    = int(nEvents/2),
                    FinalType1 = dataclasses.I3Particle.ParticleType.HNLBar,
                    FinalType2 = dataclasses.I3Particle.ParticleType.Hadrons,
                    DoublyDifferentialCrossSectionFile  = os.path.join(xs_location, "dsdxdy-nutaubar-N-nc-GRV98lo_patched_central.fits"),
                    TotalCrossSectionFile               = os.path.join(xs_location, "sigma-nutaubar-N-nc-GRV98lo_patched_central.fits"),
                    Ranged     = False)
                ]),
            ("MinimumEnergy", 	Emin * I3Units.GeV),
            ("MaximumEnergy", 	Emax * I3Units.GeV),
            ("MinimumZenith", 	Zmin * I3Units.deg),
            ("MaximumZenith", 	Zmax * I3Units.deg),
            ("PowerlawIndex", 	index),
            # ("InjectionRadius",	radius * I3Units.meter),  # these are only for ranged mode
            # ("EndcapLength",	length * I3Units.meter),  # these are only for ranged mode
            ("CylinderRadius",	radius * I3Units.meter),
            ("CylinderHeight",	length * I3Units.meter),
            ("HNL_mass",        HNL_mass * I3Units.GeV),
            ("CylinderCenter", -300. * I3Units.meter),  # just fix the z center of DC here. Can remove for simulation ouside deepcore.
            )

	# get configuration output ready for staging (slightly hacky)
    if 'I3FileStager' in tray.context:
        stager = tray.context['I3FileStager']
        generation_data_file_path = stager.GetWriteablePath(generation_data_file_path)
    
    tray.AddModule("InjectionConfigSerializer", OutputPath = str(generation_data_file_path))

    # skip empty frames
    tray.AddModule(
        skip_nonphysical_frames,
        streams=[icetray.I3Frame.DAQ]
    )

    tray.Add(
        "I3Writer",
        filename = outfile,
        Streams = [icetray.I3Frame.DAQ,],
    )

    if write_hdf5:

        outfile_hdf5 = outfile.replace('.i3.zst',".hdf5")

        # get hnl true energy
        tray.AddModule(
            store_mc_true_variables,
            "store_mc_true_variables",
            Streams=[icetray.I3Frame.DAQ],
        )

        # get event properties (LeptonInjector)
        tray.AddModule(
            store_LI_event_properties,
            "store_LI_event_properties",
            Streams=[icetray.I3Frame.DAQ],
        )
        
        tray.AddSegment(
            I3SimHDFWriter,
            output = outfile_hdf5,
            keys = HDF5_KEYS,
        )

    tray.Execute()
    tray.Finish()


    t1 = time.perf_counter() # time.clock()

    print("Time it took: {}s".format(t1-t0))
    print('done ...')
    
    return True
    
def test_helper():
    passed = main(0.1) and main(0.3) and main(0.6) and main(1.0)
    
    return (passed)
    
