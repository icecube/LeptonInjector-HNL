import numpy as np
import os

import LeptonWeighter as LW

from icecube.LeptonInjector_HNL import HNLDecayEnergy

from icecube.photospline import spglam as glam
from icecube.photospline import splinefitstable

import icecube.icetray.i3logging as logging

from icecube import dataclasses

#from icecube.oscNext.tools.misc import UpdateFrameObject
from icecube.LeptonInjector_HNL.oscnext_tools_misc import UpdateFrameObject
FullWidth_Array = np.vectorize(lambda x: HNLDecayEnergy.FullWidth(x))

#
# Globals
#

WEIGHT_DICT_KEY = "I3MCWeightDict"
WEIGHT_KEY_ONEWEIGHT = "OneWeight"
WEIGHT_KEY_LEPTONINJECTOR = "LeptonInjectorWeight"
WEIGHT_KEY_LIFETIME = "LifetimeWeight"
WEIGHT_KEY_REFERENCE = "ReferenceWeight"

LIGHTSPEED = 299792458.0  # in m/s
REDUCEDPLANCK = 6.582119569e-25  # in GeV*s


##### New generation level weighting functions #####

def load_generation_weighting_files(
    flux_file='{}/LeptonInjector/resources/flux_files/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000_2006_2014.hdf5'.format(os.environ['I3_SRC']),
    xs_location='{}/LeptonInjector/resources/cross_sections/M_0100MeV/'.format(os.environ['I3_SRC']),
):
    """
    Load files and splines needed for HNL weighting (just once).

    Parameters
    ----------
    flux_file : string
        Location of (hdf5) nusquids flux(+oscillation) file.
    xs_location : string
        Directory of cross section files.

    Returns
    -------
    nusquids_flux : LW.nuSQUIDSAtmFlux
        Loaded flux(+oscillation) file.
    leptonweighter_xsec : LW.CrossSectionFromSpline
        Loaded target cross section (differential NuTau/NuTauBar for NC/CC).
    """

    # load flux file
    nusquids_flux = LW.nuSQUIDSAtmFlux(flux_file)

    # load cross sections
    dsdxdy_nu_CC_nominal    = os.path.join(xs_location, 'dsdxdy-nutau-N-cc-GRV98lo_patched_central.fits')
    dsdxdy_nu_NC_nominal    = os.path.join(xs_location, 'dsdxdy-nutau-N-nc-GRV98lo_patched_central.fits')
    dsdxdy_nubar_CC_nominal = os.path.join(xs_location, 'dsdxdy-nutaubar-N-cc-GRV98lo_patched_central.fits')
    dsdxdy_nubar_NC_nominal = os.path.join(xs_location, 'dsdxdy-nutaubar-N-nc-GRV98lo_patched_central.fits')
    leptonweighter_xsec     = LW.CrossSectionFromSpline(dsdxdy_nu_CC_nominal, dsdxdy_nubar_CC_nominal, dsdxdy_nu_NC_nominal, dsdxdy_nubar_NC_nominal)

    return nusquids_flux, leptonweighter_xsec,


def weight_hnl_generation(
    frame,
    nusquids_flux,
    leptonweighter_xsec,
    leptoninjector_config,
    verbose=False,
):
    """
    Calculate the initial HNL weight (OneWeight and weight including flux+oscillation) and store it to the frame.

    Parameters
    ----------
    frame : icetray.I3Frame
        Current frame.
    nusquids_flux : LW.nuSQUIDSAtmFlux
        Loaded flux(+oscillation) file.
    leptonweighter_xsec : LW.CrossSectionFromSpline
        Loaded target cross section (differential NuTau/NuTauBar for NC/CC).
    leptoninjector_config : string
        Path to LeptonInjector configuration file (.lic).
    """

    # Get the weight dict
    with UpdateFrameObject(frame, WEIGHT_DICT_KEY) as weight_dict:

        # get primary
        if 'I3MCTree' in frame:
            primary = frame['I3MCTree'].get_primaries()[0]
        else:
            logging.log_error('I3MCTree missing in the frame.')
            return False
        
        # get event properties
        EventProperties = frame['EventProperties']

        # set up this particular event
        LWevent = LW.Event()
        LWevent.primary_type           = LW.ParticleType(EventProperties.initialType)
        LWevent.final_state_particle_0 = LW.ParticleType(EventProperties.finalType1)
        LWevent.final_state_particle_1 = LW.ParticleType(EventProperties.finalType2)
        LWevent.zenith                 = EventProperties.zenith
        LWevent.energy                 = EventProperties.totalEnergy
        LWevent.azimuth                = EventProperties.azimuth
        LWevent.interaction_x          = EventProperties.finalStateX
        LWevent.interaction_y          = EventProperties.finalStateY
        LWevent.total_column_depth     = EventProperties.totalColumnDepth
        LWevent.x 				       = primary.pos.x
        LWevent.y 				       = primary.pos.y
        LWevent.z 				       = primary.pos.z - EventProperties.cylinderCenter  # shift to account for shifted injection cylinder

        # load leptoninjector configuration for this event
        simulation_generation = LW.MakeGeneratorsFromLICFile(leptoninjector_config)

        # set up nominal weighter
        weighter_nominal = LW.Weighter(
            nusquids_flux,
            leptonweighter_xsec,
            simulation_generation,
            )

        if verbose:
            logging.log_info('Neutrino energy: {}'.format(LWevent.energy))
            logging.log_info('Oneweight: {}'.format(weighter_nominal.get_oneweight(LWevent)))
            logging.log_info('LI weight: {}'.format(weighter_nominal(LWevent)))

        weight_dict[WEIGHT_KEY_ONEWEIGHT] = weighter_nominal.get_oneweight(LWevent)
        weight_dict[WEIGHT_KEY_LEPTONINJECTOR] = weighter_nominal(LWevent)

    return True


def weight_hnl_lifetime_framewise(
    frame,
    U_tau4_sq=1e-03,
    c=LIGHTSPEED,
    hbar=REDUCEDPLANCK,
    add_reference_weight=True,
):
    """
    Function to calculate the lifetime weight per frame (assuming 1/L length sampling) and also add total reference weight for this mixing if wanted.

    Parameters
    ----------
    U_tau4_sq : float
        Chosen mixing strength.
    """

    # Get the weight dict
    with UpdateFrameObject(frame, WEIGHT_DICT_KEY) as weight_dict:

        ##### get needed inputs #####
        # get event properties
        EventProperties = frame['EventProperties']
        mass = EventProperties.mHNL  # HNL mass in GeV

        # get hnl energy
        hnl_true = (frame["I3MCTree"].get_daughters(frame["I3MCTree"].primaries[0]))[0]
        assert hnl_true.type == dataclasses.I3Particle.HNL or dataclasses.I3Particle.HNLBar
        energy = hnl_true.energy  # HNL kinetic energy in GeV

        # get event properties
        tau = EventProperties.lifetime  # HNL lifetime in rest frame in ns
        distance_min = EventProperties.distanceMin  # generation decay length minimum in m
        distance_max = EventProperties.distanceMax  # generation decay length maximum in m
        ##### end #####

        gamma = np.sqrt(energy**2 + mass**2) / mass  # Etot/E0
        speed = c * np.sqrt(1 - np.power(1.0 / gamma, 2))  # c * sqrt(1-1/gamma^2)

        tau_min = distance_min / (gamma * speed)
        tau_max = distance_max / (gamma * speed)

        tau = 1e-09 * tau  # convert to seconds
        tau_proper = hbar / (FullWidth_Array(mass) * U_tau4_sq)

        pdf_inverse = (1.0 / (np.log(tau_max) - np.log(tau_min))) * (
            1.0 / tau
        )  # for 1/L sampling of decay length

        pdf_exp1 = 1.0 / tau_proper
        pdf_exp2 = np.exp(-tau / tau_proper)

        pdf_exp = pdf_exp1 * pdf_exp2

        weight_lifetime = pdf_exp / pdf_inverse

        weight_dict['{}_{:.0e}'.format(WEIGHT_KEY_LIFETIME, U_tau4_sq)] = weight_lifetime

        if add_reference_weight:
            if WEIGHT_KEY_LEPTONINJECTOR in weight_dict:
                LI_weight = weight_dict[WEIGHT_KEY_LEPTONINJECTOR]
                weight_dict['{}_{:.0e}'.format(WEIGHT_KEY_REFERENCE, U_tau4_sq)] = U_tau4_sq * weight_lifetime * LI_weight

    return True
