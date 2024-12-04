from unittest import TestCase
import unittest
import os
from os.path import expandvars
import glob
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from icecube.icetray import I3Tray
from icecube import icetray, dataio, dataclasses, LeptonInjector_HNL
from icecube import tableio, hdfwriter
from I3Tray import I3Units

from icecube.LeptonInjector_HNL import hdf5_to_feather, generate_test_events_public

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

ref_base_dir = '/data/user/jbook/I3_HNL_Search/validation/testfiles_Feb5_24'
test_base_dir = os.getcwd()#os.path.expandvars('$I3_SRC/LeptonInjector/resources/tests/')

def eq_w_t(val1, val2, t):
    '''
    Tests whether 2 values are equal within some tolerance t
    Inputs:
    val1, val2: the value to be compared
    t: the absolute tolerance
    '''
    
    if (val1<1 or val2<1):
        t = 1
    return np.abs(val1-val2)<=t
def flat_w_t(data):
    '''
    Tests whether a distribution is approximately flat
    Inputs:
    data: the data you hope is flat, array-like
    '''
    nbins = 50
    n, b = np.histogram(data, bins = nbins)
    mean = np.mean(n)
    std = np.std(n)
    
    within_one = ((mean - std < n) & (n < mean + std)).sum()
    within_two = ((mean - std*2 < n) & (n < mean + std*2)).sum()
    within_three = ((mean - std*3 < n) & (n < mean + std*3)).sum()
    
    first_dev = eq_w_t(within_one/len(n), 0.68, np.sqrt(within_one)/len(n)) or within_one/len(n) > 0.68
    second_dev = eq_w_t(within_two/len(n), 0.95, np.sqrt(within_two)/len(n)) or within_two/len(n)>0.95
    third_dev = eq_w_t(within_three/len(n), 0.997, np.sqrt(within_three)/len(n)) or within_three/len(n)>0.997
    
    result = [first_dev, second_dev, third_dev]

    
    return np.all(result)


################################################################
# Set up some global variables used by multiple tests 

# When simulating sets at multiple masses, list the file prefixes for each mass here
levels = ['Gen_01', 'Gen_03', 'Gen_06', 'Gen_10']
test_levels = ['0.1', '0.3', '0.6', '1.0']

# dicts for infilepaths
ref_infiles = OrderedDict(zip(levels, [list() for level in levels]))
test_infiles = OrderedDict(zip(test_levels, [list() for level in test_levels]))

#     # get the infilepaths
#     for key, item in ref_infiles.items():
#         print(key)  
#         item.extend(
#             glob.glob(os.path.join(ref_base_dir+ '/{}*.hdf5'.format(key))) #for reference files
#         )
print('########################################################')
print("Events generated during the test will be stored in :")

for key, item in test_infiles.items(): 
    item.extend(
        glob.glob(os.path.join(test_base_dir+ '/*{}.hdf5'.format(key))) #for test files
    )
    # List the infile paths
    print(item)
print('########################################################')

# Read in the data
ref_data = OrderedDict()
test_data = OrderedDict()

# Number of events to generate at each mass point for testing
n_requested = 5000

class TestGeneration(unittest.TestCase):
#################################################################
## Generate test files and ensure the generator runs correctly ##
#################################################################

################################################################
## Collect the events we just generated (test_data)           ##
## and a sample set we know is correctly generated (ref_data) ##
## This will be used for all further tests                    ##
################################################################
    @classmethod
    def setUpClass(self):
        '''
        Tests that the provided generation script produces events.
        The events produced (when the test passes) are then used for further checks.
        The events are generates with default parameters, a seed from pn.random(), and 5000 events are requested.
        '''
        
        generated = generate_test_events_public.test_helper(n_requested)
        
        global test_infiles 
        global test_data
        
        for key, item in test_infiles.items():   
            test_data[key] = hdf5_to_feather(item, keys=keys_to_extract, outfilepath=test_base_dir+"/test_files.feather")
        return generated

    ###################################
    ## Make this into the first test ##
    ###################################

    def test_num_events(self):
        '''
        tests that the number of events requested is the number generated, within expected tolerance
        '''
        global test_data
        for mass, data in test_data.items():
            n_events = len(data['totalEnergy'])
            assert n_events<=n_requested
            #It doesn't matter which key we use - they all have the same length
            assert n_events <= n_requested*1.05 and n_events >= n_requested*0.95
    
    ###################################
    ## Test energy conservation      ##
    ###################################

    def test_energy_conservation(self):
    
        for m, data in test_data.items():
            mass = float(m)
            assert np.all(data['true_energy']==data['totalEnergy'])
            assert np.all(mass == data['mHNL'][0])

            HNL_total_E = np.sqrt(data['HNL_true_energy']**2 + mass**2)
            
            assert np.all(np.absolute(data['true_energy'] - data['casc0_true_energy'] - HNL_total_E)<10**(-6))
            
            #Test gaussian within variance
            #abs(hnl-casc1)/hnl
            energy_mask = data['HNL_true_energy']>2 #GeV
            
            casc1_residual = HNL_total_E[energy_mask] - data['casc1_true_energy'][energy_mask]                
            fractional_error = casc1_residual/HNL_total_E[energy_mask]
            
            frac_err_mask = fractional_error>0.07 #Look for events with fractional error > 7%
            percent_events_failing = len(fractional_error[frac_err_mask])/len(fractional_error)
            
            if (percent_events_failing > 0.05):
                print ('Failure reason: greater than 5% of events have 7% fractional uncertainty')
                assert False
            if (percent_events_failing > 0.03):
                print ('Warning: greater than 3% of events have 7% fractional uncertainty')
            
            else:
                if (percent_events_failing < 0.03):
                    assert True
            
            
    ###################################
    ## Test branching ratio          ##
    ###################################

    def test_branching_ratio(self):
        from icecube.LeptonInjector_HNL import HNLDecayEnergy

        #We define these histograms to do the bin counting easily
        bins = np.linspace (0.5, 15.5, 15)
        br = np.empty(shape = [4, 14])
        err = np.empty(shape = [4, 14])

        for i, m in enumerate(test_levels):
            n2, b = np.histogram(test_data[m]['decay_channel'], bins = bins)
            br[i]=n2/len(test_data[m]['decay_channel'])
            err[i] = np.sqrt(n2)/len(test_data[m]['decay_channel'])

        branches = br.transpose()
        errors = err.transpose()

        for i, x in enumerate(test_levels):
            # x = HNL mass
            x = float(x)
            full_width = HNLDecayEnergy.FullWidth(x)
            assert eq_w_t(HNLDecayEnergy.gamma_nu_nu_nu_overload(x)/full_width, branches[0][i], errors[0, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_e_e(x)/full_width, branches[1][i], errors[1, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_p0_nu(x)/full_width, branches[2][i], errors[2, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_nu_mu_mu(x)/full_width, branches[3][i], errors[3, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_eta_nu(x)/full_width, branches[4][i], errors[4, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_rho0_nu(x)/full_width, branches[5][i], errors[5, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_omega_nu(x)/full_width, branches[6][i], errors[6, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_etaprime_nu(x)/full_width, branches[7][i], errors[7, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_phi_nu(x)/full_width, branches[8][i], errors[8, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_nue_e_tau(x)/full_width, branches[9][i], errors[9, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_numu_mu_tau(x)/full_width, branches[10][i], errors[10, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_tau_pi(x)/full_width, branches[11][i], errors[11, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_tau_K(x)/full_width, branches[12][i], errors[12, i]*3)
            assert eq_w_t(HNLDecayEnergy.gamma_tau_rho(x)/full_width, branches[13][i], errors[13, i]*3)

    ###################################
    ## Test simulation distributions ## 
    ###################################


    def test_flat_cos_theta(self):
        for m in test_levels:
            assert flat_w_t(np.cos(test_data[m]['true_zenith']))

    def test_power_law(self):
        for m in test_levels:
            nEvents = len(test_data[m]['true_energy'])
            x = np.random.rand(nEvents)
            i = 2
            minE = 2
            maxE = 1000
            energyP = (1-x)*pow(minE,1-i) + x*pow(maxE,1-i)
            sampled_energy = (pow(energyP,1/(1-i)))                

            n, b = np.histogram(sampled_energy, 50)
            n1, b1 = np.histogram(test_data[m]['true_energy'], 50)
            diff_from_sampling_dist = n - n1
            assert flat_w_t(diff_from_sampling_dist)

    def test_physical(self):
        for m in test_levels:
            assert np.all(test_data[m]['physical']==True)


    if __name__ == '__main__':
        unittest.main()

















