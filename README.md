# LeptonInjector-HNL
Open-source event generator for HNL production from a tau neutrino flux, compatible with IceCube's IceTray analysis software.

# Installation Instructions
This version of LeptonInjector-HNL is compatible with <a href="https://github.com/icecube/icetray-public">icetray</a>, IceCube's analysis framework, for serialization. Download the reposity, adding LeptonInjector-HNL to the source directory. Then, replace I3Particle.cxx(.h) in dataclasses/pricate(public)/dataclasses/physics with the files of the same name in dataclasses_files. This will allow support for HNL particles within IceTray. 

The three-body kinematic distributions of HNL decays were calculated using madgraph, and must be downloaded separately from the Harvard Dataverse <a href="https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FDCUQG">here</a>. A bash script for this purpose is provided in the madgraph folder.

Finally, follow installation instructions for <a href="https://github.com/icecube/icetray-public">icetray</a>.

For help with installation, please contact Julia Book Motzkin (jbook@g.harvard.edu).


# Contact
LeptonInjector-HNL was developed for HNL searches in IceCube, by Leander Fischer (leander.fischer@desy.de) and Julia Book Motzkin (jbook@g.harvard.edu) on behalf of the IceCube Collaboration.
