BATCH -c 1 # number of cores
#SBATCH -t 0-00:30  # Runtime in D-HH:MM
#SBATCH -p shared   # Partition to run on
#SBATCH --mem=10G   # RAM
#SBATCH -o download_output.out # File for STDOUT to write to
#SBATCH -e download_error.out  # File for STDERR to write to

# Give the commands

# Dirac, 0.1 GeV
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.1GeV_vee.dat https://dataverse.harvard.edu/api/access/datafile/9264752
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.1GeV_vtveve.dat https://dataverse.harvard.edu/api/access/datafile/9264754
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.1GeV_vtvmvm.dat  https://dataverse.harvard.edu/api/access/datafile/9264761
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.1GeV_vtvtvt.dat  https://dataverse.harvard.edu/api/access/datafile/9264762

# Dirac, 0.3 GeV
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.3GeV_vee.dat  https://dataverse.harvard.edu/api/access/datafile/9264758
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.3GeV_vmm.dat  https://dataverse.harvard.edu/api/access/datafile/9264764
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.3GeV_vtveve.dat  https://dataverse.harvard.edu/api/access/datafile/9264767
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.3GeV_vtvmvm.dat  https://dataverse.harvard.edu/api/access/datafile/9264768
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.3GeV_vtvtvt.dat  https://dataverse.harvard.edu/api/access/datafile/9264759

# Dirac, 0.6 GeV
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.6GeV_vee.dat  https://dataverse.harvard.edu/api/access/datafile/9264751
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.6GeV_vmm.dat  https://dataverse.harvard.edu/api/access/datafile/9264756
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.6GeV_vtveve.dat  https://dataverse.harvard.edu/api/access/datafile/9264753
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.6GeV_vtvmvm.dat  https://dataverse.harvard.edu/api/access/datafile/9264749
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_0.6GeV_vtvtvt.dat  https://dataverse.harvard.edu/api/access/datafile/9264760

# Dirac, 1.0 GeV
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_1GeV_vee.dat  https://dataverse.harvard.edu/api/access/datafile/9264757
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_1GeV_vmm.dat  https://dataverse.harvard.edu/api/access/datafile/9264765
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_1GeV_vtveve.dat  https://dataverse.harvard.edu/api/access/datafile/9264766
wget -c -P DatosDirac/ --output-document=/HNL_Decay_Dirac_1GeV_vtvmvm.dat  https://dataverse.harvard.edu/api/access/datafile/9264755
wget -c -P DatosDirac/ --output-document=HNL_Decay_Dirac_1GeV_vtvtvt.dat  https://dataverse.harvard.edu/api/access/datafile/9264769

# Majorana, 0.1 GeV
wget -c -P DatosMajorana/ https://dataverse.harvard.edu/api/access/datafile/9264787
wget -c -P DatosMajorana/ https://dataverse.harvard.edu/api/access/datafile/926478
3
wget -c -P DatosMajorana/ https://dataverse.harvard.edu/api/access/datafile/926478
5
wget -c -P DatosMajorana/ https://dataverse.harvard.edu/api/access/datafile/9264780
