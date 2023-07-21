#!/usr/bin/env python
# coding: utf-8

import sys
import timeit
import torch
import numpy as np
import os
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals 
import pandas as pd

from molann.feature import FeatureFileReader
import molann.ann as ann 

use_gpu = False
traj_path = '/home/numerik/bzfzhang/tmp/test_pool/23-03-2023-md-traj-plot/alanine-dipeptide/water-1500ns-1/'
feature_file = 'feature_dihedrals.txt'
input_selector = 'type C or type O or type N'

traj_filename = os.path.join(traj_path, 'md_center.xtc')
top_filename = os.path.join(traj_path, './top.gro')

print (f'===Computing Devices===')
# CUDA support
if use_gpu:
    print ('Active CUDA Device: GPU', torch.cuda.current_device())
    print ('Available devices: ', torch.cuda.device_count())
    print ('CUDA name: ', torch.cuda.get_device_name(0))

print (f'=======================\n')

# load trajectory data 
universe = mda.Universe(top_filename, traj_filename)

atoms_info = pd.DataFrame(
    np.array([universe.atoms.ids, universe.atoms.names,
        universe.atoms.types, universe.atoms.masses,
        universe.atoms.resids, universe.atoms.resnames]).T, 
    columns=['id', 'name', 'type', 'mass', 'resid', 'resname']
    )

print ('==========System Info=================\n', atoms_info)

print ('\nSummary:\n', atoms_info['type'].value_counts().rename_axis('type').reset_index(name='counts'))

# print information of trajectory
print ('{} atoms, {} residues.'.format(universe.atoms.n_atoms, universe.residues.n_residues) )
print ('==========End of System Info==========\n')

print ('==============Features===================\n')
print ('Features file: {}'.format(feature_file)) 

# read features from file to define preprocessing
feature_reader = FeatureFileReader(feature_file, 'Preprocessing', universe)
feature_list = feature_reader.read()

# define the map from positions to features 
input_ag = universe.select_atoms(input_selector)
feature_mapper = ann.FeatureLayer(feature_list, input_ag, use_angle_value=True)

print ('\nFeatures in preprocessing layer:')
# display information of features used 
print (feature_mapper.get_feature_info())

print ('==============End of Features===================\n')

pp_layer = ann.PreprocessingANN(None, feature_mapper)

scripted_cv_filename = f'./scripted_dihedral_angles_cpu.pt'
torch.jit.script(pp_layer).save(scripted_cv_filename)
print (f'  script (CPU) model for CVs saved at:\n\t{scripted_cv_filename}\n', flush=True)

