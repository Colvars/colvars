#!/usr/bin/env python

# Script to process trajectories from Colvars-based bias-exchange simulations

import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt

# Add "colvartools" folder to import path
sys.path += [os.path.dirname(os.path.dirname(os.path.realpath(__file__)))]
from plot_colvars_traj import Colvars_traj, Colvar_traj



def compute_traj_masks(traj, biases):
    """
    Compute arrays with values 1 or 0 marking frames where each bias is active
    """

    traj_mask = { }
    num_frames = traj[traj.variables[0]].steps.size
    traj_mask['neutral'] = np.zeros(shape=(num_frames), dtype=np.int32)

    for bias in biases:
        if bias+'_active' in traj:
            # Check that all bias trajs have consistent numbers of frames
            assert num_frames == traj[bias+'_active'].values.size
            traj_mask[bias] = np.zeros(shape=(num_frames), dtype=np.int32)
            traj_mask[bias][traj[bias+'_active'].values > 0] = 1
            traj_mask['neutral'] = np.add(traj_mask['neutral'], traj_mask[bias])
        else:
            print("Warning: Cannot find column \""+bias+
                  "_active\" in Colvars trajectory.")
            print("    For a BE run, \"outputFeatures active\" must be set for each bias, but")
            print("    without it trajectories are assumed to be *unbiased*.")

    traj_mask['neutral'] = 1 - traj_mask['neutral']
    for bias in traj_mask.keys():
        print(bias, "has", np.sum(traj_mask[bias]), "frames")

    return traj_mask


def compute_histograms(trajs, bins=120, range=[0.0, 12.0],
                       output_prefix="histogram"):
    for bias in trajs.keys():
        for variable in trajs[bias].keys():
            hist, edges = np.histogram(
                trajs[bias][variable].values, bins=bins, range=[0.0, 12.0],
                density=True)
            np.savetxt("%s.%s-%s.dat" % (output_prefix, bias, variable),
                       np.c_[(edges[1:]+edges[:-1])/2.0, hist])


def write_constant_bias_trajs(trajs, output_prefix="colvars_traj"):
    for bias in trajs.keys():
        for variable in trajs[bias].keys():
            np.savetxt("%s.%s-%s.dat" % (output_prefix, bias, variable),
                       np.c_[trajs[bias][variable].steps,
                             trajs[bias][variable].values])


def extract_constant_bias_trajs(pattern="replica-%03d/*.colvars.traj",
                                variables=[], biases=[],
                                first=0, last=None, skip=10):

    constant_bias_trajs = { }

    replica_index = 0

    while True:

        print()

        traj_files = sorted(glob.glob(pattern % replica_index))
        if len(traj_files) == 0:
            break
        traj = Colvars_traj(traj_files, first=first, last=last, every=skip)

        # Defaults: metadynamics biases, one for each CV
        if len(variables) == 0:
            variables = [variable for variable in traj.variables
                         if variable[0:3] != 'mtd_']
        if len(biases) == 0:
            biases = [('mtd_' + v) for v in variables
                      if ('E_mtd_' + v) in traj.variables]

        traj_mask = compute_traj_masks(traj, biases)

        print("Replica", replica_index, "has", traj_mask['neutral'].size, "frames")

        for bias in ['neutral']:

            if not bias in constant_bias_trajs:
                constant_bias_trajs[bias] = {}

            for variable in traj.variables:

                traj_steps = traj[variable].steps[traj_mask[bias] == 1]
                traj_variable = traj[variable].values[traj_mask[bias] == 1]

                if not variable in constant_bias_trajs[bias]:
                    constant_bias_trajs[bias][variable] = Colvar_traj(variable)

                # TODO concatenate directly inside Colvar_traj class, without
                # accessing private members
                constant_bias_trajs[bias][variable]._step = np.concatenate([
                    constant_bias_trajs[bias][variable].steps,
                    traj_steps])
                constant_bias_trajs[bias][variable]._colvar = np.concatenate([
                    constant_bias_trajs[bias][variable].values,
                    traj_variable])

        replica_index += 1

    compute_histograms(constant_bias_trajs)


extract_constant_bias_trajs()
