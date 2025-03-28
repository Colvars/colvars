#!/usr/bin/env python

# Script to process trajectories from Colvars-based bias-exchange simulations

import os
import sys
import glob

import numpy as np

# Add "colvartools" folder to import path
sys.path += [os.path.dirname(os.path.dirname(os.path.realpath(__file__)))]
from plot_colvars_traj import Colvars_traj, Colvar_traj



def compute_traj_masks(traj):
    """
    Compute arrays with values 1 or 0 marking frames where each bias is active
    """

    traj_mask = { }
    num_frames = traj[traj.variables[0]].steps.size

    biases_active = [v for v in traj.variables if v[-7:] == '_active']

    if len(biases_active) == 0:

        print("Warning: Cannot find any columns labeled \"*_active\" in the Colvars trajectory.")
        print("    For a BE run, \"outputFeatures active\" must be set for each bias;")
        print("    without it, all trajectories are assumed to be *unbiased*.")

        traj_mask['neutral'] = np.ones(shape=(num_frames), dtype=np.int32)

    else:

        for col in biases_active:
            bias = col[:-7]
            # Check that all bias trajs have consistent numbers of frames
            assert num_frames == traj[col].values.size
            traj_mask[bias] = np.zeros(shape=(num_frames), dtype=np.int32)
            traj_mask[bias][traj[col].values > 0] = 1

    for bias in traj_mask.keys():
        prefix = "Bias \"" + bias + "\" has"
        print(prefix, np.sum(traj_mask[bias]), "frames")

    return traj_mask


def compute_histograms(trajs, bins, range,
                       biases=[], variables=[],
                       output_prefix="histogram"):
    if len(biases) == 0:
        biases = trajs.keys()
    for bias in biases:
        if len(variables) == 0:
            variables = trajs[bias].keys()
        for variable in variables:
            hist, edges = np.histogram(
                trajs[bias][variable].values, bins=bins, range=range,
                density=True)
            np.savetxt("%s.%s-%s.dat" % (output_prefix, bias, variable),
                       np.c_[(edges[1:]+edges[:-1])/2.0, hist])


def write_constant_bias_trajs(trajs, biases=[], variables=[],
                              output_prefix="colvars_traj"):
    if len(biases) == 0:
        biases = trajs.keys()
    for bias in biases:
        if len(variables) == 0:
            variables = trajs[bias].keys()
        pluto=trajs[bias][variables[0]].steps
        for i in range (0,len(variables)):
           pluto=np.c_[pluto,trajs[bias][variables[i]].values]
        pluto=pluto[pluto[:, 0].argsort()] 
        np.savetxt("%s.%s.dat" % (output_prefix, bias),pluto,
                   header="step" " " + " ".join(str(variable) for variable in variables))
           
#        with open("%s.%s.dat" % (output_prefix, bias), 'w') as f:
#            f.write(" # ")
#            for variable in variables:
#               f.write(" %s " % (variable))  
#            f.write(" \n")
#            for step in range (0,trajs[bias][variables[0]].steps.size):
#               f.write(" %s " % trajs[bias][variables[0]].steps[step])
#               for variable in variables:
#                  f.write(" %s " % (trajs[bias][variable].values[step]))
#               f.write("\n")  

def extract_constant_bias_trajs(histogram_bins, histogram_range,
                                pattern="replica-%03d/*.colvars.traj",
                                variables=[], biases=[],
                                first=0, last=None, skip=1, do_hist=False):

    constant_bias_trajs = { }

    replica_index = 0

    while True:

        # Loop over replicas
        traj_files = sorted(glob.glob(pattern % replica_index))
        if len(traj_files) == 0:
            break
        traj = Colvars_traj(traj_files, first=first, last=last, every=skip)

        print()

        # Defaults: metadynamics biases, one for each CV
        if len(variables) == 0:
            variables = [variable for variable in traj.variables
                         if variable[0:3] != 'mtd_']
        if len(biases) == 0:
            biases = [('mtd_' + v) for v in variables
                      if ('E_mtd_' + v) in traj.variables]

        print("Analyzing variables:", variables)
        print("Analyzing biases:", biases)

        print("Replica", replica_index, "has",
              traj[traj.variables[0]].steps.size, "frames")

        traj_mask = compute_traj_masks(traj)

        for bias in biases:

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

    if do_hist:
      compute_histograms(
          trajs=constant_bias_trajs,
          biases=biases, variables=variables,
          bins=histogram_bins, range=histogram_range,
          output_prefix="histogram"
      )

    write_constant_bias_trajs(
        trajs=constant_bias_trajs, 
        biases=biases, 
        variables=variables,
        output_prefix="colvars_traj"
    )

if __name__ == '__main__':

    import argparse

    parser = \
        argparse.ArgumentParser(description='Script to process trajectories from Colvars-based bias-exchange simulations.', \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--variables',
                        type=str, nargs='*',
                        help="Which collective variables to analyze "
                        "(default: all that are found in colvars.traj files",
                        default=[])

    parser.add_argument('--biases',
                        type=str, nargs='*',
                        help="Which biases to analyze"
                        "(default: all biases named \"mtd_<variable>\")",
                        default=[])

    parser.add_argument('--pattern',
                        type=str,
                        help="Filename pattern to use for each replica; "
                        "must contain a format specifier within which to place "
                        "the replica index",
                        default="replica-%03d/*.colvars.traj")

    parser.add_argument('--first',
                        dest='first',
                        type=np.int64,
                        help='First frame to read',
                        default=-1)

    parser.add_argument('--last',
                        dest='last',
                        type=np.int64,
                        help='Last frame to read',
                        default=-1)

    parser.add_argument('--skip',
                        dest='skip',
                        type=np.int64,
                        help='Read every these many frames',
                        default=1)

    parser.add_argument('--histogram-bins',
                        type=int,
                        help='Number of bins in the histograms',
                        default=120)

    parser.add_argument('--histogram-range',
                        type=float, nargs=2,
                        metavar=('xi_min', 'xi_max'),
                        help='Range of the histograms',
                        default=[0.0, 12.0])

    parser.add_argument("--do-histograms", 
                        help="calculate histograms of each variable", 
                        default=False, 
                        dest='do_hist', 
                        action='store_true')

    kwargs = vars(parser.parse_args())

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    extract_constant_bias_trajs(**kwargs)
