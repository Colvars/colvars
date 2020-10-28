#!/usr/bin/env python

# Invoke this script with --help for documentation.

# Used in this tutorial: https://colvars.github.io/multi-map/multi-map.html

# Download link: https://github.com/Colvars/colvars/blob/master/colvartools/gen_multimap.py?raw=true
# Last updated: 2020-10-28
# Contact: giacomo.fiorin@gmail.com


import os, sys, tempfile

import numpy as np

import scipy.ndimage as spnd
from scipy.interpolate import InterpolatedUnivariateSpline


try:
    import gridData
except ImportError:
    print("Module gridData not found: "
          "it is usually bundled with MDAnalysis.")
    if (__name__ == '__main__'):
        sys.exit(1)


def format_array(a, fmt="%.3f"):
    if type(a) == np.ndarray:
        return np.array2string(a,
                               formatter={'float_kind': lambda x: fmt % x})
    else:
        return (fmt % a)


def init_map_from_profile(grid, z_profile):
    """Use the given 1D density profile to initialize the 3D map"""
    grid[:,:,:] = z_profile[:]


def smooth_step_gaussian(x, xc, sigma):
    if (sigma == 0.0): return x*0.0
    return np.exp(-0.5*(x-xc)*(x-xc)/(sigma*sigma))


def get_mask_from_deformation(deformation_map, threshold=0.0, tolerance=1.0):
    """
    Compute a mask that gradually switches off the magnitude of a deformation
    map based on its displacement.

    Apply a smooth Gaussian switching function if the absolute value of the
    local deformation is less than the given threshold.

    Parameters
    ----------
    deformation_map : array_like
        Map of the deformation along the 2 bilayer dimensions.

    threshold : float
        Switch off the deformation map where the absolute value of the local
        deformation amplitude is below this value (default: 0.0, i.e. do not
        switch).

    tolerance : float
        Use this parameter as the tolerance of the switching function
        when the local deformation amplitude is less than "threshold" in
        absolute value (default: 1.0).
    """
    step = np.abs(deformation_map) - threshold
    smooth_step = np.exp(-0.5*np.square(step)/(tolerance*tolerance))
    smooth_step[step >= 0.0] = 1.0
    return smooth_step


def deform_map(grid_input, grid_output, amplitude,
               deformation_map=None,
               deformation_function=None, lengths=None):
    """
    Shift the density profile by a given table of deformation amplitudes.

    Parameters
    ----------

    grid_input : array_like
        Array with 3 dimensions, containing the input map.

    grid_output : array_like
        Array with 3 dimensions, containing the input.

    amplitude : array_like
        Displacement along the third axis of "grid", in units of grid points.
    """

    grid_output[:] = grid_input[:]
    amplitude_map = None
    if not deformation_map is None:
        amplitude_map = deformation_map * amplitude
    if not deformation_function is None:
        amplitude_map = eval(deformation_function)(shape=grid_input.shape[0:2],
                                                   amplitude=amplitude,
                                                   lengths=lengths)

    map_bend(grid=grid_output, deformation_map=amplitude_map)


def map_bend(grid, deformation_map):
    """
    Shift the density profile by a given table of deformation amplitudes.

    Parameters
    ----------

    grid : array_like
        Array with 3 dimensions, containing the input/output map; its values
        are modified in place.

    deformation_map : array_like
        Array with 2 dimensions equal to the first 2 dimensions of "grid"; its
        values are the deformation amplitudes, measured as displacements along
        the third axis of "grid", in units of grid points.
    """
    n = grid.shape
    z_local = np.zeros(shape=(n[2]))
    for ix in range(0, n[0]):
        for iy in range(0, n[1]):
            spnd.shift(input=grid[ix,iy,:], shift=deformation_map[ix,iy],
                       output=z_local)
            grid[ix,iy,:] = z_local[:]


def map_bend_gaussian(shape, amplitude, lengths, output=None):
    """
    Generate a Gaussian deformation, centered in the midpoint.

    Parameters
    ----------
    shape : array_like
        2 values, measuring the number of points along each axis.

    amplitude : float
        Height of the Gaussian in units of grid points

    lengths : array_like
        2 values, measuring the two Gaussian sigma parameters, in units of
        the interval length along each axis.

    output : array
        Array with 2 dimensions, containing the amplitudes of deformation.
        Its mean value is set to zero.
    """
    n = shape
    sigmas = lengths
    if output is None:
        output = np.zeros(shape=shape)
    for ix in range(0, n[0]):
        x = (float(ix) - 0.5*float(n[0]) + 0.5)/float(n[0])
        exp_x = np.exp(-0.5*x*x/(sigmas[0]*sigmas[0]))
        for iy in range(0, n[1]):
            y = (float(iy) - 0.5*float(n[1]) + 0.5)/float(n[1])
            exp_y = np.exp(-0.5*y*y/(sigmas[1]*sigmas[1]))
            output[ix,iy] = amplitude*exp_x*exp_y
    output -= output.mean()
    return output


def cos_x_function(x, period):
    """
    Cosine function, extended analytically outside [-period/2:period/2]
    """
    cosine = np.cos(2.0*np.pi/period*x)
    if ((x < -0.5*period) or (x > 0.5*period)): cosine = -1.0
    return cosine


def map_bend_cosine(shape, amplitude, lengths, output=None):
    """
    Generate a cosine deformation, centered on the midpoint.

    Parameters
    ----------
    shape : array_like
        Array with 2 values, measuring the number of points along each
        axis.

    amplitude : float
        Height of the cosine function in units of grid points

    lengths : array_like
        2 values, measuring the periods L of the cos(x) and cos(y) terms in
        units of the interval length along each axis; if smaller than 1.0 the
        respective cosine function is extended analytically outside of the
        [-L/2:L/2] interval.

    output : array
        Array with 2 dimensions, containing the amplitudes of the deformation.
        Its mean value is set to zero.
    """
    n = shape
    if output is None:
        output = np.zeros(shape=shape)
    period = lengths[0]
    for ix in range(0, n[0]):
        x = (float(ix) - 0.5*float(n[0]) + 0.5)/float(n[0])
        output[ix,:] = amplitude * cos_x_function(x, period)
    output -= output.mean()
    return output


def create_switching_grid(shape, map_minmax, map_domain_switching):
    """Generate a switching function that is 1 only within the given domain"""
    n = shape
    f = []
    for d in range(0, 3):
        dx = 1.0/float(n[d])
        x = np.arange(-0.5+0.5*dx, 0.5+0.5*dx, dx)
        f_x = np.ones(shape=(n[d]), dtype=float)
        min_x = map_minmax[0][d]
        max_x = map_minmax[1][d]
        switch_x = map_domain_switching[d]
        lower = np.where(x < min_x)
        f_x[lower] *= (smooth_step_gaussian(x, min_x, switch_x))[lower]
        upper = np.where(x > max_x)
        f_x[upper] *= (smooth_step_gaussian(x, max_x, switch_x))[upper]
        f += [f_x]
    f_xy = np.tensordot(f[0], f[1], axes=0)
    f_xyz = np.tensordot(f_xy, f[2], axes=0)
    return f_xyz


def read_xsc_cell(xsc_file):
    xsc = np.loadtxt(xsc_file)
    if xsc.ndim == 1:
        return xsc
    else:
        return xsc.mean(axis=0)


def get_cell_dimensions(args):
    cell = np.ones(shape=(3))
    defined = False
    if not args.xsc_file is None:
        xsc = read_xsc_cell(args.xsc_file)
        cell[:] = np.array([xsc[1], xsc[5], xsc[9]])
        defined = True
    if not args.cell is None:
        defined = True
        for i in range(3):
            if args.cell[i] == 0.0:
                continue
            if args.cell[i] <= 1.0:
                cell[i] *= args.cell[i]
            else:
                cell[i] = args.cell[i]
    if not defined:
        return None
    return cell


def create_grid_object(cell, args):
    """
    Create a new 3D map using gridData
    """
    # Set up 3D map
    cell_origin = -0.5*cell
    cell_dx = np.zeros(shape=(3))
    cell_n = np.zeros(shape=(3), dtype=int)
    if not args.grid_npoints is None:
        cell_dx = np.divide(cell, np.array(args.grid_npoints))
        cell_n = np.array(args.grid_npoints, dtype=int)
    else:
        raise Exception("Undefined number of grid points.")

    cell_origin += 0.5*cell_dx

    map_obj = gridData.Grid(np.ones(shape=cell_n, dtype=float),
                            origin=cell_origin,
                            delta=cell_dx)
    return map_obj


def init_z_profiles(origin, length, dz, args):

    z_axis = np.arange(origin+dz/2.0, origin+dz/2.0+length, dz)

    z_profiles = []
    z_profiles_coms = []

    if (len(args.input_profiles)):
        for input_file in args.input_profiles:
            # Load 1D profiles from file
            z_profile = np.loadtxt(input_file,
                                   dtype=[('z', np.float), ('p', np.float)])
            z_profiles += [z_profile]

    elif (len(args.input_centers)):
        for (input_center, input_sigma) in zip(args.input_centers,
                                               args.input_sigmas):
            # Set 1D profiles from a Gaussian
            z = np.linspace(-72.0, 72.0, 145)
            p = np.exp(-0.5*np.square(z-input_center) / \
                       np.square(input_sigma)) * \
                1.0/(np.sqrt(2.0*np.pi)*input_sigma)
            z_profile = {}
            z_profile['z'], z_profile['p'] = z, p
            z_profiles += [z_profile]

    else:
        raise Exception("No input density profiles provided.")


    for i_z, (z_profile, input_label) in enumerate(zip(z_profiles,
                                                       args.input_labels)):

        # Map the 1D profile to the Z-points of the 3D map
        z_min, z_max = z_profile['z'].min(), z_profile['z'].max()
        z_profile_spline = InterpolatedUnivariateSpline(z_profile['z'],
                                                        z_profile['p'])

        z_profile['p'] = z_profile_spline(z_axis)
        z_profile['z'] = z_axis

        # The spline is meaningless outside the original boundaries: trim it
        z_profile['p'][z_axis < z_min] = 0.0
        z_profile['p'][z_axis > z_max] = 0.0

        if (args.symmetrize_input):
            z_profile['p'] = (z_profile['p'] + \
                              z_profile_spline(-1.0 * z_axis)) / 2.0

        z_profile_com = spnd.center_of_mass(z_profile['p'])[0] * dz + \
            origin + dz/2.0
        print("Center of mass of input profile for \"%s\" =  %8.3f" % \
              (input_label, z_profile_com))
        z_profiles_coms += [z_profile_com]

    return (z_profiles, z_profiles_coms)


def select_map_domain(grid, map_domain, map_domain_switching,
                      select_inside_domain=True):
    """
    Multiply the given grid by a Gaussian-step switching function.

    Parameters
    ----------
    grid : array_like
        Array with 3 dimensions

    map_domain : array_like
        Array of the form [xmin, ymin, zmin, xmax, ymax, zmax]

    map_domain_switching : array_like
        Lengths of the switching step in each direction (3 elements)

    select_inside_domain : bool
        Preserve the map inside the domain if True, and outside the
        domain if False [default: True]
    """

    map_minmax = map_domain
    switching = create_switching_grid(grid.shape,
                                      map_minmax,
                                      map_domain_switching)
    if (not select_inside_domain):
        switching = (1.0 - switching)

    grid[:] *= switching[:]


def realign_map(grid, cell_origin, cell_dx, com):
    """
    Shift the map so that its center of mass is aligned to the given vector

    Parameters
    ----------
    grid : array_like
        Array with 3 dimensions

    cell_origin : array_like
        [x, y, z] coordinates of the unit cell's origin

        cell_origin : array_like
        [x, y, z] coordinates of the unit cell's origin


    """
    grid_com = np.multiply(spnd.center_of_mass(grid),
                           np.diag(cell_dx)) + \
                           cell_origin + np.diag(cell_dx)/2.0
    print("    Center of mass =", format_array(grid_com))
    shift_real = com - grid_com
    shift_map(grid=grid, cell_dx=cell_dx, shift_real=shift_real)
    print("    New center of mass =", format_array(grid_com + shift_real))


def shift_map(grid, cell_dx, shift_real):
    print("    Shift vector = ", format_array(shift_real))
    shift_bins = np.divide(shift_real, np.diag(cell_dx))
    spnd.shift(input=grid, shift=shift_bins, output=grid)


# Superseded in recent gridData
def write_map_to_file(map_object, output_file):
    """
    Write a gridData.Grid object to a .dx file suitable for GridForces.

    The GridForces DX parser only recognizes the "double" keyword for type,
    but different versions of gridData.Grid use "float" or ""double"".

    Parameters
    ----------
    map_object : object of gridData.Grid class

        This object is written using its native .export() function to a
        temp file, from which incompatible fields are replaced.

    output_file : str
        File name.
    """

    if gridData.__version__ >= '0.5.0':
        map_object.export(output_file, typequote='')
    else:
        # Workaround for older gridData verions
        import uuid
        tmpfile = output_file+"."+str(uuid.uuid4())+".dx"
        map_object.export(tmpfile)
        with open(tmpfile) as in_file:
            out_file = open(output_file, 'w')
            past_header = False
            for line in in_file:
                if (not past_header):
                    if (line[0:20] == 'object 3 class array'):
                        # So much for standard formats!
                        line = line.replace('float', 'double')
                        line = line.replace('"double"', 'double')
                        past_header = True
                out_file.write(line)
            out_file.close()
            os.remove(tmpfile)



def read_pdb_occ_beta(pdb_file):
    import MDAnalysis
    import warnings
    # We won't use PDB file as topology: disable this warning
    warnings.filterwarnings('ignore', category=UserWarning, message='Failed to '
                            'guess the mass for the following atom types:.*',
                            module='MDAnalysis')
    uni = MDAnalysis.Universe(pdb_file)
    if not uni == None:
        return (uni.atoms.occupancies, uni.atoms.tempfactors)
    return None


def get_gridforces_map_cmds(label, map_file, pdb_file,
                            gridforces_col='O', gridforces_charge_col='B'):
    return """
mGridForcePotFile       %s %s
mGridForceFile          %s %s
mGridForceCol           %s %s
mGridForceChargeCol     %s %s
# Disable this grid's scale factors; the Colvars module will compute the force
mGridForceScale         %s 0.0 0.0 0.0
# Continuously interpolate between periodic images of the map
mGridForceCont1         %s yes
mGridForceCont2         %s yes
mGridForceCont3         %s yes
# Do not raise an error if the grid extends beyond the simulation cell
mGridForceCheckSize     %s off
""" % (label, map_file, label, pdb_file,
       label, gridforces_col, label, gridforces_charge_col,
       label, label, label, label, label)


def get_vmd_load_map_cmds(map_file):
    return """\
mol addfile %s type dx waitfor all
""" % map_file


def multimap_colvar_def_tcl(name, map_labels, coefficients=[], map_norms=[],
                            indices=None, pdb_files=None, use_pdb_weights=False,
                            scripted_function=None, width=None):
    """
    Write a Tcl script to define the collective variable in VMD or NAMD.

    NAMD 2.14 supports mapName, whereas VMD supports only mapID.  Only one of
    the two sets is defined when writing the script for the specific back-end.

    Parameters
    ----------
    name : string
        Name of the collective variable, e.g. "mmcv"

    map_labels : list of strings
        If defined, use it as the values of the mapName Colvars keyword

    coefficients : list of floats
        If defined, use it as the values of the componentCoeff keyword

    map_norms : list of floats
        If defined, divide the componentCoeff keyword by this factor

    indices : list of integers
        If defined, use it as the values of the mapID Colvars keyword

    pdb_files : list of strings
        If defined, use it to select atoms within Colvars (atomsFile keyword)

    use_pdb_weights : bool
        If True, set the atomWeights keyword using PDB columns

    scripted_function : string
        If defined, use it as the value of the scriptedFunction keyword

    width : string
        If defined, set the width parameter of the colvar (the string may
        represent either a number or a Tcl variable expansion)

    """
    if len(coefficients) == 0 and scripted_function is None:
        raise Exception('Need either a set of coefficients or a scripting '
                        'function.')

    conf = """

# Define the Multi-Map collective variable
cv config \""""
    conf += """
colvar {

    name %s
""" % name

    if not width is None:
        conf += """\
    width %s
""" % str(width)

    if not scripted_function is None:
        conf += """\
    scriptedFunction %s
""" % scripted_function

    # Loop over the maps
    for i, label in enumerate(map_labels):

        scale = 1.0
        if len(map_norms) > 0:
            scale = 1.0/map_norms[i]

        conf += """
    mapTotal {
"""
        if not scripted_function is None:
            # Prefix the names by a numeric index to sort reliably
            conf += """\
        name %03d_%s
""" % (i, label)
        else:
            conf += """\
        name %s
        componentCoeff %13.6e
""" % (label, coefficients[i] * scale)

        if indices is None:
            conf += """\
        mapName %s
""" % label
        else:
            conf += """\
        mapID %d
""" % indices[i]

        if not pdb_files is None:
            conf += """\
        atoms {
            atomsFile %s
            atomsCol O
        }
""" % pdb_files[i]
            if use_pdb_weights:
                conf += """\
        atomWeights $weights(%s)
""" % pdb_files[i]

        conf += """\
    }
"""

    conf += """\
}
"
"""
    return conf


def singlemap_colvars_def_tcl(map_labels, coefficients=[], map_norms=[],
                              indices=None, pdb_files=None,
                              use_pdb_weights=False):
    conf = """
# Define single-map variables for diagnostics (no extra computational cost)
"""

    for i, label in enumerate(map_labels):

        scale = 1.0
        if len(map_norms) > 0:
            scale = 1.0/map_norms[i]

        conf += """
cv config \"
colvar {
    name %s
    mapTotal {
        name %s
""" % (label, label)
        if indices is None:
            # Map labels for NAMD
            conf += """\
        mapName %s
""" % label
        else:
            # Numeric IDs for VMD
            conf += """\
        mapID %d
""" % indices[i]

        if scale != 1.0:
            conf += """\
        componentCoeff %13.6e
""" % scale

        if not pdb_files is None:
            conf += """\
        atoms {
            atomsFile %s
            atomsCol O
        }
""" % pdb_files[i]
            if use_pdb_weights:
                conf += """\
        atomWeights $weights(%s)
""" % pdb_files[i]

        conf += """\
    }
}
\"
"""

    return conf


def namd_com_restraint_def(pdb_file, name='com_dist', force_constant=5.0):
    conf = """

# Define a center-of-mass restraint on all requested atoms
cv config \"
colvar {
    name %s
    width 0.1
    distance {
        group1 {
            atomsFile %s
            atomsCol O
        }
        group2 { dummyAtom (${com("x")}, ${com("y")}, ${com("z")}) }
    }
}

harmonic {
    name r_%s
    colvars %s
    centers 0.0
    forceConstant %.6f
}
\"
""" % (name, pdb_file, name, name, force_constant)
    return conf


def namd_ori_restraint_def(pdb_file, name='ori', force_constant=5.0):
    conf = """

# Define an orientation restraint on requested atoms
# WARNING: unlike the center of mass, this computation is not parallelized,
# and reducing the size of the selection (e.g. only C-alpha atoms for 
# proteins) is highly recommended to maintain reasonable parallel performance.
cv config \"
colvar {
    name %s
    width 0.01
    orientation {
        atoms {
            atomsFile %s
            atomsCol O
        }
        refPositionsFile %s
    }
}

harmonic {
    name r_%s
    colvars %s
    centers (1.0, 0.0, 0.0, 0.0)
    forceConstant %.6f
}
\"
""" % (name, pdb_file, pdb_file, name, name, force_constant)
    return conf


def namd_com_z_restraint_def(pdb_files, name='com_dist', force_constant=5.0):
    unique_files = list(set(pdb_files))
    conf = """

# Define a center-of-mass restraint on all requested atoms
cv config \"
colvar {
    name %s
    width 0.1
""" % name

    for pdb_file in unique_files:
        conf += """\
    distanceZ {
        componentCoeff %f
        main {
            atomsFile %s
            atomsCol O
        }
        ref { dummyAtom (0.0, 0.0, 0.0) }
    }
""" % (1.0/len(unique_files), pdb_file)

    conf += """\
}

harmonic {
    name r_%s
    colvars %s
    centers 0.0
    forceConstant %.6f
}
\"
""" % (name, name, force_constant)
    return conf


def multimap_colvar_function_tcl(args):
    (args.mmcv_slope, args.mmcv_intercept) = (1.0, 0.0)
    return """
proc calc_mmcv_base { args } {

    global mmcv_points
    global mmcv_length
    global mmcv_slope
    global mmcv_intercept

    set n ${mmcv_length}

    set xi 0.0
    set map_total 0.0
    set grads [list]

    for { set i 0 } { $i < ${n} } { incr i } {
        set x [lindex ${args} ${i}]
	set xi_i [lindex ${mmcv_points} ${i}]

        set xi [expr ${xi} + ${x}*${xi_i}]

        set ngrad [expr ${mmcv_slope}/${n} * ${xi_i} ]
        lappend grads ${ngrad}
    }
    set xi [expr ${mmcv_intercept} + ${mmcv_slope}*${xi}/${n}]

    if { [lindex ${args} end] == 1 } {
        return ${grads}
    }

    return ${xi}
}

proc calc_mmcv { args } {
    return [calc_mmcv_base {*}${args} 0]
}

proc calc_mmcv_gradient { args } {
    return [calc_mmcv_base {*}${args} 1]
}

# Parameters for the scripted function

set mmcv_length %d
set mmcv_points [list %s]
set mmcv_slope %.6f
set mmcv_intercept %.6f
""" % (args.n_inputs * len(args.multimap_coefficients),
       " ".join(map(str, args.multimap_coefficients)),
       args.mmcv_slope, args.mmcv_intercept)



def parse_cmdline_args(namespace=None):
    """Parse command-line arguments with argparse, return args object."""

    import argparse

    parser = \
        argparse.ArgumentParser(description='''
Generate input files for Multi-Map computations with VMD and NAMD.  Reference article: https://doi.org/10.1002/jcc.26075''',
                                usage='%(prog)s [options] (see -h or --help)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group = parser.add_argument_group(title='General options')
    group.add_argument('--input-labels',
                       nargs='*', type=str,
                       help='Label each input selection to be treated '
                       'separately (e.g. "protein", or "upper" "lower")',
                       default=[])
    group.add_argument('--input-pdb-files',
                       nargs='*', type=str,
                       help='PDB files defining each selection; if the product '
                       'of the occupancy and beta columns is non-zero, the '
                       'atom contributes to the Multi-Map variable',
                       default=None)
    group.add_argument('--input-num-particles',
                       nargs='*', type=int,
                       help='If --normalize-maps is defined, use these numbers '
                       'as denominators to normalize; otherwise, the sums of '
                       'the occupancy*beta columns of each PDB file are used',
                       default=None)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--input-map-files',
                       nargs='*', type=str,
                       help='Files containining 3D maps (OpenDX format); '
                       'if not provided, a lipid bilayer workflow is assumed '
                       'and the maps must be generated by --input-profiles '
                       'or --input-centers',
                       default=[])
    group.add_argument('--input-profiles',
                       nargs='*', type=str,
                       help='Files containing 1D lipid bilayer number density '
                       'profiles (two-column text format)',
                       default=[])
    group.add_argument('--input-centers',
                       nargs='*', type=float,
                       help='Create Gaussian 1D profiles with these centers ',
                       default=[])

    group = parser.add_argument_group(title='Bilayer options for treating 1D profiles')
    group.add_argument('--input-sigmas',
                       nargs='*', type=float,
                       help='Use these sigmas in combination with '
                       '--input-centers.',
                       default=[])
    group.add_argument('--symmetrize-input',
                       action='store_true',
                       help='Symmetrize the 1D profile around z = 0; '
                       'requires that only one profile is given.',
                       default=False)
    group.add_argument('--align-com',
                       action='store_true',
                       help='Shift each 3D map along the Z axis '+
                       'to preserve its center of mass',
                       default=False)

    # Define the bilayer deformation
    group = parser.add_argument_group(title='Bilayer options: load or define '
                                      'the 2D deformation map')
    group.add_argument('--deformation-function',
                       type=str,
                       help='Name of the function to use for deformation',
                       choices=['map_bend_cosine', 'map_bend_gaussian'])
    group.add_argument('--deformation-lengths',
                       nargs=2, type=float, metavar=('LDef_X', 'LDef_Y'),
                       help='Lengths of the deformation along the X and Y '
                       'axes if used with --deformation-function '
                       '(in fractional coordinates along X and Y axes)',
                       default=[1.0, 1.0])
    group.add_argument('--deformation-map',
                       type=str,
                       help='Name of a file containing the deformation map '
                       'in a format suitable for numpy.loadtxt '
                       '(must have the same shape as the first two '
                       'dimensions of the 3D maps); '
                       'mutually exclusive with --deformation-function.')
    group.add_argument('--deformation-amplitudes',
                       nargs='*', type=float,
                       help='Amplitudes of deformation',
                       default=[])

    group = parser.add_argument_group(title='Mask on/off certain regions of the 3D maps')
    group.add_argument('--deformation-map-mask',
                       type=str,
                       help='Name of a file containing the mask map '
                       'in the same format as --deformation-map. '
                       'May be used together with other masking functions.')
    group.add_argument('--mask-by-amplitude',
                       type=np.float,
                       help='Mask out regions that move less than this '
                       'fraction of the maximum amplitude.  '
                       'May be used together with other masking functions.',
                       default=0.0)
    group.add_argument('--map-domain',
                       nargs=6, type=float,
                       metavar=('Xmin', 'Xmax', 'Ymin', 'Ymax',
                                'Zmin', 'Zmax'),
                       help='Switch on/off the map outside these boundaries, '
                       'written in fractional coordinates of the unit cell',
                       default=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5])
    group.add_argument('--map-domain-switching',
                       nargs=3, type=float,
                       metavar=('Switch_X', 'Switch_Y', 'Switch_Z'),
                       help='Smoothing distance outside of map_domain '+
                       'written in fractional coordinates of the unit cell',
                       default=[0.05, 0.05, 0.05])
    group.add_argument('--map-domain-inside',
                       action='store_true',
                       help='Mask off the map outside the selected boundaries',
                       default=True)

    group = parser.add_argument_group(title='Recenter the 3D maps')
    group.add_argument('--shift-com',
                       nargs=3, type=float,
                       help='Shift the 3D map by this vector '
                       '(after --align-com).',
                       default=[0.0, 0.0, 0.0])

    group = parser.add_argument_group(title='Define 3D grid parameters for output maps')
    group.add_argument('--xsc-file',
                       type=str,
                       help='Read unit cell dimensions from this file '
                       '(Angstrom units, NAMD format).  '
                       'Triclinic cells are currently unsupported.')
    group.add_argument('--cell',
                       nargs=3, type=float, metavar=('L_x', 'L_y', 'L_z'),
                       help='Define the unit cell dimensions explicitly.  '
                       'If --xsc-file is also provided, values read from it '
                       'are redefined based on the values of this option as '
                       'follows: (1) if a value is zero, the corresponding '
                       'dimension from the XSC file is left unchanged; '
                       '(2) if a value is positive but smaller than 1, the '
                       'corresponding fraction of the original value is used; '
                       '(3) if a value is larger than 1, it overrides the '
                       'original value.')
    group.add_argument('--grid-npoints',
                       nargs=3, type=int, metavar=('n_x', 'n_y', 'n_z'),
                       help='Numbers of points of 3D grid (overrides '
                       '--grid-spacings); if either --input-profiles or '
                       '--deformation-map are given, default values are taken '
                       'from them.',
                       default=[0, 0, 0])

    group = parser.add_argument_group(title='Multi-Map collective variable definition')
    group.add_argument('--colvar-name',
                       type=str,
                       help='Name of the Multi-Map collective variable '
                       'being defined.',
                       default='multimap')
    group.add_argument('--multimap-coefficients',
                       nargs='*', type=float, metavar='Xi_k',
                       help='Coefficients of each map in the Multi-Map '
                       'equation; if not given, the index of the map is used, '
                       'starting from 1.')
    group.add_argument('--normalize-maps',
                       action='store_true',
                       help='Scale the coefficients of each map so that '
                       'the sum of its values over all atoms is 1 or less.',
                       default=False)
    group.add_argument('--use-pdb-weights',
                       action='store_true',
                       help='When true, multiply the occupancy and temperature '
                       'PDB columns to define the weights of each atom within '
                       'each map; when false, all weights are equal to 1. '
                       'Defaults to True if --input-pdb-files is used.',
                       default=None)
    group.add_argument('--com-restraint',
                       action='store_true',
                       help='Add the definition of a center-of-mass restraint '
                       'to the NAMD input',
                       default=True)

    group.add_argument('--ori-restraint',
                       action='store_true',
                       help='Add the definition of an orientational restraint '
                       'to the NAMD input; it is highly recommended to reduce '
                       'the size of the selection (e.g. only C-alpha atoms '
                       'for proteins) to maintain reasonable parallel '
                       'performance.  Not valid for bilayers.',
                       default=False)

    group = parser.add_argument_group(title='Output parameters')
    group.add_argument('-o', '--output-prefix',
                       type=str, required=True,
                       help='Prefix for output files')
    group.add_argument('--namd-script',
                       type=str, default=None,
                       help='Write NAMD commands to define the Multi-Map '
                       'collective variable into this file.')
    group.add_argument('--vmd-script',
                       type=str, default=None,
                       help='Write VMD commands to define the Multi-map '
                       'collective variable into this file.')

    args = parser.parse_args(namespace=namespace)

    return args


def check_prep_args(args):
    """
    Check the types of the input arguments and initialize their types
    """

    if args.align_com and (np.dot(args.shift_com, args.shift_com) > 0.0):
        raise(Exception("--align-com and --shift-com are mutually exclusive."))

    args.n_inputs = len(args.input_labels)

    if len(args.input_centers) != len(args.input_sigmas):
        raise(Exception(
            "--input-centers and --input-sigmas must have the same lengths."))

    if len(args.input_map_files) > 0:
        args.system_dim = '3d'
    else:
        args.system_dim = '2d'

    if len(args.input_map_files) > 0 and args.n_inputs > 1:
        raise(Exception("With --input-map-files, only one value is allowed "
                        "for --input-label."))

    if not args.input_num_particles is None:
        if len(args.input_num_particles) != args.n_inputs:
            raise(Exception("The length of --input-num-particles must match "
                            "the length of --input-labels."))

    if not args.input_pdb_files is None:
        if len(args.input_pdb_files) != args.n_inputs:
            raise(Exception("The length of --input-pdb-files must match "
                            "the number of input selections."))
        if args.use_pdb_weights is None:
            args.use_pdb_weights = True
        all_ones = True
        args.atom_weights = []
        for f in args.input_pdb_files:
            weights = np.multiply(*(read_pdb_occ_beta(f)))
            args.atom_weights += [weights]
            real_weigths = weights[weights != 0.0]
            if (real_weigths != 1.0).any():
                all_ones = False
        if all_ones and args.use_pdb_weights:
            print("Disabling --use-pdb-weights: all values are 1 anyway.")
            args.use_pdb_weights = False
        if args.input_num_particles is None:
            args.input_num_particles = [w.sum() for w in args.atom_weights]
    else:
        args.input_pdb_files = [None for i in range(args.n_inputs)]
        if args.use_pdb_weights is None:
            args.use_pdb_weights = False
        if args.use_pdb_weights:
            raise(Exception("--use-pdb-weights requires --input-pdb-files."))

    if not args.input_num_particles is None:
        # Convert to dictionary
        d = {}
        for i in range(args.n_inputs):
            d[args.input_labels[i]] = args.input_num_particles[i]
        args.input_num_particles = d

    if args.deformation_map and args.deformation_function:
        raise(Exception("--deformation-map and --deformation-function "
                        "are mutually exclusive."))

    args.grid_npoints = np.array(args.grid_npoints)

    if args.grid_npoints[2] == 0:
        if len(args.input_profiles) > 0:
            args.grid_npoints[2] = args.input_profiles[0].shape[0]
        else:
            args.grid_npoints[2] = 100

    if not args.deformation_map is None:
        print("Loading the deformation map from file \""+args.deformation_map+
              "\"; values of this 2D map will be multiplied by "
              "the deformation amplitudes.  "
              "The mean is automatically subtracted.")
        args.deformation_map = np.loadtxt(args.deformation_map, dtype=np.float)
        args.deformation_map_mask = np.ones(shape=args.deformation_map.shape)
        print("Subtracting the mean Z-coordinate of the map; before =",
              args.deformation_map.mean())
        args.deformation_map -= args.deformation_map.mean()
        print("Subtracting the mean Z-coordinate of the map; after  =",
              args.deformation_map.mean())
        if args.grid_npoints[0] == 0:
            args.grid_npoints[0] = args.deformation_map.shape[0]
        if args.grid_npoints[1] == 0:
            args.grid_npoints[1] = args.deformation_map.shape[1]

    if type(args.deformation_map_mask) == str:
        print("Loading the mask from file \""+args.deformation_map_mask+"\"")
        args.deformation_map_mask = np.loadtxt(args.deformation_map_mask,
                                               dtype=np.float)
        print("Mean value of the mask map =", args.deformation_map_mask.mean())

    if args.mask_by_amplitude > 0.0:
        if args.deformation_map is None:
            raise(Exception("--mask-by-amplitude requires --deformation-map."))
        print("Masking the deformation map based on its magnitude")
        max_amplitude = np.abs(args.deformation_map).max()
        mask_threshold = args.mask_by_amplitude*max_amplitude
        args.deformation_map_mask *= \
            get_mask_from_deformation(deformation_map=args.deformation_map,
                                      threshold=mask_threshold,
                                      tolerance=0.2*mask_threshold)
        print("Mean value of the mask map =", args.deformation_map_mask.mean())
        del mask_threshold, max_amplitude


    args.n_deformations = max(len(args.input_map_files),
                              len(args.deformation_amplitudes))

    if args.multimap_coefficients is None:
        args.multimap_coefficients = np.arange(1.0, args.n_deformations+1)

    args.multimap_coefficients = np.array(args.multimap_coefficients)

    if len(args.multimap_coefficients) != args.n_deformations:
        raise Exception("Numbers of arguments for --deformation-amplitudes must match the number of deformations requested.")

    if not args.deformation_lengths is None:
        args.deformation_lengths = np.array(args.deformation_lengths)


def generate_maps_from_profiles(args):
    """
    Bilayer handler
    """

    cell_sizes = get_cell_dimensions(args)
    # TODO forbid triclinic cells

    # Object to hold the undeformed map for each input
    map_input = create_grid_object(cell_sizes, args)
    print("Unit cell dimensions = ", format_array(cell_sizes))
    print("Grid sizes (# points) = ", map_input.grid.shape)
    # Define cell_dx in a consistent manner, because the shape of Grid.delta
    # changes between gridData versions (3-vector, or 3x3 matrix)
    cell_dx = map_input.delta
    if cell_dx.ndim == 1:
        cell_dx = np.diag(cell_dx)
    cell_origin = map_input.origin

    # Read 1D profiles, along the same grid of the z-axis of the 3D map
    z_profiles, z_profiles_coms = \
         init_z_profiles(origin=cell_origin[2], length=cell_sizes[2],
                         dz=cell_dx[2][2], args=args)

    maps = {}

    # Loop over all input maps (e.g. upper and lower leaflet)
    for z_profile, z_profile_com, input_label in \
        zip(z_profiles, z_profiles_coms, args.input_labels):

        maps[input_label] = []

        print("")
        print("Generating the set of maps with label: ", input_label)

        init_map_from_profile(z_profile=z_profile['p'], grid=map_input.grid)

        for i, amplitude in enumerate(args.deformation_amplitudes):

            # Object to hold the deformed map
            map_output = create_grid_object(cell_sizes, args)

            print("  Deformation amplitude =", format_array(amplitude))
            amplitude_bins = np.divide(amplitude, cell_dx[2][2])
            # print("  Amplitude in grid points =",
            # format_array(amplitude_bins))
            # if not args.deformation_function is None:
            #     print("  Deformation function =", args.deformation_function)
            # print("  Lengths (fractional units) =",
            #       format_array(args.deformation_lengths))

            deform_map(grid_input=map_input.grid,
                       grid_output=map_output.grid,
                       amplitude=amplitude_bins,
                       lengths=args.deformation_lengths,
                       deformation_function=args.deformation_function,
                       deformation_map=args.deformation_map)

            if args.align_com:
                print("  Aligning with the input profile:")
                realign_map(grid=map_output.grid,
                            cell_origin=cell_origin,
                            cell_dx=cell_dx,
                            com=[0.0, 0.0, z_profile_com])

            if not args.deformation_map_mask is None:
                factor = args.deformation_map_mask.sum() / \
                         np.float(args.deformation_map_mask.size)
                print("  Masking the map (fraction = %.1f%%)" % (factor*100))
                np.multiply(map_output.grid,
                            args.deformation_map_mask[:,:,np.newaxis],
                            out=map_output.grid)
                map_output.grid[:,:,:] *= (1.0/factor)

            maps[input_label] += [map_output]

    return maps


def write_namd_script(gridforces_script, map_labels, pdb_files, map_norms,
                      args):
    with open(args.namd_script, 'w') as namd_script:

        print("Writing NAMD script to file", args.namd_script)

        namd_script.write("""# -*- tcl -*-

## Please ensure that "mGridForce on" and "colvars on" are defined before
## sourcing this file in a NAMD script.

# The width parameters sets the amplitude of the colvar's fluctuation, and
# can be used to set the correct scale for harmonic restraint potentials
if { [info exists multimap_cv_width] == 0 } {
    set multimap_cv_width 1.0
}
""")

        if args.system_dim == '3d':
            namd_script.write("""
# Define a default position for the center-of-mass restraint
if { [info exists com_pos] == 0 } {
    set com_pos [list 0.0 0.0 0.0]
}

set com("x") [lindex ${com_pos} 0]
set com("y") [lindex ${com_pos} 1]
set com("z") [lindex ${com_pos} 2]
""")

        namd_script.write("""
# Load gridForces maps
""")
        namd_script.write(gridforces_script)

        scripted_function = None
        if not scripted_function is None:
            namd_script.write(multimap_colvar_function_tcl(args))

        mmcv_def = \
            multimap_colvar_def_tcl(name=args.colvar_name,
                                    map_labels=map_labels,
                                    map_norms=map_norms,
                                    coefficients=\
                                        np.tile(args.multimap_coefficients,
                                                args.n_inputs),
                                    scripted_function=scripted_function,
                                    use_pdb_weights=args.use_pdb_weights,
                                    width='${multimap_cv_width}')
        namd_script.write(mmcv_def)

        singles_def = \
            singlemap_colvars_def_tcl(map_labels=map_labels,
                                      map_norms=map_norms,
                                      use_pdb_weights=args.use_pdb_weights)
        namd_script.write(singles_def)

        if args.system_dim == '3d':
            # Center-of-mass restraint
            namd_script.write(namd_com_restraint_def(pdb_files[0]))
            if args.ori_restraint:
                namd_script.write(namd_ori_restraint_def(pdb_files[0]))

        if args.system_dim == '2d':
            namd_script.write(namd_com_z_restraint_def(pdb_files))



def write_vmd_script(vmd_load_map_cmds, map_labels, pdb_files, map_norms, args):
    unique_files = list(set(pdb_files))
    with open(args.vmd_script, 'w') as vmd_script_file:
        print("Writing VMD script to file", args.vmd_script)
        vmd_script_file.write("""# -*- tcl -*-

## VMD script to load map files onto the top molecule

# The width parameters sets the amplitude of the colvar's fluctuation, and
# can be used to set the correct scale for harmonic restraint potentials
if { [info exists multimap_cv_width] == 0 } {
    set multimap_cv_width 1.0
}

""")
        if args.use_pdb_weights:
            vmd_script_file.write("""
# Precompute the per-atom weights by taking them from each PDB file
foreach pdb_file { %s } {
    mol addfile $pdb_file type pdb waitfor all
    set pdb_file_sel [atomselect top "(occupancy != 0.0)"]
    set weights($pdb_file) [vecmul [$pdb_file_sel get occupancy] [$pdb_file_sel get beta]]
    ${pdb_file_sel} delete
    animate delete beg [expr [molinfo top get numframes] -1]
}
""" % ' '.join(set(unique_files)))

        vmd_script_file.write("""
# Attach Colvars to the top molecule
cv molid top

# Load volumetric maps
""")
        vmd_script_file.write(vmd_load_map_cmds)
        scripted_function = None
        if not scripted_function is None:
            vmd_script_file.write(multimap_colvar_function_tcl(args))

        n = len(map_labels)

        mmcv_def = \
            multimap_colvar_def_tcl(name=args.colvar_name,
                                    map_labels=map_labels,
                                    map_norms=map_norms,
                                    indices=range(n),
                                    pdb_files=pdb_files,
                                    coefficients=\
                                        np.tile(args.multimap_coefficients,
                                                args.n_inputs),
                                    scripted_function=scripted_function,
                                    use_pdb_weights=args.use_pdb_weights,
                                    width='${multimap_cv_width}')
        vmd_script_file.write(mmcv_def)

        singles_def = \
            singlemap_colvars_def_tcl(map_labels=map_labels,
                                      indices=range(n),
                                      map_norms=map_norms,
                                      pdb_files=pdb_files,
                                      use_pdb_weights=args.use_pdb_weights)
        vmd_script_file.write(singles_def)


def generate_multimap(args):
    """
    Process cmdline arguments, generate/load the maps, write input scripts
    """

    if not args.namd_script is None:
        # Begin writing the NAMD script
        gridforces_script = ""

    if not args.vmd_script is None:
        # Begin writing the VMD script
        vmd_load_map_cmds = ""

    map_files = {}

    if len(args.input_profiles) > 0 or len(args.input_centers) > 0:
        maps = generate_maps_from_profiles(args)
    else:
        maps = {}
        # Convert to dictionary
        for input_label in args.input_labels:
            maps[input_label] = []
            for input_map_file in args.input_map_files:
                print("Loading map from file", input_map_file)
                maps[input_label] += [gridData.Grid(input_map_file)]
            map_files[input_label] = args.input_map_files

    map_labels = {}
    pdb_files = {}
    map_norms = {}

    # Loop over all input labels (e.g. upper and lower leaflet)
    for input_label, pdb_file in zip(args.input_labels, args.input_pdb_files):

        print("")
        print("Processing the set of maps with label: ", input_label)

        if not args.namd_script is None:
            gridforces_script += """
# Maps for selection: \"%s\"
""" % input_label

        if not args.vmd_script is None:
            vmd_load_map_cmds += """
# Maps for selection: \"%s\"
""" % input_label

        if not input_label in map_files:
            map_files[input_label] = []

        map_labels[input_label] = []
        pdb_files[input_label] = []

        corr = np.zeros(shape=(args.n_deformations, args.n_deformations))
        integrals = np.zeros(shape=(args.n_deformations))

        # Loop over all maps for this selection/label
        for i, map_output in enumerate(maps[input_label]):

            map_label = input_label
            if args.n_deformations > 1:
                map_label += ("_%03d" % (i+1))
            map_labels[input_label] += [map_label]

            if len(args.input_map_files) == 0:
                if np.dot(args.shift_com, args.shift_com) > 0.0:
                    print("  Shifting the center of mass:")
                    shift_map(grid=map_output.grid, shift=args.shift_com)

                if not args.map_domain is None:
                    if args.map_domain != [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5]:
                        args.map_domain = \
                            np.array(args.map_domain).reshape((3, 2)).T
                        print("  Selecting region using the following boundaries:")
                        print('  ',
                              format_array(args.map_domain).replace(os.linesep, ' '))
                        select_map_domain(grid=map_output.grid,
                                          map_domain=args.map_domain,
                                          map_domain_switching=args.map_domain_switching,
                                          select_inside_domain=args.map_domain_inside)

            # Name of PDB file used to select atoms for this map; NAMD needs
            # it for GridForces, VMD uses it as selection within Colvars
            if pdb_file is None:
                pdb_file = args.output_prefix+'.'+input_label+'.pdb'
            pdb_files[input_label] += [pdb_file]

            if len(args.input_map_files) == 0:
                # New map files
                map_file = args.output_prefix+'.'+map_label+'.dx'
                map_files[input_label] += [map_file]
            else:
                map_file = map_files[input_label][i]

            if not args.namd_script is None:
                gridforces_script += \
                    get_gridforces_map_cmds(label=map_label,
                                            map_file=map_file,
                                            pdb_file=pdb_file)

            if not args.vmd_script is None:
                vmd_load_map_cmds += \
                    get_vmd_load_map_cmds(map_file=map_file)

            integrals += [np.sum(map_output.grid)]

            if len(args.input_map_files) == 0:
                # Compute CCC only if maps were generated here (may have
                # different grids)
                j = i
                while (j >= 0):
                    corr[i, j] = np.sum(np.multiply(maps[input_label][i].grid,
                                                    maps[input_label][j].grid))
                    corr[i, j] /= np.sum(np.square(maps[input_label][j].grid))
                    corr[j, i] = corr[i, j]
                    j -= 1

        if len(args.input_map_files) == 0:
            print("Cross-correlation matrix between maps:")
            if args.n_deformations > 2:
                print("  (elements far from the diagonal should be small)")
            print(format_array(corr))
            for m, map_file in zip(maps[input_label], map_files[input_label]):
                print("  Writing to file \""+map_file+"\"")
                write_map_to_file(m, map_file)

        if args.normalize_maps:
            map_norms[input_label] = [m.grid.max() * \
                                      args.input_num_particles[input_label] \
                                          for m in maps[input_label]]
        else:
            map_norms[input_label] = []

    all_map_labels = [x for il in args.input_labels for x in map_labels[il]]
    all_pdb_files = [x for il in args.input_labels for x in pdb_files[il]]
    all_map_norms = [x for il in args.input_labels for x in map_norms[il]]

    print("")

    if not args.namd_script is None:
        write_namd_script(gridforces_script, all_map_labels, all_pdb_files,
                          all_map_norms, args)

    if not args.vmd_script is None:
        write_vmd_script(vmd_load_map_cmds, all_map_labels, all_pdb_files,
                         all_map_norms, args)


if __name__ == '__main__':
    args = parse_cmdline_args()
    check_prep_args(args)
    generate_multimap(args)
