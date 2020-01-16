#!/usr/bin/env python

# Invoke this script with --help for documentation.

# Download link: https://github.com/Colvars/colvars/blob/master/colvartools/plot_colvars_traj.py?raw=true
# Contact: giacomo.fiorin@gmail.com

from __future__ import print_function

import os
import sys

if (sys.version_info[:2] < (2, 7)):
    # Save some explanations
    print("Python versions prior to 2.7 are no longer supported.")
    sys.exit(1)

import numpy as np


class Colvar_traj(object):
    """
    Class holding the trajectory of a collective variable.
    The series of step numbers are included as well, because collective
    variables may be added or deleted during a simulation.
    """

    _name = ""
    _step = np.empty(shape=(0))
    _colvar = np.empty(shape=(0))

    def __init__(self, name):
        """Sets the name of the variable"""
        self._name = name
        self._step = np.zeros(shape=(0), dtype=np.int64)
        self._colvar = np.zeros(shape=(0), dtype=np.float)

    def __len__(self):
        """Returns the length of the trajectory"""
        return len(self._step)

    def __str__(self):
        """String representation of the trajectory"""
        return """Trajectory of collective variable '%s':
    Steps = [%d ... %d]
    Values = [%f ... %f]""" % (self._name, self._step[0], self._step[-1],
                               self._colvar[0], self._colvar[-1])

    def _set_num_dimensions(self, n_d):
        """Set the number of components of the collective variable"""
        if (len(self) > 0):
            print("Warning: changing the number of dimensions "
                  "of collective variable \""+self._name+
                  "\" after it has already been read.")
        if (n_d > 1):
            self._colvar.resize((len(self), n_d))
        else:
            self._colvar.resize((len(self)))

    def _resize(self, n):
        """Change the number of records in the trajectory"""
        self._step.resize((n))
        if (len(self._colvar.shape) > 1):
            self._colvar.resize((n, self._colvar.shape[1]))
        else:
            self._colvar.resize((n))

    @property
    def name(self):
        """The name of the collective variable"""
        return self._name

    @property
    def num_dimensions(self):
        """Dimensionality of the collective variable (d)"""
        s = self._colvar.shape
        if (len(s) > 1):
            return s[1]
        else:
            return 1

    @property
    def num_frames(self):
        """Number of trajectory frames loaded (n)"""
        return len(self._step)

    @property
    def steps(self):
        """The array of step numbers, with shape = (n,)"""
        return self._step

    @property
    def values(self):
        """The array of collective variable values, with shape = (n, d)"""
        return self._colvar


class Colvars_traj(object):
    """Trajectories of collective variables, as read from a list of
    colvars.traj files.
    Can be accessed as a dictionary using a variable's name as key; each
    variable's trajectory is an instance of colvars_traj"""

    _keys = []
    _start = {}
    _end = {}
    _colvars = {}
    _found = {}
    _frame = -1

    def __init__(self, filenames=None, first=0, last=-1, every=1):
        """
        Initialize from the given list of colvars.traj files
        Any optional arguments are passed to read_files()
        """
        self._keys = ['step']
        self._start['step'] = 0
        self._end['step'] = -1
        self._count = 0
        self._frame = 0
        if filenames:
            self.read_files(filenames=filenames, first=first, last=last,
                            every=every)

    def __getitem__(self, key):
        return self._colvars[key]

    def __contains__(self, key):
        return key in self._colvars

    def __str__(self):
        return """Set of collective variable trajectories:
    Variables = %s
    Number of frames = %d""" % (self.variables, self.num_frames)

    @property
    def num_frames(self):
        """Number of trajectory frames loaded"""
        return self._frame

    @property
    def variables(self):
        """Names of variables defined"""
        return self._keys[1:] # The first entry is "step"

    def _parse_comment_line(self, line):
        """
        Read in a comment line from a colvars.traj file and update the names of
        collective variables according to the contents of the comment line.
        """
        new_keys = (line[1:]).split() # Skip the hash char at position 0
        if (new_keys[0] != 'step'):
            raise KeyError("Error: file format incompatible with colvars.traj")
        self._keys = new_keys
        # Find the boundaries of each column
        for i in range(1, len(self._keys)):
            self._start[self._keys[i]] = line.find(' '+self._keys[i]+' ')
            self._end[self._keys[i-1]] = line.find(' '+self._keys[i]+' ')
            self._end[self._keys[-1]] = -1

    def _parse_line(self, line):
        """
        Read in a data line from a colvars.traj file
        """

        step = np.int64(line[0:self._end['step']])

        for v in self._keys[1:]:
            text = line[self._start[v]:self._end[v]]
            v_v = np.fromstring(text.lstrip(' (').rstrip(') '), sep=',')
            n_d = len(v_v)
            if (v not in self._colvars):
                self._colvars[v] = Colvar_traj(v)
                self._colvars[v]._set_num_dimensions(n_d)
            cv = self._colvars[v]
            n = len(cv)
            cv._resize(n+1)
            cv.steps[n] = step
            cv.values[n] = v_v


    def read_files(self, filenames, list_variables=False,
                   first=0, last=-1, every=1):
        """
        Read a series of colvars.traj files.
        filenames : list of strings
            list of file names to read from
        list_variables : bool
            if true, return the list of variable names defined before
            the first data line
        first : int
            skip all records before this index (see also mol load in VMD)
        last : int
            stop reading after this record
        every : int
            read every these many records
        """
        last = np.int64(last)
        if (last == -1):
            last = np.int64(np.iinfo(np.int64).max)
        last_step = -1
        for f in [open(filename) for filename in filenames]:
            for line in f:
                if (len(line) == 0): continue
                if (line[:1] == "@"): continue # xmgr file metadata
                if (line[:1] == "#"):
                    self._parse_comment_line(line)
                    continue
                if list_variables:
                    return self.variables
                step = np.int64(line[0:self._end['step']])
                if (step == last_step): continue
                if ((self._frame >= first) and (self._frame <= last) and
                    (self._frame % every == 0)):
                    self._parse_line(line)
                self._frame += 1
                last_step = step


if (__name__ == '__main__'):

    import argparse

    parser = \
        argparse.ArgumentParser(description='Select variables from one or more '
                                'Colvars trajectory files, and write them to '
                                'a concatenated file, or optionally plot them '
                                '(1D graph) as a function of time or of '
                                'one of the variables.'+os.linesep+
                                'Plotting functions are only available when '
                                'Matplotlib is available.', \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(dest='filenames',
                        nargs='*',
                        type=str,
                        help='Space-separated list of .colvars.traj files '
                        '(will be concatenated)',
                        default=[])

    parser.add_argument('--dt',
                        type=float,
                        help='Integration time step',
                        default=2.0e-6)

    parser.add_argument('--first-frame',
                        dest='first',
                        type=np.int64,
                        help='First frame to read',
                        default=0)

    parser.add_argument('--last-frame',
                        dest='last',
                        type=np.int64,
                        help='Last frame to read',
                        default=-1)

    parser.add_argument('--skip-frames',
                        dest='skip',
                        type=np.int64,
                        help='Read every these many frames',
                        default=1)

    parser.add_argument('--variables',
                        nargs='*',
                        type=str,
                        help='Space-separated list of names of collective '
                        'variables to write or plot.',
                        default=[])

    parser.add_argument('--list-variables',
                        action='store_true',
                        help='List all names of collective variables '
                        'defined up until the first line of data.',
                        default=False)

    parser.add_argument('--output-file',
                        type=str,
                        help='Write the selected variables to a text file.  '
                        'The step number is always included as the first '
                        'column, and all variables '
                        'must be defined on the same trajectory segments.',
                        default=None)

    parser.add_argument('--plot-file',
                        type=str,
                        help='Plot into a file with this name if the '
                        'extension is one of "png", "pdf", "ps", "eps", '
                        'or "svg", or to the screen if not given '
                        'and --output-file is also not given.',
                        default=None)

    parser.add_argument('--plot-x-axis',
                        type=str,
                        help='Use this variable as X axis in the plot.',
                        default='time')

    parser.add_argument('--plot-x-label',
                        type=str,
                        help='Use this label for the X axis.',
                        default=None)

    parser.add_argument('--plot-y-label',
                        type=str,
                        help='Use this label for the Y axis.',
                        default=None)

    parser.add_argument('--plot-keys',
                        nargs='*',
                        type=str,
                        help='Alternative names for the legend',
                        default=[])

    parser.add_argument('--plot-time-factor',
                        type=int,
                        help='Average these many consecutive frames '
                        'together before plotting but after reading '
                        'from files',
                        default=1)

    args = parser.parse_args()

    if (len(args.filenames) == 0):
        raise Exception("No filenames provided.")

    colvars_traj = Colvars_traj()
    r = colvars_traj.read_files(args.filenames,
                                list_variables=args.list_variables,
                                first=args.first,
                                last=args.last,
                                every=args.skip)

    if (args.list_variables):
        print(r)
        sys.exit()

    variables = args.variables

    for var in variables:
        try:
            cv = colvars_traj[var]
        except KeyError:
            print("ERROR: could not find \""+var+"\"; below are the available "
                  "fields: ")
            print("  ", *colvars_traj.variables)
            raise KeyError

    if (len(variables) == 0): variables = colvars_traj.variables

    plot_keys = dict(zip(variables, variables))


    if (args.output_file):

        fmt = " %12d"
        columns = [colvars_traj[variables[0]].steps]
        for var in variables:
            cv = colvars_traj[var]
            for ic in range(cv.num_dimensions):
                y = cv.values
                if (cv.num_dimensions > 1):
                    y = cv.values[:,ic]
                columns += [y]
                fmt += " %21.14f"
        columns = tuple(columns)
        output_file = args.output_file
        if (output_file == '-'): output_file = sys.stdout
        np.savetxt(fname=output_file,
                   X=list(zip(*columns)),
                   fmt=str(fmt))

    else:

        plot_filename = None
        if (args.plot_file):
            plot_filename = args.plot_file
            if (not args.plot_file.lower().endswith(('png', 'pdf',
                                                     'ps', 'eps', 'svg'))):
                plot_filename = args.plot_file+".pdf"

        if (plot_filename): matplotlib.use('Agg')

        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        matplotlib.rcParams['font.size'] = 10

        fig = plt.figure(figsize=(3.0, 3.0))
        if (args.plot_x_label):
            plt.xlabel(args.plot_x_label)
        else:
            plt.xlabel(args.plot_x_axis)

        if (args.plot_y_label):
            plt.ylabel(args.plot_y_label)
        else:
            plt.ylabel('Variables')

        if (len(args.plot_keys)):
            if (len(args.plot_keys) != len(variables)):
                raise KeyError("--plot-keys must be as long "
                               "as the number of variables.")
            else:
                plot_keys = dict(zip(variables, args.plot_keys))

        time_unit = args.dt

        for var in variables:

            cv = colvars_traj[var]

            x = cv.steps
            if (args.plot_x_axis == 'time'):
                x = cv.steps * time_unit
            elif (args.plot_x_axis == 'step'):
                x = cv.steps
            else:
                x = colvars_traj[args.plot_x_axis].values

            if (args.plot_time_factor > 1):
                # Usually 1
                start_ave = colvars_traj.num_frames % args.plot_time_factor
                x = x[start_ave:].reshape(x.shape[0]//args.plot_time_factor,
                                          args.plot_time_factor).mean(1)

            for ic in range(cv.num_dimensions):
                y = cv.values
                if (cv.num_dimensions > 1):
                    y = cv.values[:,ic]

                if (args.plot_time_factor > 1):
                    start_ave = colvars_traj.num_frames % args.plot_time_factor
                    y = y[start_ave:].reshape(y.shape[0]//args.plot_time_factor,
                                              args.plot_time_factor).mean(1)

                plt.plot(x, y, '-',
                         label=plot_keys[var],
                         alpha=0.5)

        plt.legend(loc='upper left')
        plt.tight_layout()

        if (plot_filename):
            plt.savefig(plot_filename)
        else:
            plt.show()

