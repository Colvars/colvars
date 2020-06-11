import numpy as np


class colvars_grid:
    '''
    Manipulate gridded data from the Colvars Module

    The data lists may have one or several elements depending on the number of data series.
    PMF files contain one data series.
    Gradient files contain as many data series as there are dimensions

    # 3d surface plot for a PMF
    fig = plt.figure()
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(*pmf.meshgrid(), np.clip(pmf.data[0], None, 40))

    Args:
        filename (str): an optional filename to parse
    
    Attributes:
        filename (str): name of source file
        data (list of array): data for last frame
        histdata (list of array): data for all frames (array shape: nframes, nx[0], nx[1]...)
        dim (int): array dimension
        nx (list of int): data sizes
        nframes (int): n
    '''

    def __init__(self, filename=None):
        self.reset()
        if filename is not None:
            self.read(filename)


    def reset(self, ):
        '''Reset all data'''
        self.filenames = []
        self.data = []
        self.dim = 0
        self.nx = []
        self.xmin = []
        self.dx = []
        self.pbc = []
        self.histdata = []
        self.nframes = 0
        self.nsets = 0


    def summary(self):
        '''Print a short summary of the contents of the object'''
        print('Source files:', self.filenames)
        print('Grid dimension:', self.dim, ['PBC' if p else 'non-PBC' for p in self.pbc])
        print('Number of data series:', len(self.data))
        print('Number of grid points:', np.prod(self.nx), self.nx)
        print('Grid spacings:', self.dx)
        print('Number of time frames:', self.nframes)


    def axes(self):
        '''Returns the axes as a list of 1d meshes'''
        return [np.array([self.xmin[i] + (k+0.5)*self.dx[i] for k in range(self.nx[i])]) for i in range(self.dim)]


    def meshgrid(self, ):
        '''Returns a mesh grid suitable for plotting the data'''
        return np.meshgrid(*self.axes(), copy=True, indexing='ij')


    def interp_on_gradient_grid(self, data=None):
        '''Return the data array interpolated on the corresponding gradient grid.

        Valid for free energy surfaces.
        By default, acts on self.data[0]: the final state of the first data series.
        Can be used on another data array of the right shape, by passing the named
        parameter data.
        The gradient grid is shifted by half a bin width to the right, and one
        bin smaller in non-periodic dimensions.
        '''
        # Work on a copy of the data that we'll modify in place.
        if data is None:
            if len(self.data) == 1:
                data = self.data[0]
            else:
                print('Error: default is to interpolate scalar fields with a single data series, but this object has',
                len(self.data), 'data series. Pass the data parameter to customize.')
                return None

        new = np.copy(data)
        # Shift and interpolate linearly in each dimension successively
        for i, periodic in enumerate(self.pbc):
            if periodic:
                new = 0.5*(new + np.roll(new, shift=1, axis=i))
            else:
                sl1=[slice(None)]*self.dim
                sl1[i]=slice(1, None)
                sl2=[slice(None)]*self.dim
                sl2[i]=slice(0, -1)
                new = 0.5*(new[tuple(sl1)] + new[tuple(sl2)])
        return new


    def read(self, filename):
        '''Read data from a Colvars multicolumn file'''
        self.reset()
        self.filenames.append(filename)
        with open(filename) as f:
            l = f.readline().split()
            assert len(l) == 2
            self.dim = int(l[1])
            for _ in range(self.dim):
                l = f.readline().split()
                assert(len(l) == 5)
                self.xmin.append(float(l[1]))
                self.dx.append(float(l[2]))
                self.nx.append(int(l[3]))
                self.pbc.append(l[4] == "1")
            # Get number of columns
            f.readline() # empty line
            ncols = len(f.readline().split())
            f.close()

        self.nsets = ncols - self.dim
        self.histdata = [[] for _ in range(self.nsets)]
        self._append_data(filename)


    def append(self, filename):
        '''Append time frames from a Colvars multicolumn file with same data shape
        Can be used to load history files from consecutive runs.
        '''
        self.filenames.append(filename)
        with open(filename) as f:
            l = f.readline().split()
            assert len(l) == 2
            assert self.dim == int(l[1]), f'File to be appended has dimension {l[1]} instead of {self.dim}'
            xmin = []
            dx = []
            nx = []
            pbc = []
            for _ in range(self.dim):
                l = f.readline().split()
                assert(len(l) == 5)
                xmin.append(float(l[1]))
                dx.append(float(l[2]))
                nx.append(int(l[3]))
                pbc.append(l[4] == "1")

            assert xmin == self.xmin
            assert dx == self.dx
            assert nx == self.nx
            assert pbc == self.pbc

            # Get number of columns
            f.readline() # empty line
            ncols = len(f.readline().split())
            assert ncols == self.dim + self.nsets, f'File to be appended contains {ncols} columns, {self.dim + self.nsets} expected'
            f.close()
        self._append_data(filename)


    def _append_data(self, filename):
        '''Read data from given file, appending frames to the time series if already present'''
        # Timing: this routine took 6 seconds to read a file vs. 21 seconds for the previous version
        # that used np.loadtxt
        grid_size = np.prod(self.nx)
        rawdata = [[] for _ in range(self.nsets)]
        with open(filename) as f:
            for line in f:
                if line[0] == '#' or len(line) < 3:
                    continue
                l = line.split()
                # assume well-formed file for speed
                # Do not read the first self.dim columns which contain x values
                # Each data column goes into a flat list
                for i in range(self.nsets):
                    rawdata[i].append(float(l[self.dim + i]))

        # Compute number of time frames from length of first dataset (data column in multicol file)
        nf = len(rawdata[0]) / grid_size
        assert nf == int(nf), f'Data size {len(rawdata)} is not a multiple of grid size {grid_size}'
        nf = int(nf)
        self.nframes += nf

        self.data = [] # forget any previous references to non-final histdata

        # Recast data into histdata : list of (dataset) list of (time frame) array
        # data is a list of (dataset == latest time frame only) array
        for i in range(self.nsets):
            for t in range(nf):
                # Add time frame t to dataset i, converting to numpy array
                self.histdata[i].append(np.array(rawdata[i][grid_size*t:grid_size*(t+1)]))
            # Current data is ref to last frame of each data series
            self.data.append(self.histdata[i][-1])


    def _append_data_loadtxt(self, filename):
        '''An earlier version of _append_data using np.loadtxt - very slow for large datasets'''
        # Do not read the first self.dim columns which contain x values
        rawdata = np.loadtxt(filename, usecols=range(self.dim, self.dim + self.nsets), ndmin = 2)

        grid_size = np.prod(self.nx)
        nf = rawdata.shape[0] / grid_size
        assert nf == int(nf), f'Data size {rawdata.shape[0]} is not a multiple of grid size {grid_size}'
        nf = int(nf)
        self.nframes += nf

        self.data = [] # forget any previous references to non-final histdata

        # Recast data into histdata : list of (dataset) list of (time frame) array
        # data is a list of (dataset == latest time frame only) array
        for i in range(rawdata.shape[1]):
            for f in range(nf):
                # Slice np array into list of 1 array per frame
                # equivalent to a transpose().reshape() in numpy terms
                # but we have no need for huge arrays with all the data
                frame = rawdata[grid_size*f:grid_size*(f+1), i]
                # Store a clean copy
                self.histdata[i].append(np.copy(frame))
            # Current data is ref to last frame of each data series
            self.data.append(self.histdata[i][-1])


    def write(self, filename):
        '''Write data (final, not history) to a Colvars multicolumn file (2d only)'''
        if self.dim < 2 or self.dim > 3:
            print('write() only supports 2d and 3d grids at the moment')
            return
        with open(filename, 'w') as f:
            f.write(f'# {self.dim}\n')
            for i in range(self.dim):
                f.write('# {} {} {} {}\n'.format(self.xmin[i], self.dx[i],
                        self.nx[i], int(self.pbc[i])))
            f.write('')

            nx = self.nx[0]
            ny = self.nx[1]

            if self.dim == 2:
                for i in range(nx):
                    for j in range(ny):
                        x = self.xmin[0] + self.dx[0] * (i + .5)
                        y = self.xmin[1] + self.dx[1] * (j + .5)
                        f.write(f'{x} {y} ')
                        for d in self.data:
                            f.write(f'{d[i][j]} ')
                    f.write('\n')
            elif self.dim == 3:
                nz = self.nx[2]
                for i in range(nx):
                    for j in range(ny):
                        for k in range(nz):
                            x = self.xmin[0] + self.dx[0] * (i + .5)
                            y = self.xmin[1] + self.dx[1] * (j + .5)
                            z = self.xmin[2] + self.dx[2] * (k + .5)
                            f.write(f'{x} {y} {z} ')
                            for d in self.data:
                                f.write(f'{d[i][j][k]} ')
                            f.write('\n')
            f.close()

    @staticmethod
    def list2str(l):
        ''' Utility function to serialize a list into a string
        '''
        return ' '.join([str(i) for i in l])

    def write_dx(self, filename, dataset=None, frame=None):
        '''Write data to a dx file'''

        if dataset is None:
            if len(self.data) > 1:
                print(f'Dataset index must be specified for writing vector data to dx file.')
                return
            else:
                dataset = 0
        elif len(self.data) > 1:
                print(f'Writing dataset {dataset} of {len(self.data)}')


        if frame is None or frame == -1:
            # Last data point by default
            data = self.data[dataset]
            frame = self.nframes - 1
        else:
            data = self.histdata[frame][dataset]
        if self.nframes > 1:
            print(f'Writing frame {frame} of {self.nframes}')

        l2s = self.list2str
        with open(filename, 'w') as f:
            f.write(f'# DX file generated from {self.filename} by colvars_grid\n')
            f.write(f'object 1 class gridpositions counts {l2s(self.nx)}\n')
            f.write(f'origin {l2s(self.xmin)}\n')

            for i, dx in enumerate(self.dx):
                f.write(f'delta {l2s([0 for _ in range(i)])} {str(dx)} '
                        + f'{l2s([0 for _ in range(self.dim-1-i)])}\n')
            
            f.write(f'object 2 class gridconnections counts {l2s(self.nx)}\n')
            f.write(f'object 3 class array type double rank 0 items {np.prod(self.nx)} data follows\n')
    
            for i, d in enumerate(data.flat):
                # Write 10 items per line
                f.write(f'{d} ')
                if i % 10 == 9:
                    f.write('\n')
            f.close()

    # def plot(self):
    #     '''Plot the data...'''
    #     pass
    #     # TODO surfaces, 2d vector fields
    #     # choose correct rep by default (if n data == dim == 2, vector)
    #     # if n data = 1, contour in 2d, something in 3d


    def entropy(self):
        ''' Calculate the Shannon entropy of scalar data, interpreted as unnormalized distribution
            and the volume (number) of sampled bins
            If several frames are present, compute for all frames
        '''
        if len(self.data) > 1:
            print('Entropy is defined for a scalar dataset')
            return
        S_list = []
        V_list = []
        for data in self.histdata[0]:
            norm = data.sum()
            S = 0
            V = 0
            if norm > 0:
                # S will be reported as 0 if no sampling
                for d in data.flat:
                    if d > 0:
                        V += 1
                        S -= d / norm * np.log(d / norm)
            S_list.append(S)
            V_list.append(V)
        if len(S_list) == 1:
            return S_list[0], V_list[0]
        else:
            return S_list, V_list


    def convergence(self, ref=None, do_KL=False):
        '''Calculate the rmsd of data from a reference (by default, the final data)'''
        if ref is None:
            ref = self.data
        elif len(self.data) == 1 and type(ref) == np.ndarray and ref.shape == self.data[0].shape:
            # Accept single data series instead of singleton
            ref = [ref]
        elif len(ref) != len(self.data):
            print('Error: got', len(ref), 'reference data series, expected', len(self.data))
            return None
        # Adding the initial frame
        sum_sq = np.zeros(self.nframes + 1)
        #  Loop over data series
        for data_i, ref_i in zip(self.histdata, ref):
            # Initial time: our data is all zeros
            sum_sq[0] += np.sum(ref_i**2)
            # Loop over time
            for t, d_t in enumerate(data_i):
                sum_sq[t+1] += np.sum((d_t-ref_i)**2)

        npoints = len(self.data) * self.data[0].size
        rmsd = np.sqrt(sum_sq / npoints)

        if do_KL:
            # Convention: Boltzmann distribution at T=300K
            # Assume kcal/mol
            RT = 0.5962
            # Kullback-Leibler divergence
            KL = np.zeros(self.nframes + 1)
            if len(self.histdata) != 1:
                print('Can only compute KL for scalar data')
                return
            data = self.histdata[0]
            refdata = ref[0]
            refdist = np.exp(-1 * refdata / RT)
            refdist /= np.sum(refdist)
            # Initial time: our data is all zeros, dist = 1/N
            KL[0] = -1.0 / data.size * np.sum(np.log(refdist))
            # Loop over time
            for t, d_t in enumerate(data):
                dist_t = np.exp(-1 * d_t / RT)
                dist_t /= np.sum(dist_t)
                KL[t+1] = np.sum(dist_t * np.log(dist_t / refdist))
            return rmsd, KL 
        else:
            return rmsd

