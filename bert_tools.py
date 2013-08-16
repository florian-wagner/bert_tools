#!/usr/bin/env python
"""
BERT TOOLS
==========

Personal collection of useful functions for the usage with GIMLI/BERT.

Created on Thu Jun 28 11:49:15 2012
by fwagner@gfz-potsdam.de

"""

import numpy as np
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.sparse import coo_matrix


def abmn(n):
    """
       Construct all possible ABMN configurations for a given
       numer of sensors after Noel and Xu (1991)
    """
    combs = np.array(list(itertools.combinations(range(1,n+1),4)))
    perms = np.empty(((n*(n-3)*(n-2)*(n-1)/8),4),'int')
    for i in range(np.size(combs,0)):
        perms[0+i*3,:] = combs[i,:] # ABMN
        perms[1+i*3,:] = (combs[i,0],combs[i,3],combs[i,1],combs[i,2]) #AMNB
        perms[2+i*3,:] = (combs[i,0],combs[i,2],combs[i,1],combs[i,3]) #AMBN
    return perms

def area(a, b, c):
    """ Return area of triangle given the position vectors of corner points """
    area = 0.5 * np.linalg.norm(np.cross(b-a, c-a))
    return area

def das2ohm(input, output='data.ohm', verbose=True):
    """ Reads DAS-1 output and writes ohm-file for BERT """

    # Reading DAS-1 format
    if verbose: print 'Reading in', input, '... \n'
    file = open(input)

    elec_read, data_read = False, False
    elec, data = [], []
    known_tokens = {'ID': int, 'A': int, 'B': int, 'M': int, 'N': int,
                    'Appres': float, 'V/I,': float}

    for i, line in enumerate(file):
        if line.startswith('#'):
            if 'elec_start' in line: elec_read = True
            if 'elec_end' in line: elec_read = False
            if 'data_start' in line: data_read = True
            if 'data_end' in line: data_read = False
        elif line.startswith('!'):
            if data_read and 'ID' in line:
                tokens = line[1:].rsplit()
                found_tokens = {t: tokens.index(t) for t in known_tokens}
                r_token = found_tokens['V/I,']
                if tokens[r_token + 1] == 'Std.':
                    found_tokens['err'] = r_token + 1
        else:
            if elec_read:
                electrode = line.rsplit()
                id = [int(electrode[0].split(',')[-1])]
                xyz = map(float, electrode[1:4])
                elec.append(id + xyz)
            if data_read and len(line) > 180: # disregard erroneous data points:
                datum = line.rsplit()
                d = [None] * len(found_tokens)
                for token in found_tokens:
                    ix = found_tokens[token]
                    if token == 'err':
                        d[ix] = float(datum[ix])
                    elif known_tokens[token] == int:
                        d[ix] = int(datum[ix].split(',')[-1])
                    else:
                        d[ix] = known_tokens[token](datum[ix])
                data.append(d)

    file.close()

    #TODO: express error relative to resistance value

    das_tokens = ['A', 'B', 'M', 'N', 'V/I,', 'err']
    bert_tokens = ['a', 'b', 'm', 'n', 'r', 'err/Ohm']
    fmt = ['%d', '%d', '%d', '%d', '%e', '%e']

    if verbose:
        print '  Number of electrodes found:', len(elec)
        print '  Number of data found:', len(data)
        print '  Tokens found:', sorted(found_tokens.keys())
        print '  Tokens used:', das_tokens

    elec = np.asarray(elec)
    elec.sort(0)
    elec = elec[:,1:]

    # Writing BERT output
    file = open(output, 'w')
    file.write(str(len(elec)) + '\n')
    file.write('# %s %s %s \n' % ('x', 'y', 'z'))
    for pos in elec:
        for coord in pos:
            file.write('%.2f' % coord + '\t')
        file.write('\n')
    file.write(str(len(data)) + '\n')
    file.write('# %s %s %s %s %s %s \n' % tuple(bert_tokens))
    for d in data:
        for ix, t in enumerate(das_tokens):
            file.write(fmt[ix] % d[found_tokens[t]] + '\t')
        file.write('\n')

    file.close()
    if verbose: print '\nWritten data to %s.' % output

def describe(data):
    """ Print minimal statistic description of data """
    print "min:", np.min(data)
    print "mean:", np.mean(data)
    print "max:", np.max(data)

def intfile2mesh(file, mesh, method='cubic'):
    """
       Map point data from ASCII file to triangular mesh and return array with cell values.
    """
    data = np.loadtxt(file)

    return int2mesh(data, mesh, method=method)

def int2mesh(data, mesh, method='cubic'):
    """
       Map point data to triangular mesh and return array with cell values.
    """

    # extract cell centers

    cell_mids = np.zeros(( mesh.cellCount(), mesh.dim() ))
    for i, cell in enumerate( mesh.cells() ):
        cell_mids[ i ] = (cell.center()[ 0 ], cell.center()[ 1 ])

    return griddata(data[:,:2], data[:,2], cell_mids, fill_value=0, method=method)

def read_ohm(filename):
    """
       Read BERT data file (*.ohm) and return electrode positions and data array
    """
    print "Reading in %s... \n" % filename
    file = open(filename)

    elecs_count = int(file.readline().strip())
    elecs_str = file.readline()
    elecs_dim = len(elecs_str.rsplit()) - 1

    print "  Number of electrodes: %s" % elecs_count
    print "  Dimension: %s \n" % elecs_dim

    elecs_pos = np.zeros((elecs_count, elecs_dim), 'float')
    for i in range(elecs_count):
        line = file.readline()
        elecs_pos[i] = line.rsplit()

    data_count = int(file.readline().strip())
    data_str = file.readline()
    data_dim = len(data_str.rsplit()) - 1

    print "  Number of data points: %s" % data_count
    print "  Data header: %s" % data_str

    data = np.zeros((data_count, data_dim), 'float')
    for i in range(data_count):
        line = file.readline()
        data[i] = line.rsplit()

    file.close()

    return elecs_pos, data

def load_constraint_matrix(fname):
    """ Load constraint matrix in sparse format """
    global mesh
    i, j, data = np.loadtxt(fname, unpack=True, dtype=int)
    C = coo_matrix((data, (i,j)), shape=(mesh.boundaryCount(), mesh.cellCount()), dtype=int)
    return C

def loadsens(sensname):
    """
       Load sensitivity matrix from BERT binary file.
    """
    print "Loading %s... \n" % sensname
    fid = open(sensname,'rb')
    ndata = np.fromfile(fid, 'int32', 1); ndata = int(ndata[0])
    nmodel = np.fromfile(fid, 'int32', 1); nmodel = int(nmodel[0])
    S = np.empty((ndata,nmodel), 'float')

    for i in range(ndata):
        S[i,:] = np.fromfile(fid,'float',nmodel)
    print "  Number of ABMNs: %s" % ndata
    print "  Number of cells: %s \n" % nmodel
    print "%s loaded. (Size: %.2f GB)\n" % (sensname, 9.31323e-10 * S.nbytes)
    return S

def logdrop(data, lim=1e-3, normalize=False):
    """ Scale data logarithmically, maintain polarity, remove absolute values
        smaller lim and normalize everything to to the maximum value.
    """
 
    data = np.asarray(data)
    sign = np.sign(data)  # keep polarity
    data = np.abs(data) 
    data /= data.max()
    data /= lim
    data[data < 1] = 1
    data = np.log10(data)
    
    return data / data.max() * sign

def pdense(x, y, sigma, M=1000):
    """ Plot probability density of y with known stddev sigma
    """
    assert len(x) == len(y) and len(x) == len(sigma)
    N = len(x)
    # TODO: better y ranging
    ymin, ymax = min(y - 2 * sigma), max(y + 2 * sigma)
    yy = np.linspace(ymin, ymax, M)
    a = [np.exp(-((Y - yy) / s) ** 2) / s for Y, s in zip(y, sigma)]
    A = np.array(a)
    A = A.reshape(N, M)
    plt.imshow(-A.T, cmap='gray', aspect='auto',
               origin='lower', extent=(min(x), max(x), ymin, ymax))
    plt.title('Density plot')

def pole_pole(n, c=0, p=0, reciprocal=False, skip=None):
    """Return (reciprocal) complete pole-pole data set for n electrodes"""
    combs = list(itertools.combinations(range(1,n+1),2))
    confs = []
    for comb in combs:
        confs.append((comb[0], c, comb[1], p))
        if reciprocal:
            confs.append((comb[1], c, comb[0], p))
    print "%s configurations generated." % len(confs)

    confs = np.asarray(confs, dtype='int')

    if skip:
        idx = np.abs(confs[:,0] - confs[:,2]) > skip
        confs = confs[idx]
        print "%s configurations after skipping." % len(confs)

    return confs

def pole_bipole(n,c):
    """Return complete pole-bipole data set for n electrodes"""
    combs = list(itertools.combinations(range(1,n+1),2))
    confs = []
    for comb in combs:
        confs.append((comb[0], c, comb[1], p))
        if reciprocal:
            confs.append((comb[1], c, comb[0], p))
    print "%s configurations generated." % len(confs)
    #nicht fertig!! return np.asarray(confs, dtype='int')
       
def plotdata(mesh, data, cmap='Spectral_r', xlim=None, ylim=None, cmin=None,
             cmax=None, xlab='x (m)', ylab='depth (m)', clab='', title='',
             elecs=True, grid=True, bounds=False, orientation='vertical'):
    """
     Plot finite element data on triangular mesh
    """
    polys = []
    for cell in mesh.cells():
        if (cell.shape().nodeCount() == 3):
            polys.append(zip([cell.node(0).x(), cell.node(1).x(), cell.node(2).x()],
                               [cell.node(0).y(), cell.node(1).y(), cell.node(2).y()]))
        elif (cell.shape().nodeCount() == 4):
            polys.append(zip([cell.node(0).x(), cell.node(1).x(), cell.node(2).x(),
                                    cell.node(3).x()],
                               [cell.node(0).y(), cell.node(1).y(), cell.node(2).y(),
                                    cell.node(3).y()]))
        else:
            print "unknown shape to patch: " , cell.shape(), cell.shape().nodeCount()

    # Patch settings
    patches = mpl.collections.PolyCollection(polys)
    patches.set_antialiased(True)
    if grid:
        patches.set_linewidth(0.1)
        patches.set_edgecolor('0.5')
    else:
        patches.set_edgecolor('face')

    patches.set_array(data)
    patches.set_clim(cmin, cmax)
    patches.set_cmap(cmap)

    # Axes settings
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set_aspect('equal')
    axes.add_collection(patches)

    if xlim is None: xlim=(mesh.xmin(), mesh.xmax())
    if ylim is None: ylim=(mesh.ymin(), mesh.ymax())

    axes.set_xlim(xlim[0], xlim[1])
    axes.set_ylim(ylim[0], ylim[1])

    axes.set_xlabel(xlab)
    axes.set_ylabel(ylab)

    # Draw mesh boundaries
    if bounds:
        lines = []
        for bound in filter(lambda b: b.marker() > 1, mesh.boundaries()):
            lines.append(zip([bound.node(0).x(), bound.node(1).x()],
                             [bound.node(0).y(), bound.node(1).y()]))

        lineCollection = mpl.collections.LineCollection(lines)

        lineCollection.set_color('black')
        lineCollection.set_linewidth(1)
        lineCollection.set_linestyle('-')
        axes.add_collection(lineCollection)

    # Draw electrodes
    if elecs:
        el_cfg = mesh.findNodesIdxByMarker(-99)
        #print '%s electrodes found.' % len(el_cfg)
        el_pos = np.zeros([len(el_cfg), 2])
        for i, idx in enumerate(el_cfg):
            el_pos[i] = np.array([mesh.node(idx).x(), mesh.node(idx).y()])

        plt.plot(el_pos[:, 0], el_pos[:, 1], 'wo')

    # Colorbar settings
    cbar = fig.colorbar(patches, orientation=orientation)
    cbar.set_label('\n' + clab)

    plt.title(title + '\n')
    #return fig

def plotmesh(mesh, **kwargs):
    """ Plot mesh with cell attributes """
    data = np.asarray(mesh.cellAttributes())
    return plotdata(mesh, data, **kwargs)

def plotsens(mesh, data, cmap='RdBu_r', cmin=-1, cmax=1,
             clab='Normalized Sensitivity', abmn=False, cfg=False, **kwargs):
    """
     Plot sensitivities on triangular mesh
    """

    fig = plotdata(mesh, data, cmap, cmin=cmin, cmax=cmax, clab=clab, **kwargs)

    # Plot ABMNs
    if abmn:
        elecs, configs = read_ohm(abmn)
        configs = configs.astype('int')[:, :4]

        z = np.size(elecs, 1) - 1

        for k, e in enumerate(['A', 'B', 'M', 'N']):
            if elecs[:, 0][configs[cfg, k]-1] < -25: offset = -18; arrowpos = (1.1, 0.5)
            else: offset = 10; arrowpos = (-0.1, 0.5)
            plt.annotate(e, xy=(elecs[:, 0][configs[cfg, k]-1], elecs[:, z][configs[cfg, k]-1]), size=12,
                         va="center", bbox=dict(boxstyle="round", fc=(0.7, 0.7, 0.7), ec="black", alpha=0.75),
                         xytext=(offset, 0), textcoords='offset points',
                         arrowprops=dict(arrowstyle="wedge, tail_width=0.7",
                                        fc=(0.7, 0.7, 0.7), ec="black",
                                        patchA=None, relpos=arrowpos, alpha=0.75,
                                        ))
    return fig

def plothist(index,ohmfile,n=1000):
    """
    Plot histogramm of electrode usage (C or P) for best rated configurations
    """
    elecs, configs = read_ohm(ohmfile)
    configs = configs[:,:4].astype('int')

    num_elecs = len(elecs)
    el_no = np.arange(1,num_elecs + 1)

    hist_data = []
    for i in range(4):
        hist_data.append(np.bincount(configs[ index[:n], i],minlength=num_elecs + 1)[1:])

    # A+B & M+N
    C = hist_data[0] + hist_data[1]
    P = hist_data[2] + hist_data[3]

    fig = plt.figure()

    ax1 = plt.subplot(211)
    ax1.grid()
    ax1.axvspan(4.5,7.5, color='#669966', alpha=0.7)
    ax1.axvspan(23.2,26.5, color='#669966', alpha=0.7)
    ax1.axhline(C.mean(), color='#990033', linestyle='--', linewidth=1.5)
    ax1.bar(el_no - 0.4, C, edgecolor='k', color='#336699')
    ax1.set_ylabel('Current injection')
    ax1.set_axisbelow(True)

    #
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    ax2.grid()
    ax2.axvspan(4.5,7.5, color='#669966', alpha=0.7)
    ax2.axvspan(23.2,26.5, color='#669966', alpha=0.7)
    ax2.axhline(P.mean(), color='#990033', linestyle='--', linewidth=1.5)
    ax2.bar(el_no-0.4, P, color='#336699', edgecolor='k')
    ax2.set_xlim(0.5, 30.5)
    ax2.set_ylim(0, 420)
    ax2.set_xlabel('Electrode number')
    ax2.set_ylabel('Potential measurement')
    ax2.set_axisbelow(True)

    return fig

def write_configs(fname, elecs, configs):
    """
       Write abmn configurations to BERT data file (*.ohm)
    """
    print "Writing %s... \n" % fname
    fid = open(fname, 'w')
    fid.write('%d \n' % elecs.shape[0])

    if elecs.shape[1] == 2:
        fid.write('# x z \n')
    elif elecs.shape[1] == 3:
        fid.write('# x y z \n')
    else:
        print "WARNING: Wrong electrode dimension."

    for i in elecs:
        for j in i:
            fid.write('%.2f ' % j)
        fid.write('\n')

    fid.write('%d \n' % configs.shape[0])

    if configs.shape[1] > 4:
        print "\n WARNING: Specify additional attributes manually."

    fid.write('# a b m n \n')

    for i in configs:
        for j in i:
            if j in np.linspace(0, elecs.shape[0], elecs.shape[0]+1):
                fid.write('%d ' % j)
            else:
                fid.write('%f ' % j)
        fid.write('\n')

    fid.close()

def create2Dconfs(nel, ds=1):
    def create2Dint(nel, ds):
        c = 1  # counter
        ints = range(nel - (2 * ds) - 1)
        confs = []
        for i in ints:
            i += 1
            for j in range(ints[-i] + 1):
                j += 1
                a = i
                b = a + ds
                m = a + j + ds
                n = m + ds
                confs.append((a, b, m, n))
                c += 1

        return np.asarray(confs)

    if type(ds) == int:
        if ds > nel / 2:
            print "WARNING: dipole interval of %s is too large!" % ds
            confs = []
        else:
            confs = create2Dint(nel, ds)
    else:
        confs = []
        for i, d in enumerate(ds):
            if d > nel / 2:
                print "WARNING: dipole interval of %s is too large!" % d
                continue
            else:
                confs.append(create2Dint(nel, d))
        confs = np.vstack(confs)
    print "%s configurations generated." % len(confs)
    return confs

def create2Dxhconfs(nel, ds=1, inhole=False):
    def createxhint(nel, ds):

        bh1 = np.linspace(1, nel, nel, 'int')
        bh1 = bh1.astype('int')
        bh2 = np.linspace(nel + 1, 2 * nel, nel)
        bh2 = bh2.astype('int')

        confs = []

        for i in range(nel - ds):
            for j in range(nel - ds):
                a = bh1[i]
                b = a + ds
                m = bh2[j]
                n = m + ds
                confs.append((a, b, m, n))

        return np.asarray(confs)

    if type(ds) == int:
        if ds > nel - 1:
            print "WARNING: dipole interval of %s is too large!" % ds
            confs = []
        else:
            confs = createxhint(nel, ds)
    else:
        confs = []
        for i, d in enumerate(ds):
            if d > nel - 1:
                print "WARNING: dipole interval of %s is too large!" % d
                continue
            else:
                confs.append(createxhint(nel, d))
        confs = np.vstack(confs)
    if inhole:
        inconfs = create2Dconfs(nel, ds)
        confs = np.vstack((confs, inconfs, inconfs + nel))
    print "%s configurations generated." % len(confs)
    return confs
