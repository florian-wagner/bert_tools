#!/usr/bin/env python
"""
Sens_class
Created on Wed Sep 12 11:49:11 2012
by fwagner@gfz-potsdam.de
"""

import numpy as np
from progressbar import *
import scipy.spatial.distance as d
from scipy.stats import scoreatpercentile
from scipy import linalg as lin
import stopwatch as t
from bert_tools import *
import matplotlib.pyplot as plt

class Sens:
    """
    Sensitivity class
    =================

    Load and work with Jacobian array or binary matrix from BERT.
    """

    def __init__(self, sens, mesh=False, ohm=False):

        if isinstance(sens, str):
            self.sens = self._loadsens(sens)
        else:
            self.sens = sens
            self.ndata, self.ncells = self.sens.shape

        self.sum = np.sum(np.abs(self.sens), 0)
        if mesh: self.assign_mesh(mesh)
        if ohm: self.assign_ohmfile(ohm)
        self.size = 9.31323e-10 * self.sens.nbytes
        self.istransformed = None
        print self

    def _loadsens(self, fname):
        """
           Load sensitivity matrix from BERT binary file.
        """
        fid = open(fname, 'rb')
        ndata = np.fromfile(fid, 'int32', 1)
        ndata = int(ndata[0])
        ncells = np.fromfile(fid, 'int32', 1)
        ncells = int(ncells[0])
        S = np.empty((ndata, ncells), 'float')

        for i in range(ndata):
            S[i, :] = np.fromfile(fid, 'float', ncells)
        print fname + " loaded."
        self.ndata = ndata
        self.ncells = ncells
        return S

    def __str__(self):
        str = "Number of ABMNs: %s \nNumber of cells: %s \nSize: %.2f GB" % (self.ndata, self.ncells, self.size)
        return  str

    def assign_mesh(self, mesh):
        self.mesh = mesh

    def assign_ohmfile(self, fname):
        self.elecs, self.configs = read_ohm(fname)
    
    def check_rowsums(self):
        """ Check if row_sums are equal to 1 (Friedel, 2003) """
        if self.istransformed:
            print "WARNING: Jacobian is transformed."
            pass
        
        self.rowsums = self.sens.sum(axis=1)
        print "Min: %s" % self.rowsums.min()
        print "Max: %s" % self.rowsums.max()
        print "Mean: %s" % self.rowsums.mean()

    def compute_cumsum(self):
        """ Compute data-weighted cumulative sensitivity after Kemna (2000) """
        print "Computing data-weighted cumulative sensitivity..."
        t.start()
        J = np.asmatrix(self.sens)
        Wd = np.diag(np.ones(self.ndata, dtype='int'))
        self.cumsum = np.diag(J.T * Wd.T * Wd * J)
        t.stop()

    def compute_linind(self):
        """ Compute linear independency of sensitivities """
        print 'Computing linear independency...'
        t.start()
        length = len(self.sens)
        self.linind = np.zeros((length, length))
        loop = length - 1
        for i in range(loop):
            row = self.sens[i]
            for j in range(loop - i):
                secrow = self.sens[i + j + 1]
                self.linind[i, i + j + 1] = self.linind[i + j + 1, i] = np.linalg.norm(np.dot(row, secrow))/(np.linalg.norm(row) * np.linalg.norm(secrow))
                #self.linind[i, i + j + 1] = self.linind[i + j + 1, i] = d.cosine(row, secrow)
        t.stop()

    def compute_res(self, alpha=0.05):
        """ Compute formal model resolution matrix """
        print "Computing model resolution matrix..."
        t.start()
        C = np.diag(alpha * np.ones(self.ncells, dtype='int'))
        J = np.asmatrix(self.sens)
        self.R = np.diag(lin.inv(J.T * J + C) * J.T * J)
        t.stop()

    def normalize(self, perc):
        """
        Normalize sensitivity matrix and remove numerical effects (upper and lower perc)
        """
        if not hasattr(self, 'isnormalized'):
            print "Normalizing sensitivity matrix..."
            pbar = ProgressBar(widgets=[Percentage(), Bar()],
                               maxval=np.size(self.sens, 0)).start()
            for i in range(np.size(self.sens, 0)):
                lower_bound = scoreatpercentile(self.sens[i], perc)
                upper_bound = scoreatpercentile(self.sens[i], 100 - perc)
                self.sens[i][self.sens[i] < lower_bound] = lower_bound
                self.sens[i][self.sens[i] > upper_bound] = upper_bound
                cmax = np.max(np.abs(np.array((self.sens[i].min(), self.sens[i].max()))))
                cmin = -cmax
                self.sens[i] = np.divide(self.sens[i], cmax)
                pbar.update(i + 1)
            pbar.finish()

            self.isnormalized = True
        else:
            print "Sensitivity matrix already normalized."

    def plot(self, cfg, cmap='RdBu_r', cmin=-1, cmax=1,
             clab='Normalized Sensitivity', **kwargs):
        """
         Plot sensitivities on triangular mesh
        """

        fig = plotdata(self.mesh, self.sens[cfg], cmap, cmin=cmin, cmax=cmax, clab=clab, **kwargs)

        # Plot ABMNs
        if hasattr(self, 'elecs') and hasattr(self, 'configs'):
            self.configs = self.configs.astype('int')[:, :4]

            z = np.size(self.elecs, 1) - 1

            for k, e in enumerate(['A', 'B', 'M', 'N']):
                if self.elecs[:, 0][self.configs[cfg, k]-1] < -25: offset = -18; arrowpos = (1.1, 0.5)
                else: offset = 10; arrowpos = (-0.1, 0.5)
                plt.annotate(e, xy=(self.elecs[:, 0][self.configs[cfg, k]-1], self.elecs[:, z][self.configs[cfg, k]-1]), size=12,
                             va="center", bbox=dict(boxstyle="round", fc=(0.7, 0.7, 0.7), ec="black", alpha=0.75),
                             xytext=(offset, 0), textcoords='offset points',
                             arrowprops=dict(arrowstyle="wedge, tail_width=0.7",
                                            fc=(0.7, 0.7, 0.7), ec="black",
                                            patchA=None, relpos=arrowpos, alpha=0.75,
                                            ))
        else:
            print "No ohmfile assigned."
        return fig

    def plot_linind(self):
        fig = plt.imshow(self.linind)
        plt.colorbar()
        plt.gca().invert_yaxis()
        return fig

    def plot_res(self, cmap='Reds', **kwargs):
        fig = plotdata(self.mesh, self.R, cmap, cmin=0, cmax=1)
        return fig

    def write_configs(self, fname, cfgs):
        """
           Write abmn configurations to BERT data file (*.ohm)
        """
        configs = self.configs(cfgs)

        print "Writing %s... \n" % fname
        fid = open(fname, 'w')
        fid.write('%d \n' % self.elecs.shape[0])

        if self.elecs.shape[1] == 2:
            fid.write('# x z \n')
        elif self.elecs.shape[1] == 3:
            fid.write('# x y z \n')
        else:
            print "WARNING: Wrong electrode dimension."

        for i in self.elecs:
            for j in i:
                fid.write('%.2f ' % j)
            fid.write('\n')

        fid.write('%d \n' % configs.shape[0])

        if configs.shape[1] > 4:
            print "\n WARNING: Specify additional attributes manually."

        fid.write('# a b m n \n')

        for i in configs:
            for j in i:
                if j in np.linspace(1, self.elecs.shape[0], self.elecs.shape[0]):
                    fid.write('%d ' % j)
                else:
                    fid.write('%f ' % j)
            fid.write('\n')

        fid.close()
