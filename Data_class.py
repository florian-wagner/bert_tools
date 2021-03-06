#!/usr/bin/env python
"""
File: data_class.py
Author: Florian Wagner (fwagner@gfz-potsdam.de)
Description: Class for BERT data files
"""

import numpy as np
import pandas as pd

class Data:

    def __init__(self, filename, verbose=True):
        """
        Construct Data class instance from BERT's unified data file (*.ohm)
        """

        if verbose:
            print(("Reading in %s... \n" % filename))
        file = open(filename)

        self.eleccount = int(file.readline())
        elecs_str = file.readline()
        elecs_dim = len(elecs_str.rsplit()) - 1
        elecs_ix = elecs_str.rsplit()[1:]

        self.elecs = np.zeros((self.eleccount, elecs_dim), 'float')
        for i in range(self.eleccount):
            line = file.readline()
            self.elecs[i] = line.rsplit()

        self.datacount = int(file.readline())
        data_str = file.readline()
        data_dim = len(data_str.rsplit()) - 1
        data_ix = data_str.rsplit()[1:]

        self._string_ = """
        Number of electrodes: %s
        Dimension: %s
        Coordinates: %s
        Number of data points: %s
        Data header: %s
        """ % (self.eleccount, elecs_dim, elecs_str, self.datacount, data_str)

        self.data = np.zeros((self.datacount, data_dim), 'float')
        for i in range(self.datacount):
            line = file.readline()
            self.data[i] = line.rsplit()

        file.close()

        self.data = pd.DataFrame(self.data, columns=data_ix)
        self.elecs = pd.DataFrame(self.elecs, columns=elecs_ix)

        if verbose:
            print((self._string_))

    def save(self, filename):
        """ Save DataFrame to unified data format. """

        f = open(filename, 'w')
        f.write("%d\n" % self.eleccount)
        f.write("# ")
        for key in self.elecs.keys():
            f.write("%s " % key)
        f.write("\n")
        for row in self.elecs.itertuples(index=False):
            for val in row:
                f.write("%5.3f " % val)
            f.write("\n")
        f.write("%d\n" % self.datacount)
        f.write("# ")
        for key in self.data.keys():
            f.write("%s " % key)
        f.write("\n")
        for row in self.data.itertuples(index=False):
            for i, val in enumerate(row):
                if type(val) is np.float64:
                    if i < 4: # Account for ABMN TODO: make more elegant
                        f.write("%d " % val)
                    else:
                        f.write("%E " % val)
                else:
                    f.write("0.0 ")

            f.write("\n")
        f.close()

    def __str__(self):
        return self._string_
