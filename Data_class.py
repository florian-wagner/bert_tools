#!/usr/bin/env python
"""    
File: data_class.py
Author: Florian Wagner (fwagner@gfz-potsdam.de)
Description: Class for BERT data files
"""

import numpy as np
import pandas as pd

class Data:

    def __init__(self, filename):
        """
        Construct Data class instance from BERT's unified data file (*.ohm)
        """

        print "Reading in %s... \n" % filename
        file = open(filename)

        elecs_count = int(file.readline())
        elecs_str = file.readline()
        elecs_dim = len(elecs_str.rsplit()) - 1
        elecs_ix = elecs_str.rsplit()[1:]

        print "  Number of electrodes: %s" % elecs_count
        print "  Dimension: %s" % elecs_dim
        print "  Coordinates: %s \n" % elecs_str

        self.elecs = np.zeros((elecs_count, elecs_dim), 'float')
        for i in range(elecs_count):
            line = file.readline()
            self.elecs[i] = line.rsplit()

        data_count = int(file.readline())
        data_str = file.readline()
        data_dim = len(data_str.rsplit()) - 1
        data_ix = data_str.rsplit()[1:]
        
        print "  Number of data points: %s" % data_count
        print "  Data header: %s" % data_str

        self.data = np.zeros((data_count, data_dim), 'float')
        for i in range(data_count):
            line = file.readline()
            self.data[i] = line.rsplit()
        
        file.close()
        
        self.data = pd.DataFrame(self.data, columns=data_ix)
        self.elecs = pd.DataFrame(self.elecs, columns=elecs_ix)
