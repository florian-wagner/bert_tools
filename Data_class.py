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

        self.elecs = np.zeros((elecs_count, elecs_dim), 'float')
        for i in range(elecs_count):
            line = file.readline()
            self.elecs[i] = line.rsplit()

        data_count = int(file.readline())
        data_str = file.readline()
        data_dim = len(data_str.rsplit()) - 1
        data_ix = data_str.rsplit()[1:]
        
        self._string_ = """
        Number of electrodes: %s
        Dimension: %s
        Coordinates: %s 
        Number of data points: %s
        Data header: %s
        """ % (elecs_count, elecs_dim, elecs_str, data_count, data_str)
        
        self.data = np.zeros((data_count, data_dim), 'float')
        for i in range(data_count):
            line = file.readline()
            self.data[i] = line.rsplit()
        
        file.close()
        
        self.data = pd.DataFrame(self.data, columns=data_ix)
        self.elecs = pd.DataFrame(self.elecs, columns=elecs_ix)
        
        print self._string_
    
    def __str__(self):
        return self._string_
