# -*- coding: utf-8 -*-

"""
Module with functions related to laboratory work and data acquired with CATMAN.

"""

import numpy as np
import csv
import codecs


class Experiment:
    """
    Laboratory test data

    Class laboratory experiment containing methods for loading and manipulating data recorded with CATMAN software.

    """

    def __init__(self, header, data):
        self.header = header
        self.data = data

    @classmethod
    def from_file(cls, fh):
        """
        Method reading text files containing data recorded with CATMAN.

        Used to import data saved as ascii with CATMAN from the laboratory. ISO-8859-1 encoding is assumed.
        Warning: Columns in the file with the same name are overwritten, only the last one is added to the object.

        Parameters
        ----------
        fh : str
            File path
        """

        # Open the requested file.
        f = codecs.open(fh, 'r', 'ISO-8859-1')

        # Read the header
        header = list(csv.reader([next(f) for x in range(7)], delimiter='\t'))

        # Read the column headers.
        next(f)
        column_head = list(csv.reader([next(f) for x in range(29)], delimiter='\t'))

        # Read the tab separated values.
        next(f)
        values = list(csv.reader(f, delimiter='\t'))

        # Get the number of channels
        n_chan = len(values[0])

        # Create a dictionary.
        data = {}

        # Build a dictionary with the data using the column header to fetch the dict keys.
        for i in range(n_chan):
            channel = np.zeros((len(values), 1))
            name = column_head[0][i].partition(' ')[0]
            for j, row in enumerate(values):
                channel[j] = (float(row[i].replace(',', '.')))
            data[name] = channel

        # Create object
        return cls(header, data)
