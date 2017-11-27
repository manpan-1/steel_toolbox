# -*- coding: utf-8 -*-

"""
Module with functions related to laboratory work. Currently it contains tools related to 3D scanning and data
acquired with CATMAN.
"""

import numpy as np
import csv
import codecs
from stl import mesh
import os
from steel_toolbox.steel_tools import eccentricity


class Scan3D:
    """
    3D model data.

    Class of 3D objects. Can be imported from an .stl file of a .txt file of list of node coordinates.
    """

    def __init__(self, model):
        self.model = model

    # def trim_height(self, low, high):
    #     """
    #     Removes all entries from a numpy array containing points, for which the z coordinate is outside the given
    #     upper and lower limit. The numpy array is of the form [points, coordinates]z is expected at the 3rd position
    #     of each coordinate set.
    #     """
    #
    #     # Reduced array
    #     self.points = self.points[(self.points[:, 2] >= low) & (self.points[:, 2] <= high)]

    @classmethod
    def from_stl_file(cls, fh, del_original=None):
        """
        Import stl file.

        Alternative constructor, creates a Scan3D object by reading data from an .stl file. In case the file is created
        by Creaform's software (it is detected using name of the solid object as described in it's frist line), it is
        corrected accordingly before importing. The original file is renamed by adding '_old' before the extension or
        they can be deleted automatically if specified so.

        Parameters
        ----------
        fh : str
            File path.
        del_original : bool, optional
            Keep or delete the original file. Default is keep.
        """
        with open(fh, 'r') as f:
            fl = f.readlines(1)[0]
            identifier = fl.split(None, 1)[1]

        if identifier == 'ASCII STL file generated with VxScan by Creaform.\n':
            # Correct the file format
            Scan3D.repair_stl_file_structure(fh, del_original=del_original)

        model = mesh.Mesh.from_file(fh)

        return cls(model)

    @classmethod
    def from_coordinates_file(cls, fh):
        """
        Method reading text files containing x, y, z coordinates.

        Used to import data from 3D scanning files.
        """

        # Open the requested file.
        with open(fh, 'r') as f:
            # Number of points.
            n_of_points = len(f.readlines())

            # Initialise a numpy array for the values.
            all_coords = np.empty([n_of_points, 3])

            # Reset the file read cursor and loop over the lines of the file populating the numpy array.
            f.seek(0)
            for i, l in enumerate(f):
                all_coords[i] = l.split()

        return cls(all_coords)

    @staticmethod
    def repair_stl_file_structure(fh, del_original=None):
        """
        Repair header-footer of files created by Creaform's package.

        The .stl files created by Creaform's software are missing standard .stl header and footer. This method will
        create a copy of the requested file with proper header-footer using the filename (without the extension) as a
        name of the solid.

        Parameters
        ----------
        fh : str
            File path.
        del_original : bool, optional
            Keep or delete the original file. Default is keep.
        """
        if del_original is None:
            del_original = False
        solid_name = os.path.splitext(os.path.basename(fh))[0]

        start_line = "solid " + solid_name + "\n"
        end_line = "endsolid " + solid_name
        old_file = os.path.splitext(fh)[0]+'old.stl'

        os.rename(fh, old_file)
        with open(old_file) as fin:
            lines = fin.readlines()
        lines[0] = start_line
        lines.append(end_line)

        with open(fh, 'w') as fout:
            for line in lines:
                fout.write(line)

        if del_original:
            os.remove(old_file)


class Experiment:
    """
    Laboratory test data

    Class laboratory experiment containing methods for loading and manipulating data recorded with CATMAN software.
    """

    def __init__(self, header, data):
        self.header = header
        self.data = data

    def add_eccentricity(self, axis, column, moi, min_dist, thickness, young):
        """
        Calculate eccentricity.

        Adds a column in the data dictionary for the eccentricity of the load application on a given axis based on
        two opposite strain measurements.
        """

        self.data['e_' + axis] = []
        for load, strain1, strain2 in zip(self.data['Load'], self.data[column[0]], self.data[column[1]]):
            self.data['e_' + axis].append(eccentricity(
                load * 1000,
                [strain1 * 1e-6, strain2 * 1e-6],
                moi,
                min_dist + thickness / 2,
                young)
            )

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
