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
from steel_toolbox import analytic_geometry as ag
import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Scan3D:
    """
    3D model data.

    Class of 3D objects. Can be imported from an .stl file of a .txt file of list of node coordinates.
    """

    def __init__(self, scanned_data=None, grouped_data=None, centre=None, size=None):
        self.scanned_data = scanned_data
        self.grouped_data = grouped_data
        self.centre = centre
        self.size = size

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

        scanned_data = mesh.Mesh.from_file(fh)

        return cls(scanned_data=scanned_data)

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
            scanned_data = np.empty([n_of_points, 3])

            # Reset the file read cursor and loop over the lines of the file populating the numpy array.
            f.seek(0)
            for i, l in enumerate(f):
                scanned_data[i] = l.split()

        return cls(scanned_data=scanned_data)

    @classmethod
    def from_pickle(cls, fh):
        """
        Method for importing a pickle file containing x, y, z, coordinates.

        Used to import data exported from blender. The pickle file is should contain a list of lists.

        """
        with open(fh, 'rb') as f:
            return cls(scanned_data=np.array(pickle.load(f)))

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
        old_file = os.path.splitext(fh)[0] + 'old.stl'

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

    def sort_on_axis(self, axis=None):
        """
        Sort scanned data.

        The scanned points are sorted for a given axis.

        :param axis:
        :return:
        """
        if axis is None:
            axis = 0

        self.scanned_data = self.scanned_data[np.argsort(self.scanned_data[:, axis])]

    def quantize(self, axis=None, tolerance=None):
        """
        Group the scanned data.

        The points with difference on a given axis smaller than the tolerance are grouped together and stored in a list
        in the attribute `grouped_data`.

        :param axis:
        :param tolerance:
        :return:
        """
        if axis is None:
            axis = 0

        if tolerance is None:
            tolerance = 1e-4

        self.sort_on_axis(axis=axis)
        self.grouped_data = [np.r_[[self.scanned_data[0]]]]
        for point in self.scanned_data:
            if point[axis] - self.grouped_data[-1][0][axis] < tolerance:
                self.grouped_data[-1] = np.vstack([self.grouped_data[-1], point])
            else:
                self.grouped_data.append(np.r_[[point]])

    def centre_size(self):
        """
        Get the centre and the range of the data points.

        Used in combination with the plotting methods to define the bounding box.
        """
        # Bounding box of the points.
        x_min = min([i[0] for i in self.scanned_data])
        x_max = max([i[0] for i in self.scanned_data])
        y_min = min([i[1] for i in self.scanned_data])
        y_max = max([i[1] for i in self.scanned_data])
        z_min = min([i[2] for i in self.scanned_data])
        z_max = max([i[2] for i in self.scanned_data])
        x_range = abs(x_max - x_min)
        y_range = abs(y_max - y_min)
        z_range = abs(z_max - z_min)
        x_mid = (x_max + x_min) / 2
        y_mid = (y_max + y_min) / 2
        z_mid = (z_min + z_max) / 2

        self.centre = np.r_[x_mid, y_mid, z_mid]
        self.size = np.r_[x_range, y_range, z_range]

    def plot_surf(self):
        """
        Method plotting the model as a 3D surface.
        """
        # Create the x, y, z numpy arrays
        x = [i[0] for i in self.scanned_data]
        y = [i[1] for i in self.scanned_data]
        z = [i[2] for i in self.scanned_data]

        # Create a figure.
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Plot the data
        ax.plot_trisurf(x, y, z)


class PolygonalColumnSpecimen:
    """
    A column specimen of polygonal cross-section.

    Used for the scanned polygonal specimens.
    """

    def __init__(self, sides=None, edges=None, centre_line=None, thickness=None):
        if sides is None:
            sides = []

        if edges is None:
            edges = []

        self.sides = sides
        self.edges = edges
        self.centre_line = centre_line
        self.thickness = thickness

    def add_centre_line(self, point1, point2):
        """
        Calculate the centre axis of the column from 2 given points.

        :return:
        """
        self.centre_line = ag.Line3D.from_2_points(point1, point2)

    def add_single_side_from_pickle(self, filename):
        """
        Create a FlatFace instance as one side af the polygon column.

        The FlatFace instance is created from a pickle file of scanned data points.

        :param filename:
        :return:
        """
        self.sides.append(FlatFace.from_pickle(filename))

    def add_all_sides(self, n_sides, prefix, fit_planes=None, offset_to_midline=False):
        """
        Add multiple sides.

        Multiple FlatFace instances are created as sides of the polygonal column. A series of files containing scanned
        data points must be given. The files should be on the same path and have a filename structure as:
        `path/basenameXX.pkl`, where XX is an id number in ascending order starting from 01.
        Only the `path/filename` is given as input to this method.

        :param n_sides:
        :param prefix:
        :param fit_planes:
        :return:
        """
        if fit_planes is None:
            fit_planes = False

        self.sides = [FlatFace.from_pickle(prefix + '{:02d}.pkl'.format(x)) for x in range(1, n_sides + 1)]

        if fit_planes:
            [x.fit_plane() for x in self.sides]
        if offset_to_midline:
            offset = self.thickness / 2
            [x.offset_face(offset, offset_points=True) for x in self.sides]

    def add_single_edge_from_pickle(self, filename):
        """
        Create a RoundEdge instance as one edges af the polygon column.

        The RoundEdge instance is created from a pickle file of scanned data points.

        :param filename:
        :return:
        """
        self.edges.append(RoundedEdge.from_pickle(filename))

    def add_all_edges(self, n_sides, prefix, ref_lines=False):
        """
        Add multiple edges.

        Multiple RoundEdge instances are created as edges of the polygonal column. A series of files containing scanned
        data points must be given. The files should be on the same path and have a filename structure as:
        `path/basenameXX.pkl`, where XX is an id number in ascending order starting from 01.
        Only the `path/filename` is given as input to this method.

        After adding the sequential edges, if ref_lines=True, the reference lines are calculated as the intersections
        of sequential sides.

        :param n_sides:
        :param prefix:
        :param ref_lines:
        :return:
        """
        self.edges = [RoundedEdge.from_pickle(prefix + '{:02d}.pkl'.format(x)) for x in range(1, n_sides + 1)]

        if ref_lines:
            for x in range(-len(self.sides), 0):
                self.edges[x].add_ref_line(self.sides[x].ref_plane & self.sides[x + 1].ref_plane)

    def find_real_edges(self, offset_to_midline=False):
        """
        Find edge points on the scanned rounded edge.

        A series of points is returned which represent the real edge of the polygonal column. Each point is calculated
        as  the intersection of a circle and a line at different heights of the column, where the circle is best fit to
        the rounded edge scanned points and the line passing through the reference edge (see `add_all_edges`
        documentation) and the polygon's centre line.

        :return:
        """
        if offset_to_midline:
            offset = -self.thickness / 2
        else:
            offset = 0

        if isinstance(self.centre_line, ag.Line3D) and isinstance(self.edges, list):
            for x in self.edges:
                x.fit_circles(axis=2, offset=offset)
                x.intersect_data(self.centre_line)
        else:
            NotImplemented

    def plot_all(self):
        """
        Plot all data.

        :return:
        """
        max_z = max([x.scanned_data[:, 2].max() for x in self.sides])
        min_z = min([x.scanned_data[:, 2].min() for x in self.sides])
        fig1 = plt.figure()
        Axes3D(fig1)
        for i in range(-len(self.sides), 0):
            self.sides[i].plot_xy_bounded(reduced=0.003, fig=fig1)
            self.edges[i].ref_line.plot_line(fig=fig1, ends=[min_z, max_z])

    def print_report(self):
        """
        Print a report for the polygon column.

        :return:
        """
        max_z = max([x.scanned_data[:, 2].max() for x in self.sides])
        min_z = min([x.scanned_data[:, 2].min() for x in self.sides])
        for i in range(len(self.sides)):
            print('Side {} is : {}'.format(i + 1, self.sides[i].ref_plane.plane_coeff))
            print('')
            print('Edge {} (sides {}-{})\n    Direction : {}\n    Through points : \n{}\n{}'.format(
                i + 1,
                i + 1,
                i + 2,
                self.edges[i].ref_line.parallel,
                self.edges[i].ref_line.xy_for_z(min_z),
                self.edges[i].ref_line.xy_for_z(max_z))
            )
            print('')


class FlatFace(Scan3D):
    """
    Subclass of the Scan3D class, specifically for flat faces.

    Used for the individual faces of the polygonal specimens.
    """

    def __init__(self, scanned_data=None, ref_plane=None):
        self.ref_plane = ref_plane

        super().__init__(scanned_data)

    def fit_plane(self):
        self.ref_plane = ag.Plane3D.from_fitting(self.scanned_data, lay_on_xy=True)

    def offset_face(self, offset, offset_points=False):
        """
        Offset the plane and (optionally) the scanned data points.

        Useful for translating translating the scanned surface to the mid line.

        :param offset:
        :param offset_points:
        :return:
        """
        self.ref_plane.offset_plane(offset)

        if offset_points:
            point_list = [p + self.ref_plane.plane_coeff[:3] * offset for p in self.scanned_data]
            self.scanned_data = np.array(point_list)

    def plot_xy_bounded(self, fig=None, reduced=None):
        """
        Surface plotter.

        Plot the 3d points and the fitted plane.

        Parameters
        ----------
        fig : Object of class matplotlib.figure.Figure, optional
            The figure window to be used for plotting. By default, a new window is created.
        reduced: float, optional
            A reduced randomly selected subset of points is plotted (in case the data is too dense for plotting). The
            rediced size is given as a ratio of the total number of points, e.g `reduced=0.5` plots half the points. By
            default, all points are plotted.
        """

        # Average and range of the points.
        self.centre_size()

        x_lims = [self.centre[0] - self.size[0] / 2., self.centre[0] + self.size[0] / 2.]
        y_lims = [self.centre[1] - self.size[1] / 2., self.centre[1] + self.size[1] / 2.]

        x, y = np.meshgrid(x_lims, y_lims)

        # evaluate the plane function on the grid.
        z = self.ref_plane.z_return(x, y)

        # or expressed using matrix/vector product
        # z = np.dot(np.c_[xx, yy, np.ones(xx.shape)], self.plane_coeff).reshape(x.shape)
        plot_dim = max(self.size[0], self.size[1], self.size[2])

        # Quadratic:evaluate it on a grid
        # z = np.dot(np.c_[np.ones(xx.shape), xx, yy, xx * yy, xx ** 2, yy ** 2], self.plane_coeff).reshape(x.shape)

        # Get a figure to plot on
        if fig is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            ax = fig.get_axes()[0]

        # Plot the plane
        ax.plot_surface(x, y, z, rstride=1, cstride=1, alpha=0.2)

        # Make a randomly selected subset of points acc. to the input arg 'reduced=x'.
        if isinstance(reduced, float) and (0 < reduced < 1):
            i = np.random.choice(
                len(self.scanned_data[:, 0]),
                size=round(len(self.scanned_data[:, 0]) * reduced),
                replace=False
            )
        else:
            i = range(0, len(self.scanned_data[:, 0]))

        # Plot the points
        ax.scatter(self.scanned_data[i, 0], self.scanned_data[i, 1], self.scanned_data[i, 2], c='r', s=50)
        plt.xlabel('x')
        plt.ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim3d(self.centre[0] - plot_dim / 2, self.centre[0] + plot_dim / 2)
        ax.set_ylim3d(self.centre[1] - plot_dim / 2, self.centre[1] + plot_dim / 2)
        ax.set_zlim3d(self.centre[2] - plot_dim / 2, self.centre[2] + plot_dim / 2)

        # ax.axis('tight')
        plt.show()

        # Return the figure handle.
        return fig


class RoundedEdge(Scan3D):
    """
    A scanned rounded edge.
    """
    def __init__(self, scanned_data=None, ref_line=None, edge_points=None, circles=None):
        self.ref_line = ref_line
        self.edge_points = edge_points
        self.circles = circles

        super().__init__(scanned_data)

    def intersect_data(self, other):
        """
        Intersect scanned points with a surface between the reference line and a given line.

        This function is used to calculate the real edge based on the scanned rounded corner. Circles are fitted on the
        scanned points on different positions. Then the circles are intersected with the line passing through the
        reference line of the edge and another given line (e.g. the centre of the column). A list of points is
        generated which represent the real edge of rounded corner.

        :param other:
        :return:
        """
        if isinstance(other, ag.Line3D):
            self.edge_points = []
            for x in self.circles:
                z_current = x.points[0][0]
                ref_line_point = self.ref_line.xy_for_z(z_current)
                other_line_point = other.xy_for_z(z_current)
                intersection_line = ag.Line2D.from_2_points(ref_line_point[:2], other_line_point[:2])
                line_circle_intersection = x.intersect_with_line(intersection_line)
                if np.linalg.norm(line_circle_intersection[0]) > np.linalg.norm(line_circle_intersection[1]):
                    outer = line_circle_intersection[0]
                else:
                    outer = line_circle_intersection[1]
                self.edge_points.append(outer)
        else:
            NotImplemented

    def add_ref_line(self, ref_line):
        """
        Add a reference line for the edge.

        Useful when the rounded edge lies between flat faces and the theoretical edge is at their intersection.

        :param ref_line:
        :return:
        """
        if isinstance(ref_line, ag.Line3D):
            self.ref_line = ref_line
        else:
            print("ref_line must be Line3D")
            NotImplemented

    def fit_circles(self, axis=None, offset=None):
        """
        Fit a series of circles along the length of the rounded edge.

        The scanned data are first grouped together based on their z-coordinate ant then a horizontal circle is fitted
        for each group of points.

        :return:
        """
        if axis is None:
            axis = 0

        if offset is None:
            offset = 0

        self.quantize(axis=axis)
        self.circles = []
        for x in self.grouped_data:
            self.circles.append(ag.Circle2D.from_fitting(x))
            self.circles[-1].radius = self.circles[-1].radius + offset


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


def main():
    # Create a polygon column instance.
    sp1 = PolygonalColumnSpecimen(thickness=3)

    # Add a center line for the specimen.
    sp1.add_centre_line([0, 0, 0], [0, 0, 1])

    # Add all sides and edges.
    # they consist of FlatFace and RoundedEdge instances.
    sp1.add_all_sides(16, '../../sp1/sp1_side_', fit_planes=True, offset_to_midline=True)
    sp1.add_all_edges(16, '../../sp1/sp1_edge_', ref_lines=True)

    # Find a series of points for each edge based on the scanned surface.
    sp1.find_real_edges(offset_to_midline=True)

    #
    sp1.plot_all()
    sp1.print_report()

    # The following code is used to cross-check the validity of the results. Points for 2 planes are created and used
    # for testing the fitting and plotting methods.

    # Test xyz data from randomised z=1x+2y+3 for x,y values from -10 t0 +10
    def f_3(beta, xy):
        """ implicit definition of the plane"""
        return beta[0] * xy[0] + beta[1] * xy[1] + beta[2]

    def make_plane_data(beta):
        x = [x + np.random.rand() for x in np.linspace(-10, 10, 21)] + [x + np.random.rand() for x in
                                                                        np.linspace(-10, 10, 21)] + [
                x + np.random.rand() for x in np.linspace(-10, 10, 21)]
        y = [x + np.random.rand() for x in np.linspace(-10, 10, 21)] + [x + np.random.rand() for x in
                                                                        np.linspace(0, 20, 21)] + [x + np.random.rand()
                                                                                                   for x in
                                                                                                   np.linspace(10, 30,
                                                                                                               21)]
        z = f_3(beta, np.row_stack([x, y]))
        x = np.r_[x]
        y = np.r_[y]
        z = np.r_[z]
        return FlatFace(np.transpose(np.row_stack([x, y, z])))

    # Create points for the two planes.
    p1 = make_plane_data([1, 0, -4])
    p2 = make_plane_data([0, 1, -4])

    # Fit a plane on the created points.
    p1.odr_planar_fit()
    p2.odr_planar_fit()
    lp12 = p1 & p2

    # Plot the results.
    fig2 = plt.figure()
    Axes3D(fig2)
    p1.plot_xy_bounded(fig=fig2)
    p2.plot_xy_bounded(fig=fig2)
    lp12.plot_line(fig=fig2, ends=[-10, 10])

    # Return the specimen
    return sp1

# if __name__ == "__main__":
# main()
