# -*- coding: utf-8 -*-

"""
Module with functions related to laboratory work. Currently it contains tools related to 3D scanning and data
acquired with CATMAN.
"""

import numpy as np
import csv
import codecs
import scipy.linalg
from scipy import odr
from stl import mesh
import os
from steel_toolbox.steel_tools import eccentricity
import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Scan3D:
    """
    3D model data.

    Class of 3D objects. Can be imported from an .stl file of a .txt file of list of node coordinates.
    """

    def __init__(self, scanned_data):
        self.scanned_data = scanned_data

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

        return cls(scanned_data)

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

        return cls(scanned_data)

    @classmethod
    def from_pickle(cls, fh):
        """
        Method for importing a pickle file containing x, y, z, coordinates.

        Used to import data exported from blender. The pickle file is should contain a list of lists.

        """
        with open(fh, 'rb') as f:
            return cls(np.array(pickle.load(f)))

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

    def plot_surf(self):
        """
        Method plotting the model as a 3D surface.
        """
        # Create the x, y, z numpy arrays
        X = [i[0] for i in self.scanned_data]
        Y = [i[1] for i in self.scanned_data]
        Z = [i[2] for i in self.scanned_data]

        # Create a figure.
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Plot the data
        ax.plot_trisurf(X, Y, Z)


class PolygColumn:
    """
    A column specimen of polygonal cross-section.

    Used for the scanned polygonal specimens.
    """

    def __init__(self, sides=None, edges=None, centre_line=None):
        self.sides = sides
        self.edges = edges
        self.centre_line = centre_line
        if self.sides is None:
            self.sides = []

        if self.edges is None:
            self.edges = []

    def add_centre_line(self, points):
        """
        Calculate the centre axis of the column from 2 given points.
        :return:
        """
        self.centre_line = Line3D.from_2_points(points[0], points[1])

    def add_single_side_from_pickle(self, filename):
        self.sides.append(FlatFace.from_pickle(filename))

    def add_all_sides(self, n_sides, prefix, planar_fit=None, intersect_sides=None):
        if planar_fit is None:
            planar_fit = False
        if intersect_sides is None:
            intersect_sides = False

        self.sides = [FlatFace.from_pickle(prefix + '{:02d}.pkl'.format(x)) for x in range(1, n_sides + 1)]

        if planar_fit:
            [x.planar_fit(lay_on_xy=True) for x in self.sides]

        if planar_fit and intersect_sides:
            self.edges = [self.sides[x] & self.sides[x + 1] for x in range(-len(self.sides), 0)]

    def plot_all(self):
        max_z = max([x.scanned_data[:, 2].max() for x in self.sides])
        min_z = min([x.scanned_data[:, 2].min() for x in self.sides])
        fig1 = plt.figure()
        Axes3D(fig1)
        for i in range(-len(self.sides), 0):
            self.sides[i].plot_xy_bounded(reduced=0.003, fig=fig1)
            self.edges[i].plot_line(fig=fig1, ends=[min_z, max_z])

    def print_report(self):
        max_z = max([x.scanned_data[:, 2].max() for x in self.sides])
        min_z = min([x.scanned_data[:, 2].min() for x in self.sides])
        for i in range(len(self.sides)):
            print('Side {} is : {}'.format(i + 1, self.sides[i].plane_coeff))
            print('')
            print('Edge {} (sides {}-{})\n    Direction : {}\n    Through points : \n{}\n{}'.format(
                i + 1,
                i + 2,
                i + 1,
                self.edges[0].parallel,
                self.edges[0].xy_for_z(min_z),
                self.edges[0].xy_for_z(max_z))
            )
            print('')


class FlatFace(Scan3D):
    """
    Subclass of the Scan3D class, specifically for flat faces.

    Used for the individual faces of the polygonal specimens.
    """

    def __init__(self, scanned_data=None, plane_coeff=None):
        self.scanned_data = scanned_data
        self.plane_coeff = plane_coeff

        super().__init__(scanned_data)

    def __and__(self, other):
        """
        Object division returns the intersection line.

        :param other:
        :return:
        """

        if isinstance(other, FlatFace):
            # Calculate the parallel vector of the intersection line as the dot product of the vectors normal to the
            # planes.
            parallel = np.cross(np.r_[self.plane_coeff[0], self.plane_coeff[1], self.plane_coeff[2]],
                                np.r_[other.plane_coeff[0], other.plane_coeff[1], self.plane_coeff[2]])
            # Normalise the direction vector.
            parallel = unit_vector(parallel)

            # Calculate the intersection of the line with the xy plane
            a_1, a_2 = self.plane_coeff[0], other.plane_coeff[0]
            b_1, b_2 = self.plane_coeff[1], other.plane_coeff[1]
            c_1, c_2 = self.plane_coeff[2], other.plane_coeff[2]
            d_1, d_2 = self.plane_coeff[3], other.plane_coeff[3]
            # y_0 = (c_2 * a_1 - c_1 * a_2) / (b_1 * a_2 - b_2 * a_1)
            # x_0 = (-b_1 / a_1) / y_0 - c_1/a_1
            y_0 = (d_2 * a_1 - d_1 * a_2) / (b_1 * a_2 - b_2 * a_1)
            x_0 = (b_1 * y_0 + d_1) / (-a_1)
            z_0 = 0
            p_0 = np.array([x_0, y_0, z_0])

            return Line3D.from_point_and_parallel(p_0, parallel)
        else:
            return NotImplemented

    def planar_fit(self, lay_on_xy=False):

        """
        Fit a plane to 3d points.

        A regular least squares fit is performed (no error assumed in the given z-values).
        :return:
        """
        if lay_on_xy is None:
            lay_on_xy = False

        # Iterative fitting
        if lay_on_xy:
            # Perform least squares fit on the points "as is"
            beta1 = lstsq(self.scanned_data)

            # Z-axis unit vector.
            v1 = np.r_[0, 0, 1]

            # The normalised norm vector of the plane (which will be aligned to z axis)
            v2 = unit_vector(beta1[0:3])

            # Find the angle between the zz axis and the plane's normal vector, v2
            rot_ang = angle_between(v1, v2)

            # Find the rotation axis.
            rot_ax = unit_vector(np.r_[v2[1], -v2[0], 0])

            # Transform the points so that v2 is aligned to z.
            transformed = rotate_points(self.scanned_data, rot_ang, rot_ax)

            # Perform least squares.
            beta2 = lstsq(transformed)

            # Return the fitted plane to the original position of the points.
            beta2[:3] = rotate_points([beta2[:3]], -rot_ang, rot_ax)

            # Store the fitted plane in the instance.
            self.plane_coeff = beta2

        else:
            self.plane_coeff = lstsq(self.scanned_data)

    def quadratic_fit(self, scanned_data):
        """
        Fit a quadratic surface to 3d points.

        A regular least squares fit is performed (no error assumed in the given z-values).
        :param scanned_data:
        :return:
        """
        # best-fit quadratic curve
        a = np.c_[
            np.ones(self.scanned_data.shape[0]),
            self.scanned_data[:, :2],
            np.prod(self.scanned_data[:, :2], axis=1),
            self.scanned_data[:, :2] ** 2]

        self.plane_coeff, _, _, _ = scipy.linalg.lstsq(a, scanned_data[:, 2])

    def odr_planar_fit(self, rand_3_estimate=False):
        """
        Fit a plane to 3d points.

        Orthogonal distance regression is performed using the odrpack.

        :return:
        """

        def f_3(beta, xyz):
            """ implicit definition of the plane"""
            return beta[0] * xyz[0] + beta[1] * xyz[1] + beta[2] * xyz[2] + beta[3]

        # # Coordinates of the 2D points
        x = self.scanned_data[:, 0]
        y = self.scanned_data[:, 1]
        z = self.scanned_data[:, 2]
        # x = np.r_[9, 35, -13, 10, 23, 0]
        # y = np.r_[34, 10, 6, -14, 27, -10]
        # z = np.r_[100, 101, 101, 100, 101, 101]

        if rand_3_estimate:
            # initial guess for parameters
            # select 3 random points
            i = np.random.choice(len(x), size=3, replace=False)

            # Form the 3 points
            r_point_1 = np.r_[x[i[0]], y[i[0]], z[i[0]]]
            r_point_2 = np.r_[x[i[1]], y[i[1]], z[i[1]]]
            r_point_3 = np.r_[x[i[2]], y[i[2]], z[i[2]]]

            # Two vectors on the plane
            v_1 = r_point_1 - r_point_2
            v_2 = r_point_1 - r_point_3

            # normal to the 3-point-plane
            u_1 = np.cross(v_1, v_2)

            # Construct the first estimation, beta0
            d_0 = u_1[0] * r_point_1[0] + u_1[1] * r_point_1[1] + u_1[2] * r_point_1[2]
            beta0 = np.r_[u_1[0], u_1[1], u_1[2], d_0]
        else:
            self.planar_fit()
            beta0 = self.plane_coeff

        # Create the data object for the odr. The equation is given in the implicit form 'a*x + b*y + c*z + d = 0' and
        # beta=[a, b, c, d] (beta is the vector to be fitted). The positional argument y=1 means that the dimensionality
        # of the fitting is 1.
        lsc_data = odr.Data(np.row_stack([x, y, z]), y=1)
        # Create the odr model
        lsc_model = odr.Model(f_3, implicit=True)
        # Create the odr object based on the data, the model and the first estimation vector.
        lsc_odr = odr.ODR(lsc_data, lsc_model, beta0)
        # run the regression.
        lsc_out = lsc_odr.run()

        self.plane_coeff = lsc_out.beta / lsc_out.beta[3]
        # # Coefficients of the explicit plane formula 'z = A*x + B*y + C'. It stands alpha=beta/-c, where alfa=[A, B, C]
        # alpha = lsc_out.beta / (-lsc_out.beta[2])
        # self.plane_coeff = np.r_[alpha[0], alpha[1], alpha[3]]

    def centre_size(self):
        """Get the centre and the range of the data points."""

        # Check if plane data exists.
        if not isinstance(self.plane_coeff, np.ndarray):
            print('Wrong or missing plane coefficients')
            return NotImplemented

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

    def z_return(self, x, y):
        """
        Calculate z of a plane for given x, y.

        x : float or numpy.ndarray
        y : float or numpy.ndarray
        """
        if not isinstance(self.plane_coeff, np.ndarray):
            print('Wrong od missing plane coefficients')
            return NotImplemented
        alpha = (-self.plane_coeff / self.plane_coeff[2])
        z = alpha[0] * x + alpha[1] * y + alpha[3]
        return z

    def xy_return(self, z):
        """
        Intersect with a z=z0 plane.

        x : float or numpy.ndarray
        y : float or numpy.ndarray
        """
        return Line2D.from_line_coeff(
            self.plane_coeff[0],
            self.plane_coeff[1],
            self.plane_coeff[2] * z + self.plane_coeff[3]
        )

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
        z = self.z_return(x, y)

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

    def plot_z_bounded(self, z_bounds=None, fig=None, reduced=None):
        """
        Plot the plane between given upper and lower z limits.

        Intersecting the plane with a xy plane (for a given z) results in a line in two dimensions.

        :param fig:
        :param reduced:
        :return:
        """

        def quadratic(a, b, c):
            a = int(input("Enter the coefficients of a: "))
            b = int(input("Enter the coefficients of b: "))
            c = int(input("Enter the coefficients of c: "))

            d = b ** 2 - 4 * a * c  # discriminant

            if d < 0:
                x1 = None
                x2 = None
            elif d == 0:
                x1 = (-b + np.sqrt(b ** 2 - 4 * a * c)) / 2 * a
                x2 = x1
            else:
                x1 = (-b + np.sqrt((b ** 2) - (4 * (a * c)))) / (2 * a)
                x2 = (-b - np.sqrt((b ** 2) - (4 * (a * c)))) / (2 * a)

            return [x1, x2]

        # Get a figure to plot on
        if fig is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            ax = fig.get_axes()[0]

        if z_bounds is None:
            z_bounds = [-1, 1]

        # Make a randomly selected subset of points acc. to the input arg 'reduced=x'.
        if isinstance(reduced, float) and (0 < reduced < 1):
            i = np.random.choice(
                len(self.scanned_data[:, 0]),
                size=round(len(self.scanned_data[:, 0]) * reduced),
                replace=False
            )
        else:
            i = range(0, len(self.scanned_data[:, 0]))

        # Aaverage and range of the points.
        self.centre_size()

        # Calculate the radius of the circle as equal to the distance between the z limits.
        radious = abs(z_bounds[0] - z_bounds[1])

        # Line on 2d, intersection of the plane to the z=z_mid, where z_mid is the average z of the ziven points.
        l_z_mid = self.xy_return(self.centre[2])

        # Intersect the plane with a horizontal circle with centre at the average of the given points.
        # The coefficients of the 2nd order polynomial that solves the intersection of a circle to a line.
        a = self.plane_coeff[0]
        b = self.plane_coeff[1]
        c = self.plane_coeff[2]
        x_c = self.centre[0]
        y_c = self.centre[1]

        A = 1 + (a / b) ** 2
        B = 2 * (x_c - a * y_c / b - a * c / b)
        C = x_c ** 2 + y_c ** 2 - radius ** 2 + 2 * y_c * c / b + (c / b) ** 2

        roots = quadratic(A, B, C)


class Line3D:
    """A line in three dimensions"""

    def __init__(self, point=None, parallel=None):
        self.point = point
        self.parallel = parallel

    @classmethod
    def from_point_and_parallel(cls, point, parallel):
        """
        TODO: add check if the input data is numpy array and convert them to.
        :param point:
        :param parallel:
        :return:
        """
        # Normalise the given parallel vector
        parallel = unit_vector(np.r_[parallel])
        return cls(point=np.r_[point], parallel=parallel)

    @classmethod
    def from_2_points(cls, point1, point2):
        """
        TODO: add check if the input is np arrays and convert them to
        :param point1:
        :param point2:
        :return:
        """
        point1 = np.r_[point1]
        point2 = np.r_[point2]
        # Calculate and normalise the direction vector.
        parallel = unit_vector(point1 - point2)
        return cls(point=point1, parallel=parallel)

    @classmethod
    def from_pickle(cls, fh):
        """
        Import line from pickle.

        Used to import center lines for the polygonal specimens, as exported from blender.

        Parameters
        ----------
        fh: string
            Path and filename of the pickle file.
        """
        with open(fh, 'rb') as f:
            points = pickle.load(f)

        return cls.from_2_points(np.r_[points[0]], np.r_[points[1]])

    def xy_for_z(self, z_1):
        """Return x, y for a given z"""
        t = (z_1 - self.point[2]) / self.parallel[2]
        x_1 = self.parallel[0] * t + self.point[0]
        y_1 = self.parallel[1] * t + self.point[1]
        return np.r_[x_1, y_1, z_1]

    def plot_line(self, ends=None, fig=None):
        """
        Line segment plotter.

        Plot a segment of the line between two values of the parameter `t` in x=x0 + a*t

        ends : array like, optional
            The end values for the parametric form of the line segment to be plotted (array like with 2 values). Default
            is [-1, 1]
        fig : Object of class matplotlib.figure.Figure, optional
            The figure window to be used for plotting. By default, a new window is created.
        :return:
        """
        if ends is None:
            ends = np.array([-1, 1])

        x = self.point[0] + self.parallel[0] * np.r_[ends]
        y = self.point[1] + self.parallel[1] * np.r_[ends]
        z = self.point[2] + self.parallel[2] * np.r_[ends]

        if fig is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
        else:
            ax = fig.get_axes()[0]
        ax.plot(x, y, z, label='parametric curve')
        ax.legend()

        plt.show()


class Line2D:
    """A line in three dimensions"""

    def __init__(self, point=None, parallel=None, line_coeff=None):
        self.point = point
        self.parallel = parallel
        self.line_coeff = line_coeff

    @classmethod
    def from_point_and_parallel(cls, point, parallel):
        """
        TODO: add check if the input data is numpy array and convert them to.
        :param point:
        :param parallel:
        :return:
        """
        # Normalise the given parallel vector
        parallel = unit_vector(np.r_[parallel])
        line_coeff = np.r_[-parallel[1], parallel[0], parallel[1] * point[0] - parallel[0] * point[1]]
        return cls(point=np.r_[point], parallel=parallel, line_coeff=line_coeff)

    @classmethod
    def from_2_points(cls, point1, point2):
        """
        TODO: add check if the input is np arrays and convert them to
        :param point1:
        :param point2:
        :return:
        """
        point1 = np.r_[point1]
        point2 = np.r_[point2]
        # Calculate and normalise the direction vector.
        parallel = unit_vector(point1 - point2)
        return cls.from_point_and_parallel(point1, parallel)

    @classmethod
    def from_line_coeff(cls, A, B, C):
        parallel = np.r_[B, -A]
        point = [0, -(C / B)]
        line = cls.from_point_and_parallel(point, parallel)
        line.line_coeff = np.r_[A, B, C]
        return line

    @classmethod
    def from_pickle(cls, fh):
        """
        Import line from pickle.

        Used to import center lines for the polygonal specimens, as exported from blender.

        Parameters
        ----------
        fh: string
            Path and filename of the pickle file.
        """
        with open(fh, 'rb') as f:
            points = pickle.load(f)
            return cls.from_2_points(np.r_[points[0]], np.r_[points[1]])

    def x_for_y(self, y):
        """Return x a given y"""

        return (-self.line_coeff[1] * y - self.line_coeff[2]) / self.line_coeff[0]

    def y_for_x(self, x):
        """Return y a given x"""

        return (-self.line_coeff[0] * x - self.line_coeff[2]) / self.line_coeff[1]

    def plot_line(self, ends=None, fig=None):
        """
        Line segment plotter.

        Plot a segment of the line between two values of the parameter `t` in x=x0 + a*t

        ends : array like, optional
            The end values for the parametric form of the line segment to be plotted (array like with 2 values). Default
            is [-1, 1]
        fig : Object of class matplotlib.figure.Figure, optional
            The figure window to be used for plotting. By default, a new window is created.
        :return:
        """
        if ends is None:
            ends = np.array([-1, 1])

        x = self.point[0] + self.parallel[0] * np.r_[ends]
        y = self.point[1] + self.parallel[1] * np.r_[ends]

        if fig is None:
            fig = plt.figure()
            ax = fig.gca()
        else:
            ax = fig.get_axes()[0]
        ax.plot(x, y, label='parametric curve')
        ax.legend()

        plt.show()


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


def lstsq(points):
    # best-fit linear plane
    a = np.c_[points[:, 0], points[:, 1], np.ones(points.shape[0])]
    c, _, _, _ = scipy.linalg.lstsq(a, points[:, 2])  # coefficients

    # The coefficients are returned as an array beta=[a, b, c, d] from the implicit form 'a*x + b*y + c*z + d = 0'.
    # The vector is normalized so that [a, b, c] has a unit length and `d` is positive.
    return np.r_[c[0], c[1], -1, c[2]] / (np.linalg.norm([c[0], c[1], -1]) * np.sign(c[2]))


def rotate_points(points, rot_ang, rot_ax):
    """
    Rotate points for given angle around axis.

    :param points:
    :param rot_ang:
    :param rot_ax:
    :return:
    """
    # Rotation matrix
    sint = np.sin(rot_ang)
    cost = np.cos(rot_ang)
    ux, uy, uz = rot_ax

    rot_mtx = np.r_[
        [[cost + ux ** 2 * (1 - cost), ux * uy * (1 - cost) - uz * sint, ux * uz * (1 - cost) + uy * sint],
         [uy * ux * (1 - cost) + uz * sint, cost + uy ** 2 * (1 - cost), uy * uz * (1 - cost) - ux * sint],
         [uz * ux * (1 - cost) - uy * sint, uz * uy * (1 - cost) + ux * sint, cost + uz ** 2 * (1 - cost)]]
    ]

    # Transform the points.
    return np.array([np.dot(p, rot_mtx) for p in points])


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'"""

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def main():
    # Import the first two sides.
    sp1 = PolygColumn()
    sp1.add_all_sides(16, '../../sp1/sp1_side', planar_fit=True, intersect_sides=True)
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
