
import numpy as np
import scipy.linalg
from scipy import odr
import pickle
from matplotlib import pyplot as plt

class Plane3D:
    def __init__(self, plane_coeff=None):
        self.plane_coeff = plane_coeff

    def __and__(self, other):
        """
        Object division returns the intersection line.

        :param other:
        :return:
        """

        if isinstance(other, Plane3D):
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

    def offset_plane(self, offset):
        """
        Offset the plane and (optionally) the scanned data points.

        Useful for translating translating the scanned surface to the mid line.

        :param offset:
        :param offset_points:
        :return:
        """
        self.plane_coeff = np.append(self.plane_coeff[:3], self.plane_coeff[3]-offset)

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

    @classmethod
    def from_fitting(cls, points, lay_on_xy=False):
        plane_coeff = planar_fit(points, lay_on_xy=lay_on_xy)
        return cls(plane_coeff=plane_coeff)

class Circle2D:
    """A circle in two dimensions."""

    def __init__(self, radius=None, centre=None, points=None):
        self.radius = radius
        self.centre = centre
        self.points = points

    def intersect_with_line(self, line):
        """Intersect circle with line"""
        if isinstance(line, Line2D):
            a, b, c = line.line_coeff
            xc, yc = self.centre
            radius = self.radius

            alfa = 1 + (a / b)**2
            beta = 2 * (a * c / b**2 - xc + a * yc / b)
            gama = xc**2 + yc**2 - radius**2 + 2 * c * yc / b

            x_intersect = solve_quadratic(alfa, beta, gama)
            y_intersect = np.r_[-(c + a * x_intersect[0]) / b, -(c + a * x_intersect[1]) / b]

            point1 = np.r_[x_intersect[0], y_intersect[0]]
            point2 = np.r_[x_intersect[1], y_intersect[1]]
            return [point1, point2]
        else:
            NotImplemented

    @classmethod
    def from_fitting(cls, points):
        xc, yc, rad = fit_circle(points)
        return cls(radius=rad, centre=np.r_[xc, yc], points=points)

    def plot_circle(self):
        """ Draw data points, best fit circles and center for the three methods,
        and adds the iso contours corresponding to the fiel residu or residu2
        """

        plt.figure(facecolor='white')  # figsize=(7, 5.4), dpi=72,
        plt.axis('equal')

        theta_fit = np.linspace(-np.pi, np.pi, 180)

        x_fit3 = self.centre[0] + self.radius * np.cos(theta_fit)
        y_fit3 = self.centre[1] + self.radius * np.sin(theta_fit)
        plt.plot(x_fit3, y_fit3, 'r-.', label='odr fit', lw=2)

        plt.plot([self.centre[0]], [self.centre[1]], 'kD', mec='w', mew=1)

        # draw
        plt.xlabel('x')
        plt.ylabel('y')

        plt.draw()
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()

        vmin = min(xmin, ymin)
        vmax = max(xmax, ymax)

        # plot data
        plt.plot(self.points[:, 0], self.points[:, 1], 'ro', label='data', ms=8, mec='b', mew=1)
        plt.legend(loc='best', labelspacing=0.1)

        plt.xlim(xmin=vmin, xmax=vmax)
        plt.ylim(ymin=vmin, ymax=vmax)

        plt.grid()
        plt.title('Least Squares Circle')


class Line3D:
    """A line in three dimensions."""

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
    def from_line_coeff(cls, alfa, beta, gama):
        parallel = np.r_[beta, -alfa]
        point = [0, -(gama / beta)]
        line = cls.from_point_and_parallel(point, parallel)
        line.line_coeff = np.r_[alfa, beta, gama]
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


def lstsq(points):
    # best-fit linear plane
    a = np.c_[points[:, 0], points[:, 1], np.ones(points.shape[0])]
    c, _, _, _ = scipy.linalg.lstsq(a, points[:, 2])  # coefficients

    # The coefficients are returned as an array beta=[a, b, c, d] from the implicit form 'a*x + b*y + c*z + d = 0'.
    # The vector is normalized so that [a, b, c] has a unit length and `d` is positive.
    return np.r_[c[0], c[1], -1, c[2]] / (np.linalg.norm([c[0], c[1], -1]) * np.sign(c[2]))


def planar_fit(points, lay_on_xy=False):

    """
    Fit a plane to 3d points.

    A regular least squares fit is performed. If the argument lay_on_xy is given True, a regular least squares
    fitting is performed and the result is used to rotate the points so that they come parallel to the xy-plane.
    Then, a second least squares fit is performed and the result is rotated back to the initial position. This
    procedure helps overcome the problem of vertical planes.

    :return:
    """
    if lay_on_xy is None:
        lay_on_xy = False

    # Iterative fitting
    if lay_on_xy:
        # Perform least squares fit on the points "as is"
        beta1 = lstsq(points)

        # Z-axis unit vector.
        v1 = np.r_[0, 0, 1]

        # The normalised norm vector of the plane (which will be aligned to z axis)
        v2 = unit_vector(beta1[0:3])

        # Find the angle between the zz axis and the plane's normal vector, v2
        rot_ang = angle_between(v1, v2)

        # Find the rotation axis.
        rot_ax = unit_vector(np.r_[v2[1], -v2[0], 0])

        # Transform the points so that v2 is aligned to z.
        transformed = rotate_points(points, rot_ang, rot_ax)

        # Perform least squares.
        beta2 = lstsq(transformed)

        # Return the fitted plane to the original position of the points.
        beta2[:3] = rotate_points([beta2[:3]], -rot_ang, rot_ax)

        # Store the fitted plane in the instance.
        return beta2

    else:
        # Perform a least squares fitting directly on the original position and return the coefficients.
        return lstsq(points)


def quadratic_fit(points):
    """
    Fit a quadratic surface to 3d points.

    A regular least squares fit is performed (no error assumed in the given z-values).
    :param scanned_data:
    :return:
    """
    # best-fit quadratic curve
    a = np.c_[
        np.ones(points.shape[0]),
        points[:, :2],
        np.prod(points[:, :2], axis=1),
        points[:, :2] ** 2]

    beta, _, _, _ = scipy.linalg.lstsq(a, points[:, 2])
    return beta


def odr_planar_fit(points, rand_3_estimate=False):
    """
    Fit a plane to 3d points.

    Orthogonal distance regression is performed using the odrpack.

    :return:
    """

    def f_3(beta, xyz):
        """ implicit definition of the plane"""
        return beta[0] * xyz[0] + beta[1] * xyz[1] + beta[2] * xyz[2] + beta[3]

    # # Coordinates of the 2D points
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
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
        beta0 = planar_fit()

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

    return lsc_out.beta / lsc_out.beta[3]


def fit_circle(points):
    """
    Fit a circle to a set of 2D points.

    The fitting is performed using the ODR from scipy. The code is partly taken from the scipy Cookbook.

    https://github.com/mpastell/SciPy-CookBook/blob/master/originals/Least_Squares_Circle_attachments/least_squares_circle_v1d.py

    :param points:
    :return:
    """
    x, y = points[:, 0], points[:, 1]

    def calc_r(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x - xc) ** 2 + (y - yc) ** 2)

    def f_3b(beta, var):
        """ implicit definition of the circle """
        return (var[0] - beta[0]) ** 2 + (var[1] - beta[1]) ** 2 - beta[2] ** 2

    def jacb(beta, var):
        """ Jacobian function with respect to the parameters beta.
        return df_3b/dbeta
        """
        xc, yc, r = beta
        xi, yi = var

        df_db = np.empty((beta.size, var.shape[1]))
        df_db[0] = 2 * (xc - xi)  # d_f/dxc
        df_db[1] = 2 * (yc - yi)  # d_f/dyc
        df_db[2] = -2 * r  # d_f/dr

        return df_db

    def jacd(beta, var):
        """ Jacobian function with respect to the input x.
        return df_3b/dx
        """
        xc, yc, r = beta
        xi, yi = var

        df_dx = np.empty_like(var)
        df_dx[0] = 2 * (xi - xc)  # d_f/dxi
        df_dx[1] = 2 * (yi - yc)  # d_f/dyi

        return df_dx

    def calc_estimate(data):
        """ Return a first estimation on the parameter from the data  """
        xc0, yc0 = data.x.mean(axis=1)
        r0 = np.sqrt((data.x[0] - xc0) ** 2 + (data.x[1] - yc0) ** 2).mean()
        return xc0, yc0, r0

    # for implicit function :
    #       data.x contains both coordinates of the points
    #       data.y is the dimensionality of the response
    lsc_data = odr.Data(np.row_stack([x, y]), y=1)
    lsc_model = odr.Model(f_3b, implicit=True, estimate=calc_estimate, fjacd=jacd, fjacb=jacb)
    lsc_odr = odr.ODR(lsc_data, lsc_model)  # beta0 has been replaced by an estimate function
    lsc_odr.set_job(deriv=3)  # use user derivatives function without checking
    lsc_out = lsc_odr.run()

    xc_odr, yc_odr, R_odr = lsc_out.beta
    Ri_3b = calc_r(xc_odr, yc_odr)
    residu_3b = sum((Ri_3b - R_odr) ** 2)
    residu2_3b = sum((Ri_3b ** 2 - R_odr ** 2) ** 2)

    return xc_odr, yc_odr, R_odr


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


def solve_quadratic(a, b, c):
    """
    Solve a quadratic equation for real roots.

    A `None` type is returnes if the equation has no real roots.
    """
    # calculate the discriminant
    d = (b ** 2) - (4 * a * c)

    if d < 0:
        print('No real solutions.')
        x1 = None
        x2 = None
    else:
        # find two solutions
        x1 = (-b - np.sqrt(d)) / (2 * a)
        x2 = (-b + np.sqrt(d)) / (2 * a)

    return [x1, x2]
