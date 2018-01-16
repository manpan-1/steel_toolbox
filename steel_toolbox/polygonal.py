# -*- coding: utf-8 -*-

"""
A framework for the study of polygonal profiles.

"""
import numpy as np
import steel_toolbox.steel_design as sd
import steel_toolbox.lab_tests as lt
import steel_toolbox.analytic_geometry as ag
import steel_toolbox.scan_3D as s3d
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class PolygonalColumn(sd.Part):
    def __init__(self,
                 geometry=None,
                 cs_props=None,
                 material=None,
                 struct_props=None,
                 bc_loads=None,
                 specimen=None,
                 lab_test=None):

        self.specimen = specimen
        self.lab_test = lab_test

        super().__init__(
            geometry,
            cs_props,
            material,
            struct_props,
            bc_loads)

    @classmethod
    def from_geometry(
            cls,
            n_sides,
            r_circle,
            thickness,
            length,
            f_yield,
            fab_class
    ):
        geometry, cs_props, material, struct_props = PolygonalColumn.calc_properties(
            n_sides,
            r_circle,
            thickness,
            length,
            f_yield,
            fab_class
        )

        return cls(geometry, cs_props, material, struct_props)

    @classmethod
    def from_slenderness_and_thickness(
            cls,
            n_sides,
            p_classification,
            thickness,
            length,
            f_yield,
            fab_class
    ):
        """Calculate radius of a polygon cs for given sides, thickness, plate slenderness and yield strength"""

        # Epsilon for the material
        epsilon = np.sqrt(235. / f_yield)

        # Radius of the equal perimeter cylinder
        r_circle = n_sides * thickness * epsilon * p_classification / (2 * np.pi)

        geometry, cs_props, material, struct_props = PolygonalColumn.calc_properties(
            n_sides,
            r_circle,
            thickness,
            length,
            f_yield,
            fab_class
        )

        # Return radius of the cylinder of equal perimeter
        return cls(geometry, cs_props, material, struct_props)

    @classmethod
    def from_slenderness_and_radius(
            cls,
            n_sides,
            r_circle,
            p_classification,
            length,
            f_yield,
            fab_class
    ):
        """Calculate the thickness of a polygon cs for given sides, equivalent circle radius, plate slenderness and yield
        strength"""

        # Epsilon for the material
        epsilon = np.sqrt(235. / f_yield)

        # Calculate the thickness
        thickness = 2 * np.pi * r_circle / (n_sides * epsilon * p_classification)

        geometry, cs_props, material, struct_props = PolygonalColumn.calc_properties(
            n_sides,
            r_circle,
            thickness,
            length,
            f_yield,
            fab_class
        )

        # Return radius of the cylinder of equal perimeter
        return cls(geometry, cs_props, material, struct_props)

    @staticmethod
    def calc_properties(
            n_sides,
            r_circle,
            thickness,
            length,
            f_yield,
            fab_class):
        # Create material
        material = sd.Material(210000, 0.3, f_yield)
        epsilon = np.sqrt(235. / f_yield)

        # Radius of the polygon's circumscribed circle
        r_circum = (np.pi * r_circle) / (n_sides * np.sin(np.pi / n_sides))

        # Diameter
        diam_circum = 2 * r_circum

        # Central angles
        theta = 2 * np.pi / n_sides

        # Width of each side
        side_width = diam_circum * np.sin(np.pi / n_sides)

        # Polar coordinate of the polygon vertices on the cross-section plane
        phii = []
        for i_index in range(n_sides):
            phii.append(i_index * theta)

        # Polygon corners coordinates.
        x_corners = tuple(r_circum * np.cos(phii))
        y_corners = tuple(r_circum * np.sin(phii))

        # Cross-sectional properties
        nodes = [x_corners, y_corners]
        elem = [
            list(range(0, len(x_corners))),
            list(range(1, len(x_corners))) + [0],
            len(x_corners) * [thickness]
        ]

        cs_sketch = sd.CsSketch(nodes, elem)
        geometry = sd.Geometry(cs_sketch, length)
        cs_props = sd.CsProps.from_cs_sketch(cs_sketch)
        cs_props.max_dist = r_circum
        cs_props.min_dist = np.sqrt(r_circum ** 2 - (side_width / 2) ** 2)

        lmbda_y = sd.lmbda_flex(
            length,
            cs_props.area,
            cs_props.moi_1,
            kapa_bc=1.,
            e_modulus=material.e_modulus,
            f_yield=material.f_yield
        )

        lmbda_z = sd.lmbda_flex(
            length,
            cs_props.area,
            cs_props.moi_2,
            kapa_bc=1.,
            e_modulus=material.e_modulus,
            f_yield=material.f_yield
        )

        # Axial compression resistance , Npl
        n_pl_rd = n_sides * sd.n_pl_rd(thickness, side_width, f_yield)

        # Compression resistance of equivalent cylindrical shell
        n_b_rd_shell = 2 * np.pi * r_circle * thickness * sd.sigma_x_rd(
            thickness,
            r_circle,
            length,
            f_yield,
            fab_quality=fab_class,
            gamma_m1=1.
        )

        # Plate classification acc. to EC3-1-1
        p_classification = side_width / (epsilon * thickness)

        # Tube classification slenderness acc. to EC3-1-1
        t_classification = 2 * r_circle / (epsilon ** 2 * thickness)

        struct_props = sd.StructProps(
            t_classification=t_classification,
            p_classification=p_classification,
            lmbda_y=lmbda_y,
            lmbda_z=lmbda_z,
            n_pl_rd=n_pl_rd,
            n_b_rd_shell=n_b_rd_shell
        )

        return geometry, cs_props, material, struct_props

    def add_specimen(self, path):
        # Create a polygon column instance.
        sp1 = PolygSpecimen(thickness=3)

        # Add a center line for the specimen.
        sp1.add_centre_line([0, 0, 0], [0, 0, 1])

        # Add all sides and edges.
        # they consist of FlatFace and RoundedEdge instances.
        sp1.add_all_sides(16, path + 'side_', fit_planes=True, offset_to_midline=True)
        sp1.add_all_edges(16, path + 'edge_', ref_lines=True)

        # Find a series of points for each edge based on the scanned surface.
        sp1.find_real_edges(offset_to_midline=True)
        self.specimen = sp1

    def add_test(self, fh):
        self.lab_test = PolygTest.from_file(fh)


class PolygSpecimen:
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
        self.sides.append(s3d.FlatFace.from_pickle(filename))

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

        self.sides = [s3d.FlatFace.from_pickle(prefix + '{:02d}.pkl'.format(x)) for x in range(1, n_sides + 1)]

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
        self.edges.append(s3d.RoundedEdge.from_pickle(filename))

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
        self.edges = [s3d.RoundedEdge.from_pickle(prefix + '{:02d}.pkl'.format(x)) for x in range(1, n_sides + 1)]

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


class PolygTest(lt.Experiment):
    def __init__(self, header, data):
        super().__init__(header, data)

    def add_eccentricity(self, axis, column, moi, min_dist, thickness, young):
        """
        Calculate eccentricity.

        Adds a column in the data dictionary for the eccentricity of the load application on a given axis based on
        two opposite strain measurements.
        """

        self.data['e_' + axis] = []
        for load, strain1, strain2 in zip(self.data['Load'], self.data[column[0]], self.data[column[1]]):
            self.data['e_' + axis].append(self.eccentricity_from_strain(
                load * 1000,
                [strain1 * 1e-6, strain2 * 1e-6],
                moi,
                min_dist + thickness / 2,
                young)
            )

    @staticmethod
    def eccentricity_from_strain(load, strain, moi, dist, young=None):
        """
        Load eccentricity based on strain pairs.

        Calculate the eccentricity of an axial load to the neutral axis of a specimen for which pairs of strains are
        monitored with strain gauges. The eccentricity is calculated on one axis and requires the moment of inertia around
        it and a pair of strains on tow positions symmetric to the neutral axis.
        Elastic behaviour is assumed.
        """

        # Default values.
        if young is None:
            young = 210000.
        else:
            young = float(young)

        # Eccentricity.
        ecc = (strain[0] - strain[1]) * young * moi / (2 * load * dist)

        # Return
        return ecc


def semi_closed_polygon(n_sides, radius, t, tg, rbend, nbend, l_lip):
    """
    Polygon sector nodes.

    Calculates the node coordinates for a cross-section of the shape of
    a lipped polygon sector.

    Parameters
    ----------
    n_sides : int
        Number of sides of original polygon.
    radius : float
        Radius of the original polygon.
    t : float
        Thickness of the profile
    tg : float
        Thickness of the profile
    rbend : float
        Radius of the bended corners' arc
    nbend : int
        Number of nodes along the corners' arcs
    l_lip : int
        Length of the lips

    Returns
    -------
    list of lists
        Returns points for the entire profile (1st and 2nd returned values), and points for a single sector (3rd and 4th
        returned values).

    """

    # Angle corresponding to one face of the polygon
    theta = 2 * np.pi / n_sides

    # Angles of radii (measured from x-axis)
    phi = np.linspace(5 * np.pi / 6, np.pi / 6, int(n_sides / 3 + 1))

    # xy coords of the polygon's corners
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)

    # Bends

    # Distance between bending centre and corner
    lc = rbend / np.cos(theta / 2)

    # Centers of bending arcs
    xc = x[1:-1] - lc * np.cos(phi[1:-1])
    yc = y[1:-1] - lc * np.sin(phi[1:-1])

    # Angles of the edges' midlines (measured from x-axis)
    phi_mids = phi[0:-1] - theta / 2

    # xy coords of the arc's points
    xarc = [[0 for j in range(nbend + 1)] for i in range(int(n_sides / 3 - 1))]
    yarc = [[0 for j in range(nbend + 1)] for i in range(int(n_sides / 3 - 1))]
    for i in range(int(n_sides / 3 - 1)):
        for j in range(nbend + 1):
            xarc[i][j] = xc[i] + rbend * np.cos(phi_mids[i] - j * (theta / nbend))
            yarc[i][j] = yc[i] + rbend * np.sin(phi_mids[i] - j * (theta / nbend))

    # Start-end extensions
    # Bending radius
    rs = rbend / 2
    xcs = [0, 0]
    ycs = [0, 0]

    # First bend
    v1 = phi_mids[0] - np.pi / 2
    v2 = (phi[0] + phi_mids[0] - np.pi / 2) / 2
    l1 = (t + tg) / (2 * np.cos(phi[0] - phi_mids[0]))
    l2 = rs / np.sin(v2 - phi_mids[0] + np.pi / 2)
    x1 = x[0] + l1 * np.cos(v1)
    y1 = y[0] + l1 * np.sin(v1)

    # First bend centre coords
    xcs[0] = x1 + l2 * np.cos(v2)
    ycs[0] = y1 + l2 * np.sin(v2)

    # Last bend
    v1 = phi_mids[-1] + np.pi / 2
    v2 = (v1 + phi[-1]) / 2
    l1 = (t + tg) / (2 * np.cos(v1 - phi[-1] - np.pi / 2))
    l2 = rs / np.sin(v2 - phi[-1])
    x1 = x[-1] + l1 * np.cos(v1)
    y1 = y[-1] + l1 * np.sin(v1)

    # Last bend centre coords
    xcs[1] = x1 + l2 * np.cos(v2)
    ycs[1] = y1 + l2 * np.sin(v2)

    # First and last bend arc points coords
    xsarc = [[0 for j in range(nbend + 1)] for j in [0, 1]]
    ysarc = [[0 for j in range(nbend + 1)] for j in [0, 1]]
    for j in range(nbend + 1):
        xsarc[0][j] = xcs[0] + rs * np.cos(4 * np.pi / 3 + j * ((phi_mids[0] - np.pi / 3) / nbend))
        ysarc[0][j] = ycs[0] + rs * np.sin(4 * np.pi / 3 + j * ((phi_mids[0] - np.pi / 3) / nbend))
        xsarc[1][j] = xcs[1] + rs * np.cos(
            phi_mids[-1] + np.pi + j * ((phi[-1] + np.pi / 2 - phi_mids[-1]) / nbend))
        ysarc[1][j] = ycs[1] + rs * np.sin(
            phi_mids[-1] + np.pi + j * ((phi[-1] + np.pi / 2 - phi_mids[-1]) / nbend))

    # Points of the lips

    # Lip length according to bolt washer diameter

    # First lip
    xstart = [xsarc[0][0] + l_lip * np.cos(phi[0]), xsarc[0][0] + l_lip * np.cos(phi[0]) / 2]
    ystart = [ysarc[0][0] + l_lip * np.sin(phi[0]), ysarc[0][0] + l_lip * np.sin(phi[0]) / 2]

    # Last point
    xend = [xsarc[1][-1] + l_lip * np.cos(phi[-1]) / 2, xsarc[1][-1] + l_lip * np.cos(phi[-1])]
    yend = [ysarc[1][-1] + l_lip * np.sin(phi[-1]) / 2, ysarc[1][-1] + l_lip * np.sin(phi[-1])]

    # Collect the x, y values in a sorted 2xn array
    xarcs, yarcs = [], []
    for i in range(len(phi) - 2):
        xarcs = xarcs + xarc[i][:]
        yarcs = yarcs + yarc[i][:]

    x_sector = xstart + xsarc[0][:] + xarcs[:] + xsarc[1][:] + xend
    y_sector = ystart + ysarc[0][:] + yarcs[:] + ysarc[1][:] + yend

    # Copy-rotate the points of the first sector to create the entire CS
    # Rotation matrix
    rot_matrix = np.array([[np.cos(-2 * np.pi / 3), -np.sin(-2 * np.pi / 3)],
                           [np.sin(-2 * np.pi / 3), np.cos(-2 * np.pi / 3)]])

    # Dot multiply matrices
    coord1 = np.array([x_sector, y_sector])
    coord2 = rot_matrix.dot(coord1)
    coord3 = rot_matrix.dot(coord2)

    # Concatenate into a single xy array
    x_cs = np.concatenate([coord1[0], coord2[0], coord3[0]])
    y_cs = np.concatenate([coord1[1], coord2[1], coord3[1]])

    # Return matrices
    return [x_cs, y_cs, x_sector, y_sector]


def main():
    # First specimen
    number_of_sides = 16
    plate_classification = 30.
    thickness = 3.
    specimen_height = 700.
    yielt_stress = 700.
    fabrication_class = 'fcA'

    sp1 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    # Add data from the real specimen.
    # The data come from 3D scanning the fabricated specimen before performing the test.
    sp1.add_specimen('../../sp1/')

    # Add the data recorded during the test
    sp1.add_test('../data/experiments/sample_1.asc')

    # Similarly for the other 8 specimens

    number_of_sides = 16
    plate_classification = 40.
    thickness = 3.

    sp2 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp2.add_test('../data/experiments/sample_2_all_data_appended.asc')

    number_of_sides = 16
    plate_classification = 50.
    thickness = 3.

    sp3 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp3.add_test('../data/experiments/sample_3_realtest.asc')

    number_of_sides = 20
    plate_classification = 30.
    thickness = 3.

    sp4 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp4.add_test('../data/experiments/sample_4_realtest.asc')

    number_of_sides = 20
    plate_classification = 40.
    thickness = 3.

    sp5 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp5.add_test('../data/experiments/sample_5_realtest.asc')

    number_of_sides = 20
    plate_classification = 50.
    thickness = 2.

    sp6 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp6.add_test('../data/experiments/sample_6_realtest.asc')

    number_of_sides = 24
    plate_classification = 30.
    thickness = 3.

    sp7 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp7.add_test('../data/experiments/sample_7_realtest.asc')

    number_of_sides = 24
    plate_classification = 40.
    thickness = 2.

    sp8 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp8.add_test('../data/experiments/sample_8_realtest.asc')

    number_of_sides = 24
    plate_classification = 50.
    thickness = 2.

    sp9 = PolygonalColumn.from_slenderness_and_thickness(
        number_of_sides,
        plate_classification,
        thickness,
        specimen_height,
        yielt_stress,
        fabrication_class
    )

    sp9.add_test('../data/experiments/sample_9_realtest.asc')

    # Return all the specimens
    return [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9]
