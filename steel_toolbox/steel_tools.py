# Engineering and structural steel related functions
"""
Module containing methods related to structural steel design.
"""
from math import pi, sqrt, sin, cos, atan, ceil
import numpy as np
import steel_toolbox.en_tools as en


class Geometry:
    """
    Structural element geometry.

    Class for the geometric properties of a structural element.

    Parameters
    ----------
    cs_sketch : CsSketch object
        Cross-section sketch.
    length : float
        Member's length.
    """
    def __init__(self, cs_sketch, length):
        self.cs_sketch = cs_sketch
        self.length = length


class CsSketch:
    """
    Cross-section geometry.

    Parameters
    ----------
    nodes : list
        List of points.
    elem : list
        Element connectivity.
    """

    def __init__(self, nodes, elem):
        self.nodes = nodes
        self.elem = elem


class CsProps:
    """
    Cross-section properties

    Class for the mass properties of cross-sections. The properties can be calculated using the from_cs_sketch() method.

    Parameters
    ----------
    area : float
        Cross-sectional area.
    xc : float
        `x` coordinate of the gravity center.
    yc : float
        `y` coordinate of the gravity center.
    moi_xx : float
        Moment of inertia around `x` axis.
    moi_yy : float
        Moment of inertia around `y` axis.
    moi_xy : float
        Polar moment of inertia.
    theta_principal : float
        Rotation of the principal axes.
    moi_1 : float
        Moment of inertia around the major axis.
    moi_2 : float
        Moment of inertia around the minor axis.
    """

    def __init__(self,
                 area=None,
                 xc=None,
                 yc=None,
                 moi_xx=None,
                 moi_yy=None,
                 moi_xy=None,
                 theta_principal=None,
                 moi_1=None,
                 moi_2=None
                 ):

        self.area = area
        self.xc = xc
        self.yc = yc
        self.moi_xx = moi_xx
        self.moi_yy = moi_yy
        self.moi_xy = moi_xy
        self.theta_principal = theta_principal
        self.moi_1 = moi_1
        self.moi_2 = moi_2

    @classmethod
    def from_cs_sketch(cls, cs_sketch):
        """
        Cross-section calculator.

        Alternative constructor, calculates mass properties of a given sc sketch and returns a CsProps object.

        Parameters
        ----------
        cs_sketch : CsSketch object

        Notes
        -----

        """

        nele = len(cs_sketch.elem[0])
        node = cs_sketch.elem[0] + cs_sketch.elem[1]
        nnode = 0
        j = 0

        while node:
            i = [ii for ii, x in enumerate(node) if x == node[0]]
            for ii in sorted(i, reverse=True):
                del node[ii]
            if len(i) == 2:
                j += 1
            nnode += 1

        # classify the section type
        if j == nele:
            section = 'close'  # single cell
        elif j == nele - 1:
            section = 'open'  # singly-branched
        else:
            section = 'open'  # multi-branched

        # Calculate the cs-properties
        tt = []
        xm = []
        ym = []
        xd = []
        yd = []
        side_length = []
        for i in range(nele):
            sn = cs_sketch.elem[0][i]
            fn = cs_sketch.elem[1][i]
            # thickness of the element
            tt = tt + [cs_sketch.elem[2][i]]
            # compute the coordinate of the mid point of the element
            xm = xm + [mean([cs_sketch.nodes[0][sn], cs_sketch.nodes[0][fn]])]
            ym = ym + [mean([cs_sketch.nodes[1][sn], cs_sketch.nodes[1][fn]])]
            # compute the dimension of the element
            xd = xd + [(cs_sketch.nodes[0][fn] - cs_sketch.nodes[0][sn])]
            yd = yd + [(cs_sketch.nodes[1][fn] - cs_sketch.nodes[1][sn])]
            # compute the length of the element
            side_length = side_length + [sqrt(xd[i] ** 2 + yd[i] ** 2)]

        # calculate cross sectional area
        area = sum([a * b for a, b in zip(side_length, tt)])
        # compute the centroid
        xc = sum([a * b * c for a, b, c in zip(side_length, tt, xm)]) / area
        yc = sum([a * b * c for a, b, c in zip(side_length, tt, ym)]) / area

        if abs(xc / sqrt(area)) < 1e-12:
            xc = 0

        if abs(yc / sqrt(area)) < 1e-12:
            yc = 0

        # Calculate MOI
        moi_xx = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(yd, side_length, tt)],
                                          [(a - yc) ** 2 * b * c for a, b, c in
                                           zip(ym, side_length, tt)])])
        moi_yy = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(xd, side_length, tt)],
                                          [(a - xc) ** 2 * b * c for a, b, c in
                                           zip(xm, side_length, tt)])])
        moi_xy = sum(
            [sum(a) for a in zip([a * b * c * d / 12 for a, b, c, d in zip(xd, yd, side_length, tt)],
                                 [(a - xc) * (b - yc) * c * d for a, b, c, d in
                                  zip(xm, ym, side_length, tt)])])

        if abs(moi_xy / area ** 2) < 1e-12:
            moi_xy = 0

        # Calculate angle of principal axes
        if moi_xx == moi_yy:
            theta_principal = pi / 2
        else:
            theta_principal = atan(
                (-2 * moi_xy) / (moi_xx - moi_yy)) / 2

        # Change to centroid principal coordinates
        coord12 = [[a - xc for a in cs_sketch.nodes[0]],
                   [a - yc for a in cs_sketch.nodes[1]]]
        coord12 = np.array([[cos(theta_principal), sin(theta_principal)],
                            [-sin(theta_principal), cos(theta_principal)]]).dot(
            cs_sketch.nodes)

        # re-calculate cross sectional properties for the centroid
        for i in range(nele):
            sn = cs_sketch.elem[0][i]
            fn = cs_sketch.elem[1][i]
            # calculate the coordinate of the mid point of the element
            xm = xm + [mean([coord12[0][sn], coord12[0][fn]])]
            ym = ym + [mean([coord12[1][sn], coord12[1][fn]])]
            # calculate the dimension of the element
            xd = xd + [(coord12[0][fn] - coord12[0][sn])]
            yd = yd + [(coord12[1][fn] - coord12[1][sn])]

        # calculate the principal moment of inertia
        moi_1 = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(yd, side_length, tt)],
                                         [(a - yc) ** 2 * b * c for a, b, c in
                                          zip(ym, side_length, tt)])])
        moi_2 = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(xd, side_length, tt)],
                                         [(a - xc) ** 2 * b * c for a, b, c in
                                          zip(xm, side_length, tt)])])

        return cls(
            area=area,
            xc=xc,
            yc=yc,
            moi_xx=moi_xx,
            moi_yy=moi_yy,
            moi_xy=moi_xy,
            theta_principal=theta_principal,
            moi_1=moi_1,
            moi_2=moi_2
        )


class Material:
    """
    Material properties.

    Parameters
    ----------
    e_modulus : float
        Modulus of elasticity.
    poisson : float
        Poisson's ratio.
    f_yield : float
        Yield stress
    plasticity : tuple
        Plasticity table (tuple of stress-plastic strain pairs).
        By default, no plasticity is considered.

    """
    def __init__(self, e_modulus, poisson, f_yield, plasticity=None):
        self.e_modulus = e_modulus
        self.poisson = poisson
        self.f_yield = f_yield
        self.plasticity = plasticity

#TODO change default value to S235
    def plastic_table(nominal=None):
        """
        Plasticity tables.

        Tables with plastic stress-strain curve values for different steels
        given a steel name, e.g 'S355'

        Parameters
        ----------
        nominal : string [optional]
            Steel name. Default value, 'S355'

        Attributes
        ----------

        Notes
        -----

        References
        ----------

        """
        if nominal is None:
            nominal = 'S355'

        if nominal is 'S355':
            table = (
                (381.1, 0.0),
                (391.2, 0.0053),
                (404.8, 0.0197),
                (418.0, 0.0228),
                (444.2, 0.0310),
                (499.8, 0.0503),
                (539.1, 0.0764),
                (562.1, 0.1009),
                (584.6, 0.1221),
                (594.4, 0.1394),
                (5961, 1.)
            )

        if nominal is 'S650':
            table = (
                (760., 0.0),
                (770., 0.022),
                (850., 0.075),
                (900., 0.1),
                (901., 1.)
            )

        return table

    @classmethod
    def from_nominal(cls, nominal_strength=None):

        """
        Alternative constructor creating a steel material given of a given nominal strength.
        """
        if nominal_strength is None:
            f_yield = 235.
        else:
            f_yield = float(nominal_strength.replace('S', ''))
            plasticity = cls.plastic_table(nominal = nominal_strength)
        return cls(210000., 0.3, f_yield, plasticity=plasticity)


class BCs:
    def __init__(self, bcs):
        self.bcs = bcs

    @classmethod
    def from_hinged(cls):
        return cls([[1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0]])


class StructProps:
    """
    Structural properties of a member.

    Parameters
    ----------
    t_classification : float, optional
        Classification of a tube, d/(t^2*e)
    p_classification : float, optional
        Classification of a plate, c/(t*e)
    lmbda_y : float, optional
        Flexural slenderness on the strong axis.
    lmbda_z : float, optional
        Flexural slenderness on the weak axis.
    n_pl_rd : float, optional
        Plastic axial compression resistance.
    n_b_rd_shell : float, optional
        Shell buckling resistance
    """
    def __init__(self,
                 t_classification=None,
                 p_classification=None,
                 lmbda_y=None,
                 lmbda_z=None,
                 n_pl_rd=None,
                 n_b_rd_shell=None
                 ):

        self.t_classification = t_classification
        self.p_classification = p_classification
        self.lmbda_y = lmbda_y
        self.lmbda_z = lmbda_z
        self.n_pl_rd = n_pl_rd
        self.n_b_rd_shell = n_b_rd_shell


class Part:
    """
    Structural part.

    Class describing a structural part, including geometry, boundary conditions loads and resistance.

    Parameters
    ----------
    geometry : Geometry object, optional
    cs_props : CsProps object, optional
    material : Material object, optional
    struct_props : StructProps object, optional
    bc_loads: BCs object, optional
    """

    def __init__(self,
                 geometry=None,
                 cs_props=None,
                 material=None,
                 struct_props=None,
                 bc_loads=None
                 ):
        self.geometry = geometry
        self.cs_props = cs_props
        self.material = material
        self.bc_loads = bc_loads
        self.struct_props = struct_props


class PolygonalColumn(Part):
    def __init__(self,
                 geometry=None,
                 cs_props=None,
                 material=None,
                 struct_props=None,
                 bc_loads=None):

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
        epsilon = sqrt(235. / f_yield)

        # Radius of the equal perimeter cylinder
        r_circle = n_sides * thickness * epsilon * p_classification / (2 * pi)

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
        epsilon = sqrt(235. / f_yield)

        # Calculate the thickness
        thickness = 2 * pi * r_circle / (n_sides * epsilon * p_classification)

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
        material = Material(210000, 0.3, f_yield)
        epsilon = sqrt(235. / f_yield)

        # Radius of the polygon's circumscribed circle
        r_circum = (pi * r_circle) / (n_sides * sin(pi / n_sides))

        # Diameter
        diam_circum = 2 * r_circum

        # Central angles
        theta = 2 * pi / n_sides

        # Width of each side
        side_width = diam_circum * sin(pi / n_sides)

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

        cs_sketch = CsSketch(nodes, elem)
        geometry = Geometry(cs_sketch, length)
        cs_props = CsProps.from_cs_sketch(cs_sketch)
        cs_props.max_dist = r_circum
        cs_props.min_dist = sqrt(r_circum ** 2 - (side_width / 2) ** 2)

        lmbda_y = en.lmbda_flex(
            length,
            cs_props.area,
            cs_props.moi_1,
            kapa_bc=1.,
            e_modulus=material.e_modulus,
            f_yield=material.f_yield
        )

        lmbda_z = en.lmbda_flex(
            length,
            cs_props.area,
            cs_props.moi_2,
            kapa_bc=1.,
            e_modulus=material.e_modulus,
            f_yield=material.f_yield
        )

        # Axial compression resistance , Npl
        n_pl_rd = n_sides * en.n_pl_rd(thickness, side_width, f_yield)

        # Compression resistance of equivalent cylindrical shell
        n_b_rd_shell = 2 * pi * r_circle * thickness * en.sigma_x_rd(
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

        struct_props = StructProps(
            t_classification=t_classification,
            p_classification=p_classification,
            lmbda_y=lmbda_y,
            lmbda_z=lmbda_z,
            n_pl_rd=n_pl_rd,
            n_b_rd_shell=n_b_rd_shell
        )

        return geometry, cs_props, material, struct_props


# Calculate diameter of washer for a given bolt
def bolt2washer(m_bolt):
    """
    Washer diameter.

    Return the diameter of the washer for a given bolt diameter.
    The calculation is based on a function derived from linear regression
    on ENXXXXXXX[REF].

    Parameters
    ----------
    m_bolt : float
        Bolt diameter

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    d_washer = ceil(1.5893 * m_bolt + 5.1071)
    return d_washer


# mean value
def mean(numbers):
    """
    Mean value.

    Calculate the average for a list of numbers.

    Parameters
    ----------
    numbers : list

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    return float(sum(numbers)) / max(len(numbers), 1)


# Calculate xy of nodes for a given polygonal profile
# Returns points for the entire profile (1st and 2nd returned values)
# and points for a single sector (3rd and 4th returned values)
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

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    # Angle corresponding to one face of the polygon
    theta = 2 * pi / n_sides

    # Angles of radii (measured from x-axis)
    phi = np.linspace(5 * pi / 6, pi / 6, n_sides / 3 + 1)

    # xy coords of the polygon's corners
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)

    # Bends

    # Distance between bending centre and corner
    lc = rbend / np.cos(theta / 2)

    # Centers of bending arcs
    xc = x[1:-1] - lc * np.cos(phi[1:-1])
    yc = y[1:-1] - lc * np.sin(phi[1:-1])

    # Bending arc angle
    theta_b = pi - theta

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
    v1 = phi_mids[0] - pi / 2
    v2 = (phi[0] + phi_mids[0] - pi / 2) / 2
    l1 = (t + tg) / (2 * np.cos(phi[0] - phi_mids[0]))
    l2 = rs / np.sin(v2 - phi_mids[0] + pi / 2)
    x1 = x[0] + l1 * np.cos(v1)
    y1 = y[0] + l1 * np.sin(v1)

    # First bend centre coords
    xcs[0] = x1 + l2 * np.cos(v2)
    ycs[0] = y1 + l2 * np.sin(v2)

    # Last bend
    v1 = phi_mids[-1] + pi / 2
    v2 = (v1 + phi[-1]) / 2
    l1 = (t + tg) / (2 * np.cos(v1 - phi[-1] - pi / 2))
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
        xsarc[0][j] = xcs[0] + rs * np.cos(4 * pi / 3 + j * ((phi_mids[0] - pi / 3) / nbend))
        ysarc[0][j] = ycs[0] + rs * np.sin(4 * pi / 3 + j * ((phi_mids[0] - pi / 3) / nbend))
        xsarc[1][j] = xcs[1] + rs * np.cos(
            phi_mids[-1] + pi + j * ((phi[-1] + pi / 2 - phi_mids[-1]) / nbend))
        ysarc[1][j] = ycs[1] + rs * np.sin(
            phi_mids[-1] + pi + j * ((phi[-1] + pi / 2 - phi_mids[-1]) / nbend))

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
    rot_matrix = np.array([[cos(-2 * pi / 3), -sin(-2 * pi / 3)],
                           [sin(-2 * pi / 3), cos(-2 * pi / 3)]])

    # Dot multiply matrices
    coord1 = np.array([x_sector, y_sector])
    coord2 = rot_matrix.dot(coord1)
    coord3 = rot_matrix.dot(coord2)

    # Concatenate into a single xy array
    x_cs = np.concatenate([coord1[0], coord2[0], coord3[0]])
    y_cs = np.concatenate([coord1[1], coord2[1], coord3[1]])

    # Return matrices
    return x_cs, y_cs, x_sector, y_sector


def eccentricity(load, strain, moi, dist, young=None):
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
