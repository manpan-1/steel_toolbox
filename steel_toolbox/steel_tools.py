# Engineering and structural steel related functions
"""
Module containing methods related to structural steel design.
"""
from math import pi, sqrt, sin, cos, atan, ceil
import numpy as np
import en_tools as en


# Calculate cross sectional properties. Two inputs required:
# A 2d (2, n) list with x,y values of nodes and a (3, m) 2d list for elements(first-node, second-node, thickness)

def cs_prop(nodes, elem):
    """
    Cross-sectional properties

    Calculate structural cross-sectional properties of profile. The
    profile is given in a form of 2 arrays, a list of nodes and element
    connectivity.

    Parameters
    ----------
    nodes : list
        List of nodes
    elem : list
        List of connecting elements

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    # Calculate cross sectional properties
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nele = len(elem[0])
    node = elem[0] + elem[1]
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
        sn = elem[0][i]
        fn = elem[1][i]
        # thickness of the element
        tt = tt + [elem[2][i]]
        # compute the coordinate of the mid point of the element
        xm = xm + [mean([nodes[0][sn], nodes[0][fn]])]
        ym = ym + [mean([nodes[1][sn], nodes[1][fn]])]
        # compute the dimension of the element
        xd = xd + [(nodes[0][fn] - nodes[0][sn])]
        yd = yd + [(nodes[1][fn] - nodes[1][sn])]
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
    moi_x = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(yd, side_length, tt)],
                                     [(a - yc) ** 2 * b * c for a, b, c in zip(ym, side_length, tt)])])
    moi_y = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(xd, side_length, tt)],
                                     [(a - xc) ** 2 * b * c for a, b, c in zip(xm, side_length, tt)])])
    moi_xy = sum([sum(a) for a in zip([a * b * c * d / 12 for a, b, c, d in zip(xd, yd, side_length, tt)],
                                      [(a - xc) * (b - yc) * c * d for a, b, c, d in zip(xm, ym, side_length, tt)])])

    if abs(moi_xy / area ** 2) < 1e-12:
        moi_xy = 0

    # Calculate angle of principal axes
    if moi_x == moi_y:
        theta_principal = pi / 2
    else:
        theta_principal = atan((-2 * moi_xy) / (moi_x - moi_y)) / 2

    # Change to centroid principal coordinates
    coord12 = [[a - xc for a in nodes[0]], [a - yc for a in nodes[1]]]
    coord12 = np.array([[cos(theta_principal), sin(theta_principal)],
                        [-sin(theta_principal), cos(theta_principal)]]).dot(nodes)

    # re-calculate cross sectional properties for the centroid
    for i in range(nele):
        sn = elem[0][i]
        fn = elem[1][i]
        # calculate the coordinate of the mid point of the element
        xm = xm + [mean([coord12[0][sn], coord12[0][fn]])]
        ym = ym + [mean([coord12[1][sn], coord12[1][fn]])]
        # calculate the dimension of the element
        xd = xd + [(coord12[0][fn] - coord12[0][sn])]
        yd = yd + [(coord12[1][fn] - coord12[1][sn])]

    # calculate the principal moment of inertia
    moi_1 = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(yd, side_length, tt)],
                                     [(a - yc) ** 2 * b * c for a, b, c in zip(ym, side_length, tt)])])
    moi_2 = sum([sum(a) for a in zip([a ** 2 * b * c / 12 for a, b, c in zip(xd, side_length, tt)],
                                     [(a - xc) ** 2 * b * c for a, b, c in zip(xm, side_length, tt)])])

    # Return values
    return area, xc, yc, moi_x, moi_y, moi_xy, moi_1, moi_2, theta_principal


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


def semi_closed_polygon(n_sides, radius, t, tg, rbend, nbend, l_lip):
    """
    Polygon sector nodes.

    Calculates the node coordinates for a a lipped polygon sector cross-section.


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

    Notes
    -----

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

    # First bend centre coordinates
    xcs[0] = x1 + l2 * np.cos(v2)
    ycs[0] = y1 + l2 * np.sin(v2)

    # Last bend
    v1 = phi_mids[-1] + pi / 2
    v2 = (v1 + phi[-1]) / 2
    l1 = (t + tg) / (2 * np.cos(v1 - phi[-1] - pi / 2))
    l2 = rs / np.sin(v2 - phi[-1])
    x1 = x[-1] + l1 * np.cos(v1)
    y1 = y[-1] + l1 * np.sin(v1)

    # Last bend centre coordinates
    xcs[1] = x1 + l2 * np.cos(v2)
    ycs[1] = y1 + l2 * np.sin(v2)

    # First and last bend arc points coordinates
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


def closed_polygon(
    n_sides,
    r_circle,
    p_classification,
    column_length=None,
    f_yield=None,
    fab_class=None
):
    """
    Geometry and properties of a polygonal cross-section.

    Calculate the cross-section geometrical characteristics (thickness, circumference, vertex coordinates) and
    and resistances according to EC1-1 - EC-1-5 and EC1-6 (plated cross section and cylindrical shell respectively).

    Parameters
    ----------
    n_sides : int
        The number of sides, >2
    r_circle : float
    p_classification : float
        Slenderness C/et of the polygon facets.
    column_length : float, optional
        Length of the member. The length is affecting only the shell resistance calculation acc. to EC1-6.
        Default is 2*pi*r.
    f_yield : float, optional
        Yield stress.
        Default is 355 MPa
    fab_class : {'fcA', 'fcB', 'fcC'}
        The fabrication quality class.
        Default is 'fcA', which implies u-max = 0.006*l_gx.

    Returns
    -------
    list : {float, float, float, list of floats, list of floats, float, float}
        r_circum, thickness, t_classification, x_corners, y_corners, N_pl_Rd, N_b_Rd_shell
        Where:
        r_circum is the radius of the circumscribed circle
        thickness of the profile
        t_classification is the tube slenderness, D/(et), as described in EC3-1-1
        x_corners, y_corners lists of x-y pairs of the polygon vertices
        N_pl_Rd is the resistance of the profile acc. to EC1-1 - EC3-1-5, plated profile
        N_b_Rd_shell is the shell buckling resistance acc. to EC3-1-6

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings.
       Brussels: CEN, 2005.
    .. [2] Eurocode 3: Design of steel structures - Part 1-5: Plated structural elements. Brussels: CEN, 2005.
    .. [3] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures.
        Brussels: CEN, 2006.
    """
    # Default vlaues
    if p_classification is None:
        p_classification = 42.
    else:
        p_classification = float(p_classification)

    if column_length is None:
        column_length = 2 * pi * r_circle
    else:
        column_length = float(column_length)

    if f_yield is None:
        f_yield = 355.
    else:
        f_yield = float(f_yield)

    if fab_class is None:
        fab_class = 'fcA'

    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)

    # Radius of the polygon's circumscribed circle
    r_circum = (pi * r_circle) / (n_sides * sin(pi / n_sides))

    # Diameter
    diam_circum = 2 * r_circum

    # Central angles
    theta = 2 * pi / n_sides

    # Width of each side
    w_side = diam_circum * sin(pi / n_sides)

    # Thickness the given class (classification as plated, not tube)
    thickness = (diam_circum * sin(theta / 2)) / (p_classification * epsilon)

    # Tube classification slenderness acc. to EC3-1-1
    t_classification = 2 * r_circle / (epsilon ** 2 * thickness)

    # Polar coordinate of ths polygon vertices on the cross-section plane
    phii = []
    for i_index in range(n_sides):
        phii.append(i_index * theta)

    # Polygon corners coordinates
    x_corners = r_circum * np.cos(phii)
    y_corners = r_circum * np.sin(phii)

    # Axial compression resistance , Npl
    n_pl_rd = n_sides * en.n_pl_rd(thickness, w_side, f_yield)

    # Compression resistance of equivalent cylindrical shell
    n_b_rd_shell = 2 * pi * r_circle * thickness * en.sigma_x_rd(
        thickness,
        r_circle,
        column_length,
        f_yield,
        fab_quality=fab_class,
        gamma_m1=1.
    )

    # Return values
    return r_circum, thickness, t_classification, tuple(x_corners), tuple(y_corners), n_pl_rd, n_b_rd_shell


def plastic_table(nominal=None):
    """
    Plasticity tables.

    Tables with plastic stress-strain curve values for different steels
    given a steel name.

    Parameters
    ----------
    nominal : string [optional]
        Steel name. Default value, 'S355'

    Notes
    -----

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
