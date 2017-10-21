# -*- coding: utf-8 -*-
"""

A collection of structural steel related functions based on EN1993.

"""

from math import sqrt, pi, sin, cos, atan


# Define classes
# Cross-section class
class CrossSection:
    """
    Cross-section properties.

    Geometric and inertia parameters for structural purpose cross-section.

    Parameters
    ----------
    area : float
        [mm^2] Cross-sectional area.
    moi_y : float
        [mm^4] Moment of inertia around y-axis.
        y-axis on the centre of gravity but not necessarily principal.
    moi_z : float

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    def __init__(self, area, moi_y, moi_z):
        self.area = area
        self.moi_y = moi_y
        self.moi_z = moi_z


class Member:
    def __init__(self, profile, length, f_yield):
        self.profile = profile
        self.length = length
        self.f_yield = f_yield


# SIMPLY SUPPORTED PLATE

def n_pl_rd(
    thickness,
    width,
    f_yield,
    psi=None
):
    # Docstring
    """
    Plastic design resistance of a plate.

    Calculates the resistance of a plate according to EN1993-1-1 and
    EN1993-1-5. The plate is assumed simply supported.

    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    f_yield : float
        [MPa] Yield stress
    psi : float, optional
        [_] Ratio of the min over max stress for a linear distribution,
        (sigma_min / sigma_max)
        Default = 1, which implies a uniform distribution

    Returns
    -------
    float
        [N] Plastic design resistance

    Notes
    -----
    To be extended to include cantilever plate (outstand members)

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings.
        Brussels: CEN, 2005.
    .. [2] Eurocode 3: Design of steel structures - Part 1-5: Plated structural elements. Brussels: CEN, 2005.

    """

    # Convert inputs to floats
    thickness, width, f_yield = float(thickness), float(width), float(f_yield)

    # Default value for psi
    if psi is None:
        psi = 1.
    else:
        psi = float(psi)

    # Calculate kapa_sigma
    k_sigma = 8.2 / (1.05 + psi)

    # Aeff calculation.
    # Reduction factor for the effective area of the profile acc. to EC3-1-5
    classification = width / (thickness * sqrt(235 / f_yield))
    lambda_p = classification / (28.4 * sqrt(k_sigma))
    if lambda_p > 0.673 and plate_class(thickness, width, f_yield) == 4:
        rho = (lambda_p - 0.055 * (3 + psi)) / lambda_p ** 2
    else:
        rho = 1.

    # Effective area
    a_eff = rho * thickness * width

    # Axial compression resistance , Npl
    nn_pl_rd = a_eff * f_yield

    # Return value
    return nn_pl_rd


def plate_class(
    thickness,
    width,
    f_yield
):
    # Docstring
    """
    Plate classification.

    Returnes the class for a given plate, according to EN1993-1-1.
    Currently works for simply supported plates under pure compression.

    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    f_yield : float
        [MPa] Yield stress

    Returns
    -------
    int
        [_] Class number

    Notes
    -----
    To be extended to include the rest of the cases of Table 5.3 [1].
    Members under combined axial and bending and outstand members.

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005

    """

    # Convert inputs to floats
    width, thickness, f_yield = float(width), float(thickness), float(f_yield)

    # Calculate classification
    classification = width / (thickness * sqrt(235. / f_yield))
    if classification <= 33.:
        p_class = 1
    elif classification <= 38.:
        p_class = 2
    elif classification <= 42.:
        p_class = 3
    else:
        p_class = 4

    # Return value
    return p_class


def sigma_cr_plate(
    thickness,
    width,
    psi=None
):
    # Docstring
    """
    Critical stress of a plate.

    Calculates the critical stress for a simply supported plate.

    Parameters
    ----------
    thickness : float
        [mm] Plate thickness
    width : float
        [mm] Plate width
    psi : float, optional
        [_] Ratio of the min over max stress for a linear distribution,
        (sigma_min / sigma_max)
        Default = 1, which implies a uniform distribution

    Returns
    -------
    float
        [MPa] Plate critical stress

    Notes
    -----
    To be extended to include cantilever plate (outstand members)

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-5: Plated structural elements. Brussels: CEN, 2005.

    """
    # Convert inputs to floats
    thickness, width = float(thickness), float(width)

    # Default value for psi
    if psi is None:
        psi = 1.
    else:
        psi = float(psi)

    # Calculate kapa_sigma
    k_sigma = 8.2 / (1.05 + psi)

    # Elastic critical stress acc. to EN3-1-5 Annex A
    sigma_e = 190000 * (thickness / width) ** 2
    sigma_cr = sigma_e * k_sigma

    # Return value
    return sigma_cr


# CYLINDRICAL SHELLS

def sigma_x_rd(
    thickness,
    radius,
    length,
    f_y_k,
    fab_quality=None,
    gamma_m1=None
):
    # Docstring
    """
    Meridional design buckling stress.

    Calculates the meridional buckling stress for a cylindrical shell
    according to EN1993-1-6 [1].

    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length
    f_y_k : float
        [MPa] Characteristic yield strength
    fab_quality : str, optional
        [_] Fabrication quality class. Accepts: 'fcA', 'fcB', 'fcC'
        The three classes correspond to .006, .010 and .016 times the
        width of a dimple on the shell.
        Default = 'fcA', which implies excelent fabrication
    gamma_m1 : int, optional
        [_] Partial safety factor
        Default = 1.1

    Returns
    -------
    float
        [MPa] Meridional buckling stress

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures.
        Brussels: CEN, 2006.

    """

    # Default values
    if fab_quality is None:
        fab_quality = 'fcA'

    if gamma_m1 is None:
        gamma_m1 = 1.1
    else:
        gamma_m1 = float(gamma_m1)

    # Fabrication quality class acc. to table D2
    if fab_quality is 'fcA':
        q_factor = 40.
    elif fab_quality is 'fcB':
        q_factor = 25.
    elif fab_quality is 'fcC':
        q_factor = 16.
    else:
        print('Invalid fabrication class input. Choose between \'fcA\', \'fcB\' and \'fcC\' ')
        return

    # Critical meridinal stress, calculated on separate function
    sigma_cr = sigma_x_rcr(thickness, radius, length)

    # Shell slenderness
    lmda = sqrt(f_y_k / sigma_cr[0])
    delta_w_k = (1. / q_factor) * sqrt(radius / thickness) * thickness
    alpha = 0.62 / (1 + 1.91 * (delta_w_k / thickness) ** 1.44)
    beta = 0.6
    eta = 1.
    if sigma_cr[1] is 'long':
        # For long cylinders, a formula is suggested for lambda, EC3-1-6 D1.2.2(4)
        # Currently, the general form is used. to be fixed.
        lmda_0 = 0.2
        # lmda_0 = 0.2 + 0.1 * (sigma_e_M / sigma_e)
    else:
        lmda_0 = 0.2

    lmda_p = sqrt(alpha / (1. - beta))

    # Buckling reduction factor, chi
    if lmda <= lmda_0:
        chi = 1.
    elif lmda < lmda_p:
        chi = 1. - beta * ((lmda - lmda_0) / (lmda_p - lmda_0)) ** eta
    else:
        chi = alpha / (lmda ** 2)

    # Buckling stress
    sigma_rk = chi * f_y_k
    sigma_rd = sigma_rk / gamma_m1

    # Return value
    return sigma_rd


def n_cr_shell(
    thickness,
    radius,
    length
):
    # Docstring
    """
    Critical compressive load for cylindrical shell.

    Calculates the critical load for a cylindrical shell under pure
    compression and assumes uniform stress distribution. Calculation
    according to EN1993-1-6 [1], Annex D.

    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length

    Returns
    -------
    float
        [N] Critical load

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures.
        Brussels: CEN, 2006.

    """

    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)

    # Elastic critical load acc to EN3-1-6 Annex D
    nn_cr_shell = 2 * pi * radius * thickness * sigma_x_rcr(thickness, radius, length)[0]

    # Return value
    return nn_cr_shell


def sigma_x_rcr(
    thickness,
    radius,
    length
):
    # Docstring
    """
    Critical meridional stress for cylindrical shell.

    Calculates the critical load for a cylindrical shell under pure
    compression and assumes uniform stress distribution. Calculation
    according to EN1993-1-6 [1], Annex D.

    Parameters
    ----------
    thickness : float
        [mm] Shell thickness
    radius : float
        [mm] Cylinder radius
    length : float
        [mm] Cylnder length

    Returns
    -------
    float
        [N] Critical load

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures.
        Brussels: CEN, 2006.

    """
    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)

    # Elastic critical load acc. to EN3-1-6 Annex D
    omega = length / sqrt(radius * thickness)
    if 1.7 <= omega <= 0.5 * (radius / thickness):
        c_x = 1.
        length_category = 'medium'
    elif omega < 1.7:
        c_x = 1.36 - (1.83 / omega) + (2.07 / omega ** 2)
        length_category = 'short'
    else:
        # c_x_b is read on table D.1 of EN3-1-5 Annex D acc. to BCs
        # BC1 - BC1 is used on the Abaqus models (both ends clamped, see EN3-1-5 table 5.1)
        c_x_b = 6.
        c_x_n = max((1 + 0.2 * (1 - 2 * omega * thickness / radius) / c_x_b), 0.6)
        c_x = c_x_n
        length_category = 'long'

    # Calculate critical stress, eq. D.2 on EN3-1-5 D.1.2.1-5
    sigma_cr = 0.605 * 210000 * c_x * thickness / radius

    # Return value
    return sigma_cr, length_category


def fabclass_2_umax(fab_class=None):
    # Docstring
    """
    Max dimple displacement.

    Returns the maximum displacement for a dimple imperfection on a cylindrical shell. The values are taken from table
    8.4 of EN1993-1-6[1] for a given fabrication quality class, A, B or C.

    Parameters
    ----------
    fab_class : {'fcA', 'fcB', 'fcC'}
        The fabrication quality class.

    Returns
    -------
    float
        u_max / l, where u_max is the maximum deviation and l the dimple's size (circumferencial or meridional)

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures.
        Brussels: CEN, 2006.

    """
    # default values
    if fab_class is None:
        fab_class = 'fcA'

    # Assign imperfection amplitude, u_max acc. to the fabrication class
    if fab_class is 'fcA':
        u_max = 0.006
    elif fab_class is 'fcB':
        u_max = 0.010
    else:
        u_max = 0.016

    # Return values
    return u_max


# OVERALL BUCKLING
def n_cr_flex(
    length,
    moi_y,
    kapa_bc=None,
    e_modulus=None
):
    # Docstring
    """
    Euler's critical load.

    Calculates the critical load for flexural buckling of a given column.
    A single direction is considered. If more directions are required
    (e.g the two principal axes), the function has to be called multiple
    times. For torsional mode critical load use n_cr_tor(), and for
    flexural-torsional critical load use n_cr_flex_tor()

    Parameters
    ----------
    length : float
        [mm] Column length.
    moi_y : float
        [mm^4] Moment of inertia.
    kapa_bc : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column.
    e_modulus : float, optional
        [MPa] Modulus of elasticity.
        Default = 210000., typical value for steel.

    Returns
    -------
    float
        [N] Critical load.

    """
    # default values
    if kapa_bc is None:
        kapa_bc = 1.
    else:
        kapa_bc = float(kapa_bc)

    if e_modulus is None:
        e_modulus = 210000.
    else:
        e_modulus = float(e_modulus)

    # Euler's critical load
    nn_cr_flex = (pi ** 2) * e_modulus * moi_y / (kapa_bc * length) ** 2

    # Return the result
    return nn_cr_flex


def n_cr_tor(
    length,
    area,
    moi_y0,
    moi_z0,
    moi_torsion,
    moi_warp,
    y_0=None,
    z_0=None,
    e_modulus=None,
    poisson=None,
):
    # Docstring
    """
    Torsional elastic critical load

    Calculates the torsional elastic critical load for a hinged column.
    The input values are refering to the principal axes. For flexural
    buckling (Euler cases) use n_cr_flex. For the combined
    flexural-torsional modes use n_cr_flex_tor.

    Parameters
    ----------
    length : float
        [mm] Column length.
    area : float
        [mm^2] Cross-sectional area.
    moi_y0 : float
        [mm^4] Moment of inertia around `y`-axis.
        `y`-axis on the centre of gravity but not necessarily principal.
    moi_z0 : float
        [mm^4] Moment of inertia around `z`-axis.
        `z`-axis on the centre of gravity but not necessarily principal.
    moi_torsion : float
        [mm^4] Saint Venant constant.
    moi_warp : float
        [mm^6] Torsion constant.
    y_0 : float, optional
        [mm] Distance on `y`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    z_0 : float, optional
        [mm] Distance on `z`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    e_modulus : float, optional
        [MPa] Modulus of elasticity.
        Default = 210000., general steel.
    poisson : float, optional
        [_] Young's modulus of elasticity.
        Default = 0.3, general steel.

    Returns
    -------
    float
        [N] Flexural-torsional critical load.

    Notes
    -----
    The torsional critical load is calculated as:

    .. math:: N_{cr, tor} = {GJ + {\pi^2EI_w\over{L^2}}\over{r^2}}

    Where:
        :math:`E`    : Elasticity modulus

        :math:`G`    : Shear modulus

        :math:`J`    : Torsional constant (Saint Venant)

        :math:`I_w`  : Warping constant

        :math:`r^2=(moi_y + moi_z)/A + x_0^2 + y_0^2`

        :math:`x_0, y_0`  : Shear centre coordinates on the principal coordinate system


    References
    ----------
    ..[1]N. S. Trahair, Flexural-torsional buckling of structures, vol. 6. CRC Press, 1993.
    ..[2]NS. Trahair, MA. Bradford, DA. Nethercot, and L. Gardner, The behaviour and design of steel structures to EC3, 4th edition. London; New York: Taylor & Francis, 2008.

    """
    # default values
    if y_0 is None:
        y_0 = 0
    else:
        y_0 = float(y_0)

    if z_0 is None:
        z_0 = 0
    else:
        z_0 = float(z_0)

    if e_modulus is None:
        e_modulus = 210000.
    else:
        e_modulus = float(e_modulus)

    if poisson is None:
        poisson = 0.3
    else:
        poisson = float(poisson)

    # Shear modulus
    g_modulus = e_modulus / (2 * (1 + poisson))

    # Polar radius of gyration.
    i_pol = sqrt((moi_y0 + moi_z0) / area)
    moi_zero = sqrt(i_pol ** 2 + y_0 ** 2 + z_0 ** 2)

    # Calculation of critical torsional load.
    nn_cr_tor = (1 / moi_zero ** 2) * (g_modulus * moi_torsion + (pi ** 2 * e_modulus * moi_warp / length ** 2))

    # Return the result
    return nn_cr_tor


def n_cr_flex_tor(
    length,
    area,
    moi_y,
    moi_z,
    moi_yz,
    moi_torsion,
    moi_warp,
    y_sc=None,
    z_sc=None,
    e_modulus=None,
    poisson=None,
):
    # Docstring
    """
    Flexural-Torsional elastic critical load

    Calculates the critical load for flexural-torsional buckling of a
    column with hinged ends. The returned value is the minimum of the
    the three flexural-torsional and the indepedent torsional mode, as
    dictated in EN1993-1-1 6.3.1.4 [1]. (for further details, see Notes).

    Parameters
    ----------
    length : float
        [mm] Column length.
    area : float
        [mm^2] Cross-sectional area.
    moi_y : float
        [mm^4] Moment of inertia around `y`-axis.
        `y`-axis on the centre of gravity but not necessarily principal.
    moi_z : float
        [mm^4] Moment of inertia around `z`-axis.
        `z`-axis on the centre of gravity but not necessarily principal.
    moi_yz : float
        [mm^4] Product of inertia.
    moi_torsion : float
        [mm^4] Saint Venant constant.
    moi_warp : float
        [mm^6] Torsion constant.
    y_sc : float, optional
        [mm] Distance on `y`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    z_sc : float, optional
        [mm] Distance on `z`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    e_modulus : float, optional
        [MPa] Modulus of elasticity.
        Default = 210000., general steel.
    poisson : float, optional
        [_] Young's modulus of elasticity.
        Default = 0.3, general steel.

    Returns
    -------
    float
        [N] Flexural-torsional critical load.

    Notes
    -----
    The flexural-torsional critical loads are calculated as a combination
    of the three independent overall buckling modes:
    i)   flexural around the major axis,
    ii)  flexural around the minor axis,
    iii) Torsional buckling (around x-axis).

    First, the cs-properties are described on the principal axes. Then
    the three independent  modes are calculated. The combined
    flexural-torsional modes are calculated as the roots of a 3rd order
    equation, as given in [1], [2]. The minimum of the torsional and the
    three combined modes is returned (the two independent flexural modes
    are not considered; for critical load of pure flexural mode use
    'n_cr_flex').

    References
    ----------
    ..[1]N. S. Trahair, Flexural-torsional buckling of structures, vol. 6. CRC Press, 1993.
    ..[2]NS. Trahair, MA. Bradford, DA. Nethercot, and L. Gardner, The behaviour and design of steel structures to EC3, 4th edition. London; New York: Taylor & Francis, 2008.

    """
    # default values
    if y_sc is None:
        y_sc = 0
    else:
        y_sc = float(y_sc)

    if z_sc is None:
        z_sc = 0
    else:
        z_sc = float(z_sc)

    if e_modulus is None:
        e_modulus = 210000.
    else:
        e_modulus = float(e_modulus)

    if poisson is None:
        poisson = 0.3
    else:
        poisson = float(poisson)

    # Angle of principal axes
    if abs(moi_y - moi_z) < 1e-20:
        theta = pi / 4
    else:
        theta = -atan((2 * moi_yz) / (moi_y - moi_z)) / 2

    # Distance of the rotation centre to the gravity centre on the
    # principal axes coordinate system
    y_0 = y_sc * cos(-theta) - z_sc * sin(-theta)
    z_0 = z_sc * cos(-theta) + y_sc * sin(-theta)

    # Moment of inertia around principal axes.
    moi_y0 = (moi_y + moi_z) / 2 + sqrt(((moi_y - moi_z) / 2) ** 2 + moi_yz ** 2)
    moi_z0 = (moi_y + moi_z) / 2 - sqrt(((moi_y - moi_z) / 2) ** 2 + moi_yz ** 2)

    # Polar radius of gyration.
    i_pol = sqrt((moi_y0 + moi_z0) / area)
    moi_zero = sqrt(i_pol ** 2 + y_0 ** 2 + z_0 ** 2)

    # Independent critical loads for flexural and torsional modes.
    n_cr_max = (pi ** 2 * e_modulus * moi_y0) / (length ** 2)
    n_cr_min = (pi ** 2 * e_modulus * moi_z0) / (length ** 2)
    n_tor = n_cr_tor(
        length,
        area,
        moi_y0,
        moi_z0,
        moi_torsion,
        moi_warp=moi_warp,
        y_0=y_0,
        z_0=z_0,
        e_modulus=e_modulus,
        poisson=poisson
    )

    # Coefficients of the 3rd order equation for the critical loads
    # The equation is in the form aaaa * N ^ 3 - bbbb * N ^ 2 + cccc * N - dddd
    aaaa = moi_zero ** 2 - y_0 ** 2 - z_0 ** 2
    bbbb = ((n_cr_max + n_cr_min + n_tor) * moi_zero ** 2) - (n_cr_min * y_0 ** 2) - (n_cr_max * z_0 ** 2)
    cccc = moi_zero ** 2 * (n_cr_min * n_cr_max) + (n_cr_min * n_tor) + (n_tor * n_cr_max)
    dddd = moi_zero ** 2 * n_cr_min * n_cr_max * n_tor

    det_3 = (
        4 * (-bbbb ** 2 + 3 * aaaa * cccc) ** 3 + (2 * bbbb ** 3 - 9 * aaaa * bbbb * cccc + 27 * aaaa ** 2 * dddd) ** 2
    )

    if det_3 < 0:
        det_3 = -1. * det_3
        cf = 1j
    else:
        cf = 1

    # Critical load
    # The following n_cr formulas are the roots of the 3rd order equation of the global critical load

    n_cr_1 = bbbb/(3.*aaaa) - (2**(1./3)*(-bbbb**2 + 3*aaaa*cccc))/                     \
        (3.*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                      \
        (cf*sqrt(det_3)))**(1./3)) + (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +  \
        (cf*sqrt(det_3)))**(1./3)/(3.*2**(1./3)*aaaa)

    n_cr_2 = bbbb/(3.*aaaa) + ((1 + (0 + 1j)*sqrt(3))*(-bbbb**2 + 3*aaaa*cccc))/        \
        (3.*2**(2./3)*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +            \
        (cf*sqrt(det_3)))**(1./3)) - ((1 - (0 + 1j)*sqrt(3))*                           \
        (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                               \
        (cf*sqrt(det_3)))**(1./3))/(6.*2**(1./3)*aaaa)

    n_cr_3 = bbbb/(3.*aaaa) + ((1 - (0 + 1j)*sqrt(3))*(-bbbb**2 + 3*aaaa*cccc))/        \
        (3.*2**(2./3)*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +            \
        (cf*sqrt(det_3)))**(1./3)) - ((1 + (0 + 1j)*sqrt(3))*                           \
        (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                               \
        (cf*sqrt(det_3)))**(1./3))/(6.*2**(1./3)*aaaa)


    # Lowest root is the critical load
    nn_cr_flex_tor = min(abs(n_cr_1), abs(n_cr_2), abs(n_cr_3), n_tor)

    # Return the critical load
    return nn_cr_flex_tor


def lmbda_flex(
    length,
    area,
    moi_y,
    kapa_bc=None,
    e_modulus=None,
    f_yield=None
):
    # Docstring
    """
    Flexural slenderness.

    Calculates the slenderness of a columne under pure compression.
    Euler's critical load is used.

    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    moi_y : float
        [mm^4] Moment of inertia
    kapa_bc : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    e_modulus : float, optional
        [MPa] Modulus of elasticity
        Default = 210000., typical value for steel
    f_yield : float, optional
        [MPa] yield stress.
        Default = 380., brcause this value was used extencively while the
        function was being written. To be changed to 235.

    Returns
    -------
    float
        [_] Member slenderness

    """
    # default values
    if kapa_bc is None:
        kapa_bc = 1.
    else:
        kapa_bc = float(kapa_bc)

    if e_modulus is None:
        e_modulus = 210000.
    else:
        e_modulus = float(e_modulus)

    if f_yield is None:
        f_yield = 380.
    else:
        f_yield = float(f_yield)

    # Calculate Euler's critical load
    n_cr = n_cr_flex(
        length,
        moi_y,
        e_modulus=e_modulus,
        kapa_bc=kapa_bc
    )

    # Flexural slenderness EN3-1-1 6.3.1.3 (1)
    lmbda_flexx = sqrt(area * f_yield / n_cr)

    # Return the result
    return lmbda_flexx


def imp_factor(b_curve):
    # Docstring
    """
    Imperfection factor.

    Returns the imperfection factor for a given buckling curve.
    The values are taken from Table 6.1 of EN1993-1-1 [1]

    Parameters
    ----------
    b_curve : {'a0', 'a', 'b', 'c', 'd'}
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].

    Returns
    -------
    float
        [_] Imperfection factor.

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings.
        Brussels: CEN, 2005.

    """
    switcher = {
        'a0': 0.13,
        'a': 0.21,
        'b': 0.34,
        'c': 0.49,
        'd': 0.76,
    }
    return switcher.get(b_curve, "nothing")


def chi_flex(
    length,
    area,
    moi_y,
    f_yield,
    b_curve,
    kapa_bc=None
):
    # Docstring
    """
    Flexural buckling reduction factor.

    Claculates the reduction factor, chi, according to EN1993-1-1 6.3.1.2

    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    moi_y : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are {'a0', 'a', 'b', 'c', 'd'}
    kapa_bc : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column

    Returns
    -------
    float
        [_] Reduction factor.

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings.
        Brussels: CEN, 2005.

    """
    if kapa_bc is None:
        kapa_bc = 1.

    lmda = lmbda_flex(
        length=length,
        area=area,
        moi_y=moi_y,
        kapa_bc=kapa_bc,
        e_modulus=None,
        f_yield=f_yield
    )

    alpha = imp_factor(b_curve)

    phi = (1 + alpha * (lmda - 0.2) + lmda ** 2) / 2.

    chi = 1 / (phi + sqrt(phi ** 2 - lmda ** 2))

    if chi > 1.:
        chi = 1.

    return chi


def n_b_rd(
    length,
    area,
    moi_y,
    f_yield,
    b_curve,
    kapa_bc=None,
    gamma_m1=None
):
    # Docstring
    """
    Flexural buckling resistance.

    Verifies the resistance of a column against flexural buckling
    according to EN1993-1-1 6.3.1.1.

    Parameters
    ----------
    length : float
        [mm] Column length
    area : float
        [mm^2] Cross section area
    moi_y : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are: {'a0', 'a', 'b', 'c', 'd'}
    kapa_bc : float, optional
        [_] Length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    gamma_m1 : float, optional
        [_] Partial safety factor.
        Default = 1.

    Returns
    -------
    float
        [N] Buckling resistance.

    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings.
        Brussels: CEN, 2005.

    """
    if kapa_bc is None:
        kapa_bc = 1.

    if gamma_m1 is None:
        gamma_m1 = 1.

    chi = chi_flex(length,
                   area,
                   moi_y,
                   f_yield,
                   b_curve,
                   kapa_bc=kapa_bc)

    nn_b_rd = area * f_yield * chi / gamma_m1

    return nn_b_rd


# CONNECTIONS
def bolt_grade2stress(bolt_grade):
    # Docstring
    """
    Convert bolt grade to yield and ultimate stress.

    Standard designation for bolt grade as a decimal is converted to yield and ultimate stress values in MPa. In the
    standard bolt grade designation, the integer part of the number represents the ultimate stress in MPa/100 and the
    decimal part is the yield stress as a percentage of the ultimate (e.g 4.6 is f_u = 400, f_y = 400 * 0.6 = 240).

    Parameters
    ----------
    bolt_grade : float

    Returns
    -------
    tuple : (f_ultimate, f_yield)
    """
    # Calculation using divmod
    f_ultimate = 100 * divmod(bolt_grade, 1)[0]
    f_yield = round(f_ultimate * divmod(bolt_grade, 1)[1])

    # Return values
    return f_ultimate, f_yield


def shear_area(bolt_size, shear_threaded=None):
    # Docstring
    """
    Shear area of a bolt.

    Returns the area to be used for the calculation of shear resistance of a bolt, either the gross cross-section of the
    bolt (circle area) or the reduced area of the threaded part of the bolt.

    Parameters
    ----------
    bolt_size : float
        Bolt's diameter.
    shear_threaded : bool, optional
        Designates if the shear plane is on the threaded portion or not.
        Default in False, which implies shearing of the non-threaded portion

    Returns
    -------
    float

    Notes
    -----
    Currently, the threaded area is based on an average reduction of the shank area. To be changed to analytic formula.
    """
    # Default
    if shear_threaded is None:
        shear_threaded = False

    # Calculate area
    if shear_threaded:
        a_shear = 0.784 * (pi * bolt_size ** 2 / 4)
    else:
        a_shear = pi * bolt_size ** 2 / 4

    # Return
    return a_shear


def f_v_rd(
           bolt_size,
           bolt_grade,
           shear_threaded=None,
           gamma_m2=None
           ):
    # Docstring
    """
    Bolt's shear resistance.

    Calculates the shear resistance of single bolt for one shear plane as given in table 3.4 of EC3-1-8.

    bolt_size : float
        Diameter of the non-threaded part (nominal bolt size e.g. M16 = 16)
    bolt_grade : float
        Bolt grade in standard designation format (see documentation of bolt_grade2stress())
    shear_threaded : bool, optional
        Designates if the shear plane is on the threaded portion or not.
        Default in False, which implies shearing of the non-threaded portion
    gamma_m2 : float, optional
        Safety factor.
        Default value is 1.25

    Returns
    -------
    float

    """

    # Defaults
    bolt_size = float(bolt_size)
    if shear_threaded is None:
        shear_threaded = False

    if gamma_m2 is None:
        gamma_m2 = 1.25
    else:
        gamma_m2 = float(gamma_m2)

    # av coefficient
    if shear_threaded and bolt_grade == (4.6 or 8.6):
        a_v = 0.5
    else:
        a_v = 0.6

    # Get ultimate stress for bolt
    f_ub = bolt_grade2stress(bolt_grade)[0]

    # Shear area
    a_shear = shear_area(bolt_size, shear_threaded)

    # Shear resistance
    ff_v_rd = a_v * f_ub * a_shear / gamma_m2

    # Return value
    return ff_v_rd
