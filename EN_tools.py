"""

A collection of structural steel related functions based on EN1993.

"""
from math import sqrt, pi, log, sin, cos, atan

# Define classes
# Cross-section class
class Cross_section:
    """
    Cross-section properties.
    
    Geometric and inertia parameters for structural purpose cross-section.
    
    Parameters
    ----------
    area : float
        [mm^2] Cross-sectional area.
    I_y : float
        [mm^4] Moment of inertia around y-axis.
        y-axis on the centre of gravity but not necessarily principal.
    I_z : float
        [mm^4] Moment of inertia around z-axis.
        z-axis on the centre of gravity but not necessarily principal.
    I_yz : float
        [mm^4] Product of inertia.
    I_torsion : float
        [mm^4] Saint Venant constant.
    I_warp : float
        [mm^6] Torsion constant.
    y_sc : float, optional
        [mm] Distance on y-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    z_sc : float, optional
        [mm] Distance on z-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    
    Attributes
    ----------
    
    Notes
    -----

    References
    ----------

    """
    def __init__(self, area, I_y, I_z):
        self.area = area
        self.I_y = I_y
        self.I_z = I_z


class Member:
    def __init__(self, profile, length, f_yield):
        self.profile = profile
        self.length = length
        self.f_yield = f_yield


###### SIMPLY SUPPORTED PLATE ######

def N_pl_Rd(
            thickness,
            width,
            f_yield,
            psi
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
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.
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

    # Critical compression load
    N_cr_plate = thickness * width * sigma_cr_plate(thickness, width, psi=psi)
    
    # Aeff calculation.
    # Reduction factor for the effective area of the profile acc. to EC3-1-5
    classification = width / (thickness * sqrt(235 / f_yield))
    lambda_p = classification / (28.4 * sqrt(k_sigma))
    if lambda_p > 0.673 and plate_class(width, thickness, f_yield) == 4:
        rho = (lambda_p - 0.055 * (3 + psi)) / lambda_p ** 2
    else:
        rho = 1.
    
    # Effective area
    A_eff = rho * thickness * width
    
    # Axial compression resistance , Npl
    N_pl_Rd = A_eff * f_yield
    
    # Return value
    return N_pl_Rd


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
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.
    
    """
    

    # Convert inputs to floats
    width, thickness, f_yield = float(width), float(thickness), float(f_yield)
    
    # Calculate classification
    classification = width / (thickness * sqrt(235. / f_yield))
    if  classification <= 33.:
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
    sigma_E =  190000 * (thickness / width) ** 2
    sigma_cr = sigma_E * k_sigma
    
    # Return value
    return sigma_cr


###### CYLINDRICAL SHELLS ######

def sigma_x_Rd(
               thickness,
               radius,
               length,
               f_y_k,
               fab_quality = None,
               gamma_M1 = None
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
    gamma_M1 : int, optional
        [_] Partial safety factor
        Default = 1.
    
    Returns
    -------
    float
        [MPa] Meridional buckling stress
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """

    # Default values
    if fab_quality is None:
        fab_quality = 'fcA'
    elif not((fab_quality is 'fcA') or (fab_quality is 'fcB') or (fab_quality is 'fcC')):
        print('Invalid fabrication class input. Choose between \'fcA\', \'fcB\' and \'fcC\' ')
    
    if gamma_M1 is None:
        gamma_M1 = 1.
    else:
        gamma_M1 = float(gamma_M1)
    
    # Fabrication quality class acc. to table D2
    if fab_quality is 'fcA':
        Q_factor = 40.
    elif fab_quality == 'fcB':
        Q_factor = 25.
    elif fab_quality == 'fcC':
        Q_factor = 16.
    
    # Critical meridinal stress, calculated on separate function
    sigma_cr = sigma_x_Rcr(thickness, radius, length)
    
    # Shell slenderness
    lmda = sqrt(f_y_k / sigma_cr[0])
    delta_w_k = (1. / Q_factor) * sqrt(radius / thickness) * thickness
    alpha = 0.62 / (1 + 1.91 * (delta_w_k / thickness) ** 1.44)
    beta = 0.6
    eta = 1.
    if sigma_cr[1] is 'long':
        # For long cylinders, a formula is suggested fo lambda, EC3-1-6 D1.2.2(4)
        # Currently, the general form is used. to be fixed.
        lmda_0 = 0.2
        #lmda_0 = 0.2 + 0.1 * (sigma_E_M / sigma_E)
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
    sigma_Rk = chi * f_y_k
    sigma_Rd = sigma_Rk / gamma_M1
    
    # Return value
    return sigma_Rd


def N_cr_shell(
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
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """

    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc to EN3-1-6 Annex D
    N_cr_shell = 2 * pi * radius * thickness * sigma_x_Rcr(thickness, radius, length)[0]
    
    # Return value
    return N_cr_shell


def sigma_x_Rcr(
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
    .. [1] Eurocode 3: Design of steel structures - Part 1-6: Strength and stability of shell structures. Brussels: CEN, 2006.

    """
    # Convert inputs to floats
    thickness, radius, length = float(thickness), float(radius), float(length)
    
    # Elastic critical load acc. to EN3-1-6 Annex D
    omega = length / sqrt(radius * thickness)
    if 1.7 <= omega and omega <= 0.5 * (radius / thickness):
        C_x = 1.
        length_category = 'medium'
    elif omega < 1.7:
        C_x = 1.36 - (1.83 / omega) + (2.07 / omega ** 2)
        length_category = 'short'
    else:
        # C_x_b is read on table D.1 of EN3-1-5 Annex D acc. to BCs
        # BC1 - BC1 is used on the Abaqus models (both ends clamped, see EN3-1-5 table 5.1)
        C_x_b = 6.
        C_x_N = max((1 + 0.2 * (1 - 2 * omega * thickness / radius) / C_x_b), 0.6)
        C_x = C_x_N
        length_category = 'long'
    
    # Calculate critical stress, eq. D.2 on EN3-1-5 D.1.2.1-5
    sigma_cr = 0.605 * 210000 * C_x * thickness / radius
    
    # Return value
    return sigma_cr, length_category


###### OVERALL BUCKLING ######
def N_cr_flex(
    length,
    I_y,
    kapa_BC = None,
    E_modulus = None
    ):
    
    # Docstring
    """
    Euler's critical load.
    
    Calculates the critical load for flexural buckling of a given column.
    A single direction is considered. If more directions are required
    (e.g the two principal axes), the function has to be called multiple
    times. For torsional mode critical load use N_cr_tor(), and for
    flexural-torsional critical load use N_cr_flex_tor()
        
    Parameters
    ----------
    length : float
        [mm] Column length.
    I_y : float
        [mm^4] Moment of inertia.
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column.
    E_modulus : float, optional
        [MPa] Modulus of elasticity.
        Default = 210000., typical value for steel.
    
    Returns
    -------
    float
        [N] Critical load.
    
    """
    # default values
    if kapa_BC is None:
        kapa_BC = 1.
    else:
        kapa_BC = float(kapa_BC)
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    # Euler's critical load
    N_cr_flex = (pi ** 2) * E_modulus * I_y / (kapa_BC * length) ** 2
    
    # Return the result
    return N_cr_flex


def N_cr_tor(
    length,
    area,
    I_y0,
    I_z0,
    I_torsion,
    I_warp,
    y_0 = None,
    z_0 = None,
    E_modulus = None,
    poisson = None,
    ):
    
    # Docstring
    """
    Torsional elastic critical load
    
    Calculates the torsional elastic critical load for a hinged column.
    The input values are refering to the principal axes. For flexural
    buckling (Euler cases) use N_cr_flex. For the combined 
    flexural-torsional modes use N_cr_flex_tor.
    
    Parameters
    ----------
    length : float
        [mm] Column length.
    area : float
        [mm^2] Cross-sectional area.
    I_y : float
        [mm^4] Moment of inertia around `y`-axis.
        `y`-axis on the centre of gravity but not necessarily principal.
    I_z : float
        [mm^4] Moment of inertia around `z`-axis.
        `z`-axis on the centre of gravity but not necessarily principal.
    I_torsion : float
        [mm^4] Saint Venant constant.
    I_warp : float
        [mm^6] Torsion constant.
    y_0 : float, optional
        [mm] Distance on `y`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    z_0 : float, optional
        [mm] Distance on `z`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    E_modulus : float, optional
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
        
        :math:`r^2=(I_y + I_z)/A + x_0^2 + y_0^2` 
        
        :math:`x_0, y_0`  : Shear centre coordinates on the principal 
        coordinate system
        

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
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    if poisson is None:
        poisson = 0.3
    else:
        poisson = float(poisson)
    
    # Shear modulus
    G_modulus = E_modulus / (2 * (1 + poisson))
    
    # Polar radius of gyration.
    i_pol = sqrt((I_y0 + I_z0) / area)
    i_zero = sqrt(i_pol ** 2 + y_0 ** 2 + z_0 ** 2)
    
    # Calculation of critical torsional load.
    N_cr_tor = (1 / i_zero ** 2) * (G_modulus * I_torsion + (pi ** 2 * E_modulus * I_warp / length ** 2))
    
    # Return the result
    return N_cr_tor
    

def N_cr_flex_tor(
    length,
    area,
    I_y,
    I_z,
    I_yz,
    I_torsion,
    I_warp,
    y_sc = None,
    z_sc = None,
    E_modulus = None,
    poisson = None,
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
    I_y : float
        [mm^4] Moment of inertia around `y`-axis.
        `y`-axis on the centre of gravity but not necessarily principal.
    I_z : float
        [mm^4] Moment of inertia around `z`-axis.
        `z`-axis on the centre of gravity but not necessarily principal.
    I_yz : float
        [mm^4] Product of inertia.
    I_torsion : float
        [mm^4] Saint Venant constant.
    I_warp : float
        [mm^6] Torsion constant.
    y_sc : float, optional
        [mm] Distance on `y`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    z_sc : float, optional
        [mm] Distance on `z`-axis of the shear center to the origin.
        Default = 0, which implies symmetric profile
    E_modulus : float, optional
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
    'N_cr_flex').

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
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    if poisson is None:
        poisson = 0.3
    else:
        poisson = float(poisson)
    
    # Shear modulus
    G_modulus = E_modulus / (2 * (1 + poisson))
        
    # Angle of principal axes
    if abs(I_y - I_z) < 1e-20:
        theta = pi/4
    else:
        theta = -atan((2 * I_yz) / (I_y - I_z)) / 2
    
    # Distance of the rotation centre to the gravity centre on the
    # principal axes coordinate system
    y_0 = y_sc * cos(-theta) - z_sc * sin(-theta)
    z_0 = z_sc * cos(-theta) + y_sc * sin(-theta)
    
    # Moment of inertia around principal axes.
    I_y0 = (I_y + I_z)/2 + sqrt(((I_y - I_z) / 2)**2 + (I_yz)**2)
    I_z0 = (I_y + I_z)/2 - sqrt(((I_y - I_z) / 2)**2 + (I_yz)**2)
    
    # Polar radius of gyration.
    i_pol = sqrt((I_y0 + I_z0) / area)
    i_zero = sqrt(i_pol ** 2 + y_0 ** 2 + z_0 ** 2)
    
    # Independent critical loads for flexural and torsional modes.
    N_cr_max = (pi ** 2 * E_modulus * I_y0) / (length ** 2)
    N_cr_min = (pi ** 2 * E_modulus * I_z0) / (length ** 2)
    N_tor = N_cr_tor(
        length,
        area,
        I_y0,
        I_z0,
        I_torsion,
        I_warp = I_warp,
        y_0 = y_0,
        z_0 = z_0,
        E_modulus = E_modulus,
        poisson = poisson
        )
    
    # Coefficients of the 3rd order equation for the critical loads
    # The equation is in the form aaaa * N ^ 3 - bbbb * N ^ 2 + cccc * N - dddd
    aaaa = i_zero ** 2 - y_0 ** 2 - z_0 ** 2
    bbbb = ((N_cr_max + N_cr_min + N_tor) * i_zero ** 2) - (N_cr_min * y_0 ** 2) - (N_cr_max * z_0 ** 2)
    cccc = i_zero ** 2 * (N_cr_min * N_cr_max) + (N_cr_min * N_tor) + (N_tor * N_cr_max)
    dddd = i_zero ** 2 * N_cr_min * N_cr_max * N_tor
    
    DD = (4 * (-bbbb ** 2 + 3 * aaaa * cccc) ** 3 + (2 * bbbb ** 3 - 9 * aaaa * bbbb * cccc + 27 * aaaa ** 2 * dddd) ** 2)
    if (DD < 0):
       DD = -1. * DD
       cf = 1j
    else:
       cf = 1
    
    # Critical load
    # The following N_cr formulas are the roots of the 3rd order equation of the global critical load
    N_cr_1 = bbbb/(3.*aaaa) - (2**(1./3)*(-bbbb**2 + 3*aaaa*cccc))/                     \
        (3.*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                      \
        (cf*sqrt(DD)))**(1./3)) + (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +     \
        (cf*sqrt(DD)))**(1./3)/(3.*2**(1./3)*aaaa)
    
    N_cr_2 = bbbb/(3.*aaaa) + ((1 + (0 + 1j)*sqrt(3))*(-bbbb**2 + 3*aaaa*cccc))/        \
        (3.*2**(2./3)*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +            \
        (cf*sqrt(DD)))**(1./3)) - ((1 - (0 + 1j)*sqrt(3))*                              \
        (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                               \
        (cf*sqrt(DD)))**(1./3))/(6.*2**(1./3)*aaaa)
    
    N_cr_3 = bbbb/(3.*aaaa) + ((1 - (0 + 1j)*sqrt(3))*(-bbbb**2 + 3*aaaa*cccc))/        \
        (3.*2**(2./3)*aaaa*(2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +            \
        (cf*sqrt(DD)))**(1./3)) - ((1 + (0 + 1j)*sqrt(3))*                              \
        (2*bbbb**3 - 9*aaaa*bbbb*cccc + 27*aaaa**2*dddd +                               \
        (cf*sqrt(DD)))**(1./3))/(6.*2**(1./3)*aaaa)
    
    # Lowest root is the critical load
    N_cr_flex_tor = min(abs(N_cr_1), abs(N_cr_2), abs(N_cr_3), N_tor)
    
    # Return the critical load
    return N_cr_flex_tor


def lmbda_flex(
    length,
    area,
    I_y,
    kapa_BC = None,
    E_modulus = None,
    f_yield = None
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
    I_y : float
        [mm^4] Moment of inertia
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    E_modulus : float, optional
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
    if kapa_BC is None:
        kapa_BC = 1.
    else:
        kapa_BC = float(kapa_BC)
    
    if E_modulus is None:
        E_modulus = 210000.
    else:
        E_modulus = float(E_modulus)
    
    if f_yield is None:
        f_yield = 380.
    else:
        f_yield = float(f_yield)
    
    # Calculate Euler's critical load
    N_cr = N_cr_flex(
        length,
        I_y,
        E_modulus = E_modulus,
        kapa_BC = kapa_BC
        )
    
    # Flexural slenderness EN3-1-1 6.3.1.3 (1)
    lmbda_flex = sqrt(area * f_yield / N_cr)
    
    # Return the result
    return lmbda_flex


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
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

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
             I_y,
             f_yield,
             b_curve,
             kapa_BC = None
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
    I_y : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are {'a0', 'a', 'b', 'c', 'd'}
    kapa_BC : float, optional
        [_] length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    
    Returns
    -------
    float
        [_] Reduction factor.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

    """
    if kapa_BC is None:
        kapa_BC = 1.
    
    lmda = lmbda_flex(
        length = length,
        area = area,
        I_y = I_y,
        kapa_BC = kapa_BC,
        E_modulus = None,
        f_yield = f_yield
        )
    
    alpha = imp_factor(b_curve)
    
    phi = (1 + alpha * (lmda - 0.2) + lmda**2) / 2.
    
    chi = 1 / (phi + sqrt(phi**2 - lmda**2))
    
    if chi > 1.:
        chi = 1.
    
    return chi


def N_b_Rd(
           length, 
           area, 
           I_y, 
           f_yield, 
           b_curve, 
           kapa_BC = None, 
           gamma_M1 = None
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
    I_y : float
        [mm^4] Moment of inertia
    f_yield : float
        [MPa] Yield stress.
    b_curve : str
        [_] Name of the buckling curve as obtained from Table 6.2 of [1].
        Valid options are: {'a0', 'a', 'b', 'c', 'd'}
    kapa_BC : float, optional
        [_] Length correction for the effect of the boundary conditions.
        Default = 1, which implies simply supported column
    gamma_M1 : float, optional
        [_] Partial safety factor.
        Default = 1.
    
    Returns
    -------
    float
        [N] Buckling resistance.
    
    References
    ----------
    .. [1] Eurocode 3: Design of steel structures - Part 1-1: General rules and rules for buildings. Brussels: CEN, 2005.

    """
    if kapa_BC is None:
        kapa_BC = 1.
    
    if gamma_M1 is None:
        gamma_M1 = 1.
    
    chi = chi_flex(length,
        area,
        I_y,
        f_yield,
        b_curve,
        kapa_BC=kapa_BC)
        
    N_b_Rd = area * f_yield * chi / gamma_M1
    
    return N_b_Rd


###### CIRCULAR PLATE  #####
# The following formulas are taken from table 18-3 of:
# W. D. Pilkey, Formulas for stress, strain, and structural matrices. New York: Wiley, 1994.
#
## Hinged perimeter
### Concentrated force applied on circle
#############(not finished, there is an error, to be revised)###########
#
#def circular_plate(a_L, a_1, W_load, r):
#    poisson = 0.3
#    
#    alfa = r / a_L
#    beta = a_1 / a_L
#    
#    C_1 = (3 + poisson) / (1 + poisson) * (1 + beta ** 2) + 2 * beta ** 2 * log(beta)
#    C_2 = (1 - poisson) / (1 + poisson) * (1 - beta ** 2) + 2 * log(beta)
#    C_3 = (3 + poisson) / (1 + poisson) - beta ** 2 * (1 - poisson) / (1 + poisson)
#    
#    if alfa <= beta:
#        M_r = (1 / 4) * W_load * a_L * (1 + poisson) * C_2
#        M_phi = M_r
#    else:
#        M_r = (1 / 4) * (W_load * a_1) / (alfa ** 2) * ((1 - poisson) * (1 - alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
#        M_phi = (1 / 4) * (W_load * a_1) / (alfa ** 2) * (2 * (1 - poisson) * alfa ** 2 - (1 - poisson) * (1 + alfa ** 2) * beta ** 2 - 2 * (1 + poisson) * alfa ** 2 * log(alfa))
#    
#    return M_r, M_phi