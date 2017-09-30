# Extra python toolset for Abaqus
# version: 1
# Contains:
# open_odb
# max_result

# code to be evaluated for max field output value
# https://gist.github.com/crmccreary/1074551
import sys
import os
import odbAccess
import numpy as np
import math
import string
from abaqusConstants import *

# Retudn a list with all the divisors of a numbers
def divisorGenerator(n):
    large_divisors = []
    for i in xrange(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor


# mean value
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)


# Calculate diameter of washer for a given bolt
def bolt2washer(M_bolt):
    d_washer = math.ceil(1.5893*M_bolt+5.1071)
    return d_washer


# definition of a method to search for a keyword position
def GetBlockPosition(model,blockPrefix):
    import string
    pos = 0
    for block in model.keywordBlock.sieBlocks:
        if string.lower(block[0:len(blockPrefix)])==string.lower(blockPrefix):
            return pos
        pos=pos+1
    return -1


# A more sophisticated open odb function
def open_odb(odbPath):
    base, ext = os.path.splitext(odbPath)
    odbPath = base + '.odb'
    new_odbPath = None
    if odbAccess.isUpgradeRequiredForOdb(upgradeRequiredOdbPath=odbPath):
        print('odb %s needs upgrading' % (odbPath,))
        path,file_name = os.path.split(odbPath)
        file_name = base + "_upgraded.odb"
        new_odbPath = os.path.join(path,file_name)
        odbAccess.upgradeOdb(existingOdbPath=odbPath, upgradedOdbPath=new_odbPath)
        odbPath = new_odbPath
    odb = odbAccess.openOdb(path=odbPath, readOnly=True)
    return odb


# Calculate cross sectional properties. Two inputs required:
# A 2d (2, n) list with x,y values of nodes and a (3, m) 2d list for elements(first-node, second-node, thickness)
def cs_prop(nodes, elem):
    # Calculate cross sectional properties
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nele = len(elem[0])
    node = elem[0]+elem[1]
    nnode = 0
    j = 0
    
    while node:
        i = [ii for ii,x in enumerate(node) if x==node[0]]
        for ii in sorted(i, reverse=True):
            del node[ii]
        if len(i)==2:
            j += 1
        nnode += 1
    
    # classify the section type
    if j == nele:
        section = 'close' #single cell
    elif j == nele-1:
        section = 'open' #singly-branched
    else:
        section = 'open' #multi-branched

    # Calculate the cs-properties
    tt = []
    xm = []
    ym = []
    xd = []
    yd = []
    L = []
    for i in range(nele):
        sn = elem[0][i]
        fn = elem[1][i]
        # thickness of the element
        tt = tt+[elem[2][i]]
        # compute the coordinate of the mid point of the element
        xm = xm + [mean([nodes[0][sn], nodes[0][fn]])]
        ym = ym + [mean([nodes[1][sn], nodes[1][fn]])]
        # compute the dimension of the element
        xd = xd + [(nodes[0][fn]-nodes[0][sn])]
        yd = yd + [(nodes[1][fn]-nodes[1][sn])]
        # compute the length of the element
        L = L + [math.sqrt(xd[i]**2+yd[i]**2)]

    # calculate cross sectional area
    A = sum([a*b for a,b in zip(L, tt)])
    # compute the centroid 
    xc = sum([a*b*c for a,b,c in zip(L,tt,xm)])/A
    yc = sum([a*b*c for a,b,c in zip(L,tt,ym)])/A
    
    if abs(xc/math.sqrt(A)) < 1e-12:
        xc = 0
    
    if abs(yc/math.sqrt(A)) < 1e-12:
        yc = 0
    
    # Calculate MOI
    Ix = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(yd,L,tt)], [(a-yc)**2*b*c for a,b,c in zip(ym,L,tt)])])
    Iy = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(xd,L,tt)], [(a-xc)**2*b*c for a,b,c in zip(xm,L,tt)])])
    Ixy = sum([sum(a) for a in zip([a*b*c*d/12 for a,b,c,d in zip(xd,yd,L,tt)], [(a-xc)*(b-yc)*c*d for a,b,c,d in zip(xm,ym,L,tt)])])
    
    if abs(Ixy/A**2) < 1e-12:
        Ixy = 0
    
    # Calculate angle of principal axes
    if Ix == Iy:
        theta_principal = math.pi/2
    else:
        theta_principal = math.atan((-2*Ixy)/(Ix-Iy))/2
    
    # Change to centroid principal coordinates
    coord12 = [[a-xc for a in nodes[0]],[a-yc for a in nodes[1]]]
    coord12 = np.array([[math.cos(theta_principal), math.sin(theta_principal)],[-math.sin(theta_principal), math.cos(theta_principal)]]).dot(nodes)
    
    # re-calculate cross sectional properties for the centroid 
    for i in range(nele):
        sn = elem[0][i]
        fn = elem[1][i]
        # calculate the coordinate of the mid point of the element
        xm = xm + [mean([coord12[0][sn], coord12[0][fn]])]
        ym = ym + [mean([coord12[1][sn], coord12[1][fn]])]
        # calculate the dimension of the element
        xd = xd + [(coord12[0][fn]-coord12[0][sn])]
        yd = yd + [(coord12[1][fn]-coord12[1][sn])]
    
    # calculate the principal moment of inertia
    I1 = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(yd,L,tt)], [(a-yc)**2*b*c for a,b,c in zip(ym,L,tt)])])
    I2 = sum([sum(a) for a in zip([a**2*b*c/12 for a,b,c in zip(xd,L,tt)], [(a-xc)**2*b*c for a,b,c in zip(xm,L,tt)])])
    
    # Return values
    return A, xc, yc, Ix, Iy, Ixy, I1, I2, theta_principal


# Calculate xy of nodes for a given polygonal profile
# Returns points for the entire profile (1st and 2nd returned values)
# and points for a single sector (3rd and 4th returned values)
def polygon_sector(n, R, t, tg, rbend, nbend, l_lip):
    # Angle corresponding to one face of the polygon
    theta = 2*math.pi/n;
    
    # Angles of radii (measured from x-axis)
    phi=np.linspace(5*math.pi/6, math.pi/6, n/3+1)
    
    # xy coords of the polygon's corners
    x = R*np.cos(phi);
    y = R*np.sin(phi);
    
    ## Bends
    
    # Distance between bending centre and corner
    lc = rbend/np.cos(theta/2)
    
    # Centers of bending arcs
    xc  = x[1:-1] - lc*np.cos(phi[1:-1])
    yc  = y[1:-1] - lc*np.sin(phi[1:-1])
    
    # Bending arc angle
    theta_b = math.pi - theta
    
    # Angles of the edges' midlines (measured from x-axis)
    phi_mids = phi[0:-1] - theta/2 
    
    # xy coords of the arc's points
    xarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))]
    yarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))] 
    for i in range(int(n/3 -1)):
        for j in range(nbend+1):
            xarc[i][j] = xc[i] + rbend*np.cos(phi_mids[i]-(j)*(theta/nbend))
            yarc[i][j] = yc[i] + rbend*np.sin(phi_mids[i]-(j)*(theta/nbend))
    
    ## Start-end extensions
    # Bending radius
    rs = rbend/2
    xcs = [0, 0]
    ycs = [0, 0]
    
    # First bend
    v1 = phi_mids[0]-math.pi/2
    v2 = (phi[0]+phi_mids[0]-math.pi/2)/2
    l1 = (t+tg)/(2*np.cos(phi[0]-phi_mids[0]))
    l2 = rs/np.sin(v2-phi_mids[0]+math.pi/2)
    x1 = x[0]+l1*np.cos(v1)
    y1 = y[0]+l1*np.sin(v1)
    
    # First bend centre coords
    xcs[0] = x1+l2*np.cos(v2)
    ycs[0] = y1+l2*np.sin(v2)
    
    # Last bend
    v1 = phi_mids[-1]+math.pi/2
    v2 = (v1+phi[-1])/2
    l1 = (t+tg)/(2*np.cos(v1-phi[-1]-math.pi/2))
    l2 = rs/np.sin(v2-phi[-1])
    x1 = x[-1]+l1*np.cos(v1)
    y1 = y[-1]+l1*np.sin(v1)
    
    # Last bend centre coords
    xcs[1] = x1+l2*np.cos(v2)
    ycs[1] = y1+l2*np.sin(v2)
    
    # First and last bend arc points coords
    xsarc = [[0 for j in range(nbend+1)] for j in [0,1]]
    ysarc = [[0 for j in range(nbend+1)] for j in [0,1]] 
    for j in range(nbend+1):
        xsarc[0][j] = xcs[0] + rs*np.cos(4*math.pi/3+(j)*((phi_mids[0]-math.pi/3)/nbend))
        ysarc[0][j] = ycs[0] + rs*np.sin(4*math.pi/3+(j)*((phi_mids[0]-math.pi/3)/nbend))
        xsarc[1][j] = xcs[1] + rs*np.cos(phi_mids[-1]+math.pi+(j)*((phi[-1]+math.pi/2-phi_mids[-1])/nbend))
        ysarc[1][j] = ycs[1] + rs*np.sin(phi_mids[-1]+math.pi+(j)*((phi[-1]+math.pi/2-phi_mids[-1])/nbend))
    
    
    ## Points of the lips
    
    # Lip length according to bolt washer diameter
    
    # First lip
    xstart = [xsarc[0][0] + l_lip*np.cos(phi[0]), xsarc[0][0] + l_lip*np.cos(phi[0])/2]
    ystart = [ysarc[0][0] + l_lip*np.sin(phi[0]), ysarc[0][0] + l_lip*np.sin(phi[0])/2]
    
    
    # Last point
    xend = [xsarc[1][-1] + l_lip*np.cos(phi[-1])/2, xsarc[1][-1] + l_lip*np.cos(phi[-1])]
    yend = [ysarc[1][-1] + l_lip*np.sin(phi[-1])/2, ysarc[1][-1] + l_lip*np.sin(phi[-1])]
    
    ## Collect the x, y values in a sorted 2xn array
    xarcs, yarcs=[],[]
    for i in range(len(phi)-2):
        xarcs=xarcs+xarc[i][:]
        yarcs=yarcs+yarc[i][:]
    
    x_sector = xstart+xsarc[0][:]+xarcs[:]+xsarc[1][:]+xend
    y_sector = ystart+ysarc[0][:]+yarcs[:]+ysarc[1][:]+yend
    
    # Copy-rotate the points of the first sector to create the entire CS
    # Rotation matrix
    Rmat = np.array([[math.cos(-2*math.pi/3), -math.sin(-2*math.pi/3)], [math.sin(-2*math.pi/3), math.cos(-2*math.pi/3)]])
    
    # Dot multiply matrices
    coord1 = np.array([x_sector, y_sector])
    coord2 = Rmat.dot(coord1)
    coord3 = Rmat.dot(coord2)
    
    # Concatenate into a single xy array
    x_cs = np.concatenate([coord1[0], coord2[0], coord3[0]])
    y_cs = np.concatenate([coord1[1], coord2[1], coord3[1]])
    
    # Return matrices
    return x_cs, y_cs, x_sector, y_sector


# Look for the max value in a field output
# Function name to be changed to field_max(). (look for instances on other script e.g. semi-closed script)
def field_max(odb, result):
    result_field, result_invariant = result
    _max = -1.0e20
    for step in odb.steps.values():
        print 'Processing Step:', step.name
        for frame in step.frames:
            if frame.frameValue > 0.0:
                allFields = frame.fieldOutputs
                if (allFields.has_key(result_field)):
                    stressSet = allFields[result_field]
                    for stressValue in stressSet.values:
                        if result_invariant:
                            if hasattr(stressValue, result_invariant.lower()):
                                val = getattr(stressValue,result_invariant.lower())
                            else:
                                raise ValueError('Field value does not have invariant %s' % (result_invariant,))
                        else:
                            val = stressValue.data
                        if ( val > _max):
                            _max = val
                else:
                    raise ValueError('Field output does not have field %s' % (results_field,))
    return _max


# Look for the max value in a history output
# TO BE FIXED. THE REFERENCE POINTS (rp1key, ho1key etc.) ARE NOT GENERIC.
# Fetch maximum load, displacement and LPF for a riks analysis.
# The method assumes that a) the the odb is located in the current directory
def history_max(odb_name, step_name):
    myOdb = odbAccess.openOdb(path=odb_name+'.odb')
    RIKSstep = myOdb.steps[step_name]
    rp1key = RIKSstep.historyRegions.keys()[1]
    ho1key = RIKSstep.historyRegions[rp1key].historyOutputs.keys()[0]
    rp2key = RIKSstep.historyRegions.keys()[2]
    ho2key = RIKSstep.historyRegions[rp2key].historyOutputs.keys()[0]
    asskey = RIKSstep.historyRegions.keys()[0]
    hoasse = RIKSstep.historyRegions[asskey].historyOutputs.keys()[-1]
    load_hist = RIKSstep.historyRegions[rp1key].historyOutputs[ho1key].data
    disp_hist = RIKSstep.historyRegions[rp2key].historyOutputs[ho2key].data
    lpf_hist = RIKSstep.historyRegions[asskey].historyOutputs[hoasse].data
    maxpos = load_hist.index(max(load_hist,key=lambda x:x[1]))
    load = load_hist[maxpos][1]
    disp = -disp_hist[maxpos][1]
    lpf = lpf_hist[maxpos][1]
    odbAccess.closeOdb(myOdb)
    return lpf, load, disp


def fetch_eigenv(odb_name, step_name, n_eigen):
    bckl_odb = odbAccess.openOdb(path=odb_name+'.odb')
    bckl_step = bckl_odb.steps[step_name]
    
    # Gather the eigenvalues
    eigenvalues = ()
    eigen_string = ""
    for J_eigenvalues in range(1, n_eigen + 1):
        current_eigen = float(bckl_step.frames[J_eigenvalues].description[-11:])
        eigenvalues = eigenvalues + (current_eigen, )
        eigen_string = eigen_string + "%.3E "%(current_eigen)
    
    # Close the odb
    odbAccess.closeOdb(bckl_odb)
    
    # Return variables
    return eigenvalues, eigen_string


def plastic_table(nominal = None):
    
    if nominal is None:
        nominal = 'S355'
    else:
        nominal = str(nominal)
    
    if nominal is 'S355':
        table=(
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
        table=(
            (760., 0.0),
            (770., 0.022),
            (850., 0.075),
            (900., 0.1),
            (901., 1.)
            )
    
    return table