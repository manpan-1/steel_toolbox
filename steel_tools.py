# Engineering and structural steel related functions
'''
'''

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


# Calculate diameter of washer for a given bolt
def bolt2washer(M_bolt):
    d_washer = math.ceil(1.5893*M_bolt+5.1071)
    return d_washer
