"""
Author: Yaset Acevedo
Created: 10/10/14
Update Date: 11/9/15
Version: 2.1
"""

#Import modern print module
from __future__ import print_function
from __future__ import division
import warnings
import math

class dreidingAtom:
    'Common base class for all atoms'
    # Initialize global class atom variables
    atomCount = 0
    # Format for elementList: [['B', atomID1, atomID2], ['H',atomID1, atomID2]]
    # Initialize list of elements to output in a particular order.
    # Useful for visualization programs
    elementList = [['H'], ['B'], ['C'], ['N'], ['O'], ['Ni'], ['CB']]
    #elementList.append()
    
    bondCount = 0
    globalBondType = []
    
    angleCount = 0
    angleType = {}
    
    dihedralCount = 0
    dihedralList = []
    
    improperCount = 0
    
    # Initialize box and boundary variables
    periodic = False
    xlo = None
    xhi = None
    ylo = None
    yhi = None
    zlo = None
    zhi = None
    xy = None
    xz = None
    yz = None

    def __init__(self, element, x, y, z):
        # Record element and x,y,z information
        self.element = element
        self.elementType = element
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        dreidingAtom.atomCount = dreidingAtom.atomCount + 1
        self.atomID = dreidingAtom.atomCount
        
        # Create dummy variable for possible bonds, bond orders, bond types
        self.bondList = []
        self.bondOrderList = []
        self.bondTypeList = []
        # Create dummy variable for possible angles
        self.angleList = []
        
        # Set initial charge to zero
        self.charge = 0
        
        # Assign covalent radius of element in Angstrom
        # More specific radius will be used once valency has been determined
        if self.element == 'H':
            self.radius = 0.420
        elif self.element == 'B':
            self.radius = 0.800
        elif self.element == 'C':
            self.radius = 0.770
        elif self.element == 'N':
            self.radius = 0.702
        elif self.element == 'O':
            self.radius = 0.800
        elif self.element == 'Ni':
            self.radius = 1.3
        else:
            # Will not bond if given a huge radius
            warnings.warn("Unknown element type")
            self.radius = 1000
            
        # Assign max number of bonds
        if self.element == 'H':
            self.maxBond = 1
        elif self.element == 'B':
            self.maxBond = 3
        elif self.element == 'C':
            self.maxBond = 3
        elif self.element == 'N':
            self.maxBond = 3
        elif self.element == 'O':
            self.maxBond = 2
        elif self.element == 'Ni':
            self.maxBond = 4
        else:
            # Will not bond if maxBond is 0
            warnings.warn("Unknown element type")
            self.maxBond = 0
            
        # Organize atoms into lists of same elements
        for (offset, group) in enumerate(dreidingAtom.elementList):
            if group[0] == self.getElement():
                dreidingAtom.elementList[offset].append(self.getID())
                break
        else:
            dreidingAtom.elementList.append([self.getElement()])
            dreidingAtom.elementList[-1].append(self.getID())
            
    def shiftXYZ(self, shiftx, shifty, shiftz):
        #Move the atom by the prescribed amount
        self.x = self.x + shiftx
        self.y = self.y + shifty
        self.z = self.z + shiftz

    # Requires periodic box information
    def wrapXYZ(self):
        # Create copy of xyz coordinates
        x1 = self.getX()
        y1 = self.getY()
        z1 = self.getZ()
        temp_pos = {'x':x1, 'y':y1, 'z':z1}

        # Calculate z multiple, y multiple, then x multiple
        # z multiple
        z_multiple = temp_pos['z'] / (dreidingAtom.zhi - dreidingAtom.zlo)
        temp_yz = dreidingAtom.yz * z_multiple
        temp_xz = dreidingAtom.xz * z_multiple

        # Redefine temporary position
        temp_pos['x'] = temp_pos['x'] - dreidingAtom.xz * z_multiple
        temp_pos['y'] = temp_pos['y'] - dreidingAtom.yz * z_multiple
        temp_pos['z'] = temp_pos['z'] - (dreidingAtom.zhi - dreidingAtom.zlo) * z_multiple

        # y multiple
        y_multiple = temp_pos['y'] / (dreidingAtom.yhi - dreidingAtom.ylo)
        temp_xy = dreidingAtom.xy * y_multiple

        # Redefine temporary position
        temp_pos['x'] = temp_pos['x'] - dreidingAtom.xy * y_multiple
        temp_pos['y'] = temp_pos['y'] - (dreidingAtom.yhi - dreidingAtom.ylo) * y_multiple

        # z multiple
        x_multiple = temp_pos['x'] / (dreidingAtom.xhi - dreidingAtom.xlo)

        # Redefine temporary position
        temp_pos['x'] = temp_pos['x'] - (dreidingAtom.xhi - dreidingAtom.xlo) * x_multiple

        # Output final temporary position. Should be (0,0,0)
        #text = '%.5f, %.5f. %.5f' % (temp_pos['x'], temp_pos['y'], temp_pos['z'])
        #print(text)

        # Output scaled coordinates
        #text = '%.5f, %.5f. %.5f' % (x_multiple, y_multiple, z_multiple)
        #print(text)

        # Place atom into center periodic box
        if x_multiple < 0:
            x_multiple = x_multiple + 1
        elif x_multiple > 1:
            x_multiple = x_multiple - 1
        if y_multiple < 0:
            y_multiple = y_multiple + 1
        elif y_multiple > 1:
            y_multiple = y_multiple - 1
        if z_multiple < 0:
            z_multiple = z_multiple + 1
        elif z_multiple > 1:
            z_multiple = z_multiple - 1

        # Translate atom and assign to atom
        temp_pos = translateVector1A(temp_pos, x_multiple)
        temp_pos = translateVector2B(temp_pos, y_multiple)
        temp_pos = translateVector3C(temp_pos, z_multiple)

        self.x = temp_pos['x']
        self.y = temp_pos['y']
        self.z = temp_pos['z']

        # Output final position. Should be (0,0,0)
        #self.displayAtom()
            
    def assignBond(self, atom2):
        # Get id of atom being bonded to
        atom2ID = atom2.getID()
        
        # Make sure that atom has not already been bonded to. Add bond if 
        # it has not already been bonded.
        if self.bondList.count(atom2ID) == 0:
            self.bondList.append(atom2ID)
            
            # Always make initial guess of bond order 1
            self.bondOrderList.append(1)
            # Always make initial assumption of only one bond type
            self.bondTypeList.append(1)
            
            # Return boolean true if bond was successful
            return True
            
        else:
            # Return boolean false if bond was not successful
            return False
            
    def breakAllBonds(self):
        # Get list of bonds
        bondList = self.getBondList()        
        
        # Cycle through all bonds
        for bondID in bondList:
            
            # Get bond atom object and element type
            bondAtom = atomList[bondID-1]
            
            # Break bond for self and other atom
            self.breakBond(bondAtom)
            bondAtom.breakBond(self)
                
    def breakBond(self, atom2):
        # Must be done before bond order and bond type is set
        # Get id of atom being bonded to
        atom2ID = atom2.getID()
        
        # Remove specified atom
        self.bondList.remove(atom2ID)
        
        # Remove bond order guess and bond type guess
        # Must be done before bond order and bond type is set
        self.bondOrderList.pop()
        self.bondTypeList.pop()
            
    def assignElementType(self):
        # For uniformity with elements that might be two characters long, add 
        # an underscore to the elements with one character
        if len(self.elementType) == 1:
            self.elementType = self.elementType + '_'
            
        # Find number of electron pairs from bonding
        bondElectronPairs = len(self.getBondList())            
        
        # Count number of bonds and determine hybridization
        if bondElectronPairs == 4:
            self.elementType = self.elementType + '3'
        elif bondElectronPairs == 3:
            # Check for the special case where carbon/nitrogen is in a resonance 
            # structure
            if self.elementType == 'C_' or self.elementType == 'N_':
                self.elementType = self.elementType + 'R'
            else:
                self.elementType = self.elementType + '2'
        elif bondElectronPairs == 2:
            # Check for the special case where there are extra electron pairs
            if self.elementType == 'O_':
                self.elementType = self.elementType + '3'
            # Check for the special case where nitrogen is in a resonance 
            # structure
            elif self.elementType == 'N_':
                self.elementType = self.elementType + 'R'
            else:
                self.elementType = self.elementType + '1'
        elif bondElectronPairs == 1:
            # Check for the special case where there are extra electron pairs
            if self.elementType == 'O_':
                self.elementType = self.elementType + '2'
            else:
                self.elementType = self.elementType + '1'
        else:
            if bondElectronPairs == 0:
                self.elementType = self.elementType + str(len(self.getBondList()))
            else:
                self.elementType = self.elementType + str(len(self.getBondList())-1)
            #print(self.elementType)
            #self.displayAtom()
            
        return self.elementType
        
    def assignAngles(self, atomList):
        # Generate internal list of angles if applicable
        # In addition, add to global list of types and assign angle type
        if len(self.getBondList()) > 1:
            # See if angle type has already been added to the angleType list
            if (self.getElementType() in dreidingAtom.angleType) == False:
                # If angle has not been added, add to the dictionary with new 
                # number
                n = len(dreidingAtom.angleType)
                dreidingAtom.angleType[self.getElementType()] = n+1
            
            # Select appropriate angleType number for next step
            n = dreidingAtom.angleType[self.getElementType()]
            
            # Create list of angle triplets
            for (offset, bond1) in enumerate(self.bondList[:-1]):
                for bond2 in self.bondList[(offset+1):]:
                    #For Ni, do not use 180 degree angles
                    if (self.getElement() == 'Ni'):
                        bond1Atom = atomList[bond1-1]
                        bond2Atom = atomList[bond2-1]
                        angle = calcAngle(self, bond1Atom, bond2Atom)
                        
                        # Angle is closer to 180 than 90, do not add angle
                        if angle > 135:
                            continue
                        
                    self.angleList.append([n, bond1, self.getID(), bond2])
                    dreidingAtom.angleCount = dreidingAtom.angleCount + 1
   
    def assignCharge(self, atomList):
        # Set atom charge depending on element type and surrounding bonds
        # See notes on 2/6/14
        elementType = self.getElementType()
        
        if elementType[0] == 'O':
            self.charge = -0.517649
        elif elementType[0] == 'B':
            self.charge = 0.754126
        elif elementType == 'C_1':
            # If atom is a C_1 and next to a C_R, assign C1 charge            
            neighborList = self.next2ElementType('C_R', atomList)
            if neighborList[0] == True:
                self.charge = -0.173075
                
            # If atom is a C_1 and next to two C_1 atoms, assign C2 charge            
            neighborList = self.next2ElementType('C_1', atomList)
            if len(neighborList) > 2:
                self.charge = 0.010663
                
        elif elementType == 'C_R':
            # If atom is a C_R and next to a O_3, assign C3 charge            
            neighborList = self.next2ElementType('O_3', atomList)
            if neighborList[0] == True:
                self.charge = 0.292127
                
            # If atom is a C_R and next to a B_2, assign C5 charge            
            neighborList = self.next2ElementType('B_2', atomList)
            if neighborList[0] == True:
                self.charge = -0.294596
                
            # If atom is a C_R and next to a C_1, assign C7 charge            
            neighborList = self.next2ElementType('C_1', atomList)
            if neighborList[0] == True:
                self.charge = 0.300653
                
            # If atom is a C_R and next to three C_R atoms, assign C9 charge            
            neighborList = self.next2ElementType('C_R', atomList)
            if len(neighborList) > 3:
                self.charge = 0.029955
                
            # Look for neighbors of neighbors to determine relative position
            # in comparision to other C_R atoms
            elif self.charge == 0:
                for bondAtom in neighborList[1:]:
                    # If atom is a C_R, next to a C_R, next to O_3, 
                    # assign C4 charge
                    tempNeighborList = bondAtom.next2ElementType('O_3', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = -0.286856
                        
                    # If atom is a C_R, next to a C_R, next to B_2, 
                    # assign C6 charge
                    tempNeighborList = bondAtom.next2ElementType('B_2', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = 0.045103
                        
                    # If atom is a C_R, next to a C_R, next to C_1, 
                    # assign C8 charge
                    tempNeighborList = bondAtom.next2ElementType('B_2', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = -0.222841
                        
        elif elementType == 'H_1':
            # If atom is H_1, use neighbor carbon to determine self charge        
            neighbor = self.next2ElementType('C_R', atomList)
            neighborCarbon = neighbor[1]
            
            # Use neighboring carbon to determine relative position
            neighborList = neighborCarbon.next2ElementType('C_R', atomList)
            for bondAtom in neighborList[1:]:
                    # If atom is a C_R, next to a C_R, next to O_3, 
                    # assign H10 charge
                    tempNeighborList = bondAtom.next2ElementType('O_3', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = 0.174619
                        
                    # If atom is a C_R, next to a C_R, next to B_2, 
                    # assign H11 charge
                    tempNeighborList = bondAtom.next2ElementType('B_2', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = 0.064254
                        
                    # If atom is a C_R, next to a C_R, next to C_1, 
                    # assign H12 charge
                    tempNeighborList = bondAtom.next2ElementType('B_2', atomList)
                    if tempNeighborList[0] == True:
                        self.charge = 0.121582
                        
        elif self.charge == 0:
            warnings.warn("Unknown element type for charge determination")

    def assignManualCharge(self, chargeVal):
        # Set atom charge manually.
        # Originally designed to use external list derived from gaussian
        self.charge = float(chargeVal)

        if self.charge == 0:
            warnings.warn("Unknown element type for charge determination")
        
    # Determines whether atom is next to a certain element type and returns
    # the bondAtom index for the element Type
    def next2ElementType(self, testElementType, atomList):
        # Get list of bonds
        bondList = self.getBondList()
        
        # Initialize list of elment type neighbors with first index being true
        # or false and the following indices being the neighbor IDs
        elementTypeNeighborList = [False]
        
        # Cycle through all bonds
        for bondID in bondList:
            # Get bond atom object and element type
            bondAtom = atomList[bondID-1]
            bondElementType = bondAtom.getElementType()
            
            # Check if bond atom is of the relevant element type. Add to list
            if bondElementType == testElementType:
                elementTypeNeighborList.append(bondAtom)
        
        # Return results, True if any neighbors found, false otherwise
        if len(elementTypeNeighborList) > 0:
            elementTypeNeighborList[0] = True
        return elementTypeNeighborList
   
    def getID(self):
        return self.atomID
        
    def getBondList(self):
        # Return list of bonds [Note: these should already be presorted]
        return self.bondList
        
    def getBondOrderList(self):
        # Return bond orders in the same sequence as bondList
        return self.bondOrderList
        
    def getBondTypeList(self):
        # Return bond types in the same sequence as bondList
        return self.bondTypeList
        
    def getMaxBond(self):
        # Return maximum number of bonds allowed
        return self.maxBond
        
    def getAngleList(self):
        # Return list of angles
        return self.angleList        
    
    def getCount(self):
        return dreidingAtom.atomCount
        
    def getElement(self):
        return self.element
        
    def getElementType(self):
        return self.elementType
              
    def getX(self):
        return self.x
        
    def getY(self):
        return self.y
        
    def getZ(self):
        return self.z
        
    def getRadius(self):
        return self.radius
        
    def getCharge(self):
        return self.charge
        
    def displayAtom(self):
        print('Element: ', self.element,'; (x,y,z): ', self.x, ', ', self.y, \
              ', ', self.z, '; ID: ', self.atomID, sep='')
              
    def __eq__(self, other):
        return calcDistance(self, other) == 0
              

# Return list of elements sorted by element
def getElementList():
    return dreidingAtom.elementList
    
# Return list of bond types
def getGlobalBondTypeList():
    return dreidingAtom.globalBondType
    
# Return list of angle types
def getAngleTypeList():
    return dreidingAtom.angleType
    
# Establish layered box boundary conditions for this data file
# Primarily used for repeating COF system
def setLayerBox(periodic, a, b, angle, n, space):
    # Set whether periodic boundary conditions are used
    dreidingAtom.periodic = periodic
    
    # For lammmps, trigonal vectors established by using xy xz yz
    # A = (xhi-xlo,0,0); 
    # B = (xy,yhi-ylo,0); 
    # C = (xz,yz,zhi-zlo)
    dreidingAtom.xlo = 0
    dreidingAtom.xhi = dreidingAtom.xlo + a
    dreidingAtom.ylo = 0
    dreidingAtom.yhi = dreidingAtom.ylo + b*math.sin(math.radians(angle))
    dreidingAtom.zlo = -5
    dreidingAtom.zhi = 5
    dreidingAtom.xy = b*math.cos(math.radians(angle))
    dreidingAtom.xz = 0
    dreidingAtom.yz = 0
    
    # Check if this is a many layer system
    if n != -1:
        dreidingAtom.zlo = 0
        dreidingAtom.zhi = n * space
        
# Establish triclinic box boundary conditions for this data file
# Currently designed for C60 on pentacene system
def setTriclinicBox(periodic, a, b, c, alpha, beta, gamma):
    # Set whether periodic boundary conditions are used
    dreidingAtom.periodic = periodic
    
    # For lammmps, trigonal vectors established by using xy xz yz
    # A = (xhi-xlo,0,0); 
    # B = (xy,yhi-ylo,0); 
    # C = (xz,yz,zhi-zlo)
    dreidingAtom.xlo = 0
    dreidingAtom.ylo = 0
    dreidingAtom.zlo = 0
    
    # Formula for converting (a,b,c,alpha,beta,gamma) to (lx,ly,lz,xy,xz,yz)
    # taken from online lammps help
    dreidingAtom.xhi = a
    dreidingAtom.xy = b*math.cos(math.radians(gamma))
    dreidingAtom.xz = c*math.cos(math.radians(beta))
    dreidingAtom.yhi = math.sqrt(b**2 - dreidingAtom.xy**2)
    dreidingAtom.yz = (b*c*math.cos(math.radians(alpha)) - dreidingAtom.xy * dreidingAtom.xz)/ dreidingAtom.yhi
    dreidingAtom.zhi = math.sqrt(c**2 - dreidingAtom.xz**2 - dreidingAtom.yz**2)
  
# Calculate distance between atoms
def calcDistance(atom1, atom2):
    # Create dummy box for all calculated distances    
    distance = []    
    
    # Import xyz coordinates for both atoms
    x1 = atom1.getX()
    y1 = atom1.getY()
    z1 = atom1.getZ()
    pos1 = {'x':x1, 'y':y1, 'z':z1}
    
    x2 = atom2.getX()
    y2 = atom2.getY()
    z2 = atom2.getZ()
    pos2 = {'x':x2, 'y':y2, 'z':z2}
    
    # Calculate distance within one box
    distance.append(math.sqrt((pos1['x']-pos2['x'])**2 + (pos1['y']-pos2['y'])**2 + (pos1['z']-pos2['z'])**2))
    
    # If periodic conditions have been established, try calculating distance 
    # between atoms that have been repeated across box boundaries
    # Need to translate 26 times for a 3D space
    if dreidingAtom.periodic:
        # Iterate over -1, 0, and 1 translation for vectors 1/A, 2/B, and 3/C
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    # Translate as many times as required
                    testPos = translateVector1A(pos1,i)
                    testPos = translateVector2B(testPos,j)
                    testPos = translateVector3C(testPos,k)
                    
                    # Calculate new distance
                    distance.append(math.sqrt((testPos['x']-pos2['x'])**2 + (testPos['y']-pos2['y'])**2 + (testPos['z']-pos2['z'])**2))

    # Keep only smallest distance
    distance = min(distance) 
    return distance

# Calculate non-periodic distance between atoms
# Useful when checking if a bond crosses the periodic boundary
# Equivalent to calcDistance() when dreidingAtom.periodic = False
def calcNonPeriodicDistance(atom1, atom2):
    # Import xyz coordinates for both atoms
    x1 = atom1.getX()
    y1 = atom1.getY()
    z1 = atom1.getZ()
    pos1 = {'x':x1, 'y':y1, 'z':z1}

    x2 = atom2.getX()
    y2 = atom2.getY()
    z2 = atom2.getZ()
    pos2 = {'x':x2, 'y':y2, 'z':z2}

    # Calculate distance within one box
    distance = math.sqrt((pos1['x']-pos2['x'])**2 + (pos1['y']-pos2['y'])**2 + (pos1['z']-pos2['z'])**2)
    return distance
    
# Calculate angle between three atoms. atom1 is center atom.
# Angle is degrees from 0-180
def calcAngle(atom1, atom2, atom3):
    # Create variable for minimum distance  
    minDistance12 = 0  
    minDistance13 = 0  
    
    # Import xyz coordinates for all atoms
    x1 = atom1.getX()
    y1 = atom1.getY()
    z1 = atom1.getZ()
    pos1 = {'x':x1, 'y':y1, 'z':z1}
    
    x2 = atom2.getX()
    y2 = atom2.getY()
    z2 = atom2.getZ()
    pos2 = {'x':x2, 'y':y2, 'z':z2}

    x3 = atom3.getX()
    y3 = atom3.getY()
    z3 = atom3.getZ()
    pos3 = {'x':x3, 'y':y3, 'z':z3}
    
    # Calculate distance between atoms 1-2 and atoms 1-3 within one box
    minDistance12 = math.sqrt((pos1['x']-pos2['x'])**2 + (pos1['y']-pos2['y'])**2 + (pos1['z']-pos2['z'])**2)
    minDistance13 = math.sqrt((pos1['x']-pos3['x'])**2 + (pos1['y']-pos3['y'])**2 + (pos1['z']-pos3['z'])**2)
    
    #Calculate vectors relative to atom1
    v12 = {'x':(pos2['x']-pos1['x']), 'y':(pos2['y']-pos1['y']), 'z':(pos2['z']-pos1['z'])}
    v13 = {'x':(pos3['x']-pos1['x']), 'y':(pos3['y']-pos1['y']), 'z':(pos3['z']-pos1['z'])}    
    
    # If periodic conditions have been established, try calculating distance 
    # between atoms that have been repeated across box boundaries
    # Need to translate 26 times for a 3D space
    if dreidingAtom.periodic:
        # Iterate over -1, 0, and 1 translation for vectors 1/A, 2/B, and 3/C
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    # Translate as many times as required for both atom2 
                    # and atom3
                    testPos2 = translateVector1A(pos2,i)
                    testPos2 = translateVector2B(testPos2,j)
                    testPos2 = translateVector3C(testPos2,k)
                    
                    testPos3 = translateVector1A(pos3,i)
                    testPos3 = translateVector2B(testPos3,j)
                    testPos3 = translateVector3C(testPos3,k)
                    
                    # Calculate new distances
                    testDistance12 = math.sqrt((pos1['x']-testPos2['x'])**2 + (pos1['y']-testPos2['y'])**2 + (pos1['z']-testPos2['z'])**2)
                    testDistance13 = math.sqrt((pos1['x']-testPos3['x'])**2 + (pos1['y']-testPos3['y'])**2 + (pos1['z']-testPos3['z'])**2)
                    
                    # Check if either distance is closer than previously 
                    # calculated distance. If so, change the vector v12 or v13
                    if testDistance12 < minDistance12:
                        v12 = {'x':(testPos2['x']-pos1['x']), 'y':(testPos2['y']-pos1['y']), 'z':(testPos2['z']-pos1['z'])}
                        minDistance12 = testDistance12
                        
                    if testDistance13 < minDistance13:
                        v13 = {'x':(testPos3['x']-pos1['x']), 'y':(testPos3['y']-pos1['y']), 'z':(testPos3['z']-pos1['z'])}
                        minDistance13 = testDistance13


    # Calculate angle using vectors 1 and 2. Vectors should be the shortest
    # possible
    # v*w = |v||w|cos(theta)
    try:
        cos_angle = ((v12['x']*v13['x'] + v12['y']*v13['y'] + v12['z']*v13['z']) /
                    (minDistance12 * minDistance13))
        angle = math.degrees(math.acos(cos_angle))
    except:
        print('Error in calcAngle: math error in calculating angle')
        angle = 0
    
    # Normalize angle to be between 0-180 degrees
    if angle > 180:
        angle = 360 - angle
        
    return angle

# Translate x,y,z coordinates across vector 1/A, as many times as 'multiple'
def translateVector1A(pos, multiple):
    x_a = pos['x'] + (dreidingAtom.xhi - dreidingAtom.xlo) * multiple
    y_a = pos['y'] + 0
    z_a = pos['z'] + 0
    return {'x':x_a, 'y':y_a, 'z':z_a}

# Translate x,y,z coordinates across vector 2/B, as many times as 'multiple'
def translateVector2B(pos, multiple):
    x_b = pos['x'] + dreidingAtom.xy * multiple
    y_b = pos['y'] + (dreidingAtom.yhi - dreidingAtom.ylo) * multiple
    z_b = pos['z'] + 0
    return {'x':x_b, 'y':y_b, 'z':z_b}
    
# Translate x,y,z coordinates across vector 3/C, as many times as 'multiple'
def translateVector3C(pos, multiple):
    x_c = pos['x'] + dreidingAtom.xz * multiple
    y_c = pos['y'] + dreidingAtom.yz * multiple
    z_c = pos['z'] + (dreidingAtom.zhi - dreidingAtom.zlo) * multiple
    return {'x':x_c, 'y':y_c, 'z':z_c}
    
def calcAVGBondLength(atom1, atom2):
    radius1 = atom1.getRadius()
    radius2 = atom2.getRadius()
    avgBondLength = radius1 + radius2
    return avgBondLength
    
# Creates a bond between to atoms by adding to their bond list
def createBond(atom1, atom2):
    if atom1.assignBond(atom2) and atom2.assignBond(atom1):
        dreidingAtom.bondCount = dreidingAtom.bondCount + 1
        #print(dreidingAtom.bondCount, atom1.getID(), atom2.getID())
        
# Change bond order between two atoms
def changeBondOrder(atom1, atom2, newBondOrder):
    # Get bond and ID properties from atom 1 and 2
    bondList1 = atom1.getBondList()
    atom1_ID = atom1.getID()
    bondList2 = atom2.getBondList()
    atom2_ID = atom2.getID()
    
    # Change bond orders within both atoms to match
    index1 = bondList1.index(atom2_ID)
    atom1.bondOrderList[index1] = newBondOrder
    index2 = bondList2.index(atom1_ID)
    atom2.bondOrderList[index2] = newBondOrder
    
# Change bond type between two atoms
def changeBondType(atom1, atom2, newBondType):
    # Get bond and ID properties from atom 1 and 2
    bondList1 = atom1.getBondList()
    atom1_ID = atom1.getID()
    bondList2 = atom2.getBondList()
    atom2_ID = atom2.getID()
    
    # Change bond types within both atoms to match
    index1 = bondList1.index(atom2_ID)
    atom1.bondTypeList[index1] = newBondType
    index2 = bondList2.index(atom1_ID)
    atom2.bondTypeList[index2] = newBondType
    
# Finds the closest atom(s) to atom1. numberCloseAtoms tells how many closest 
# atoms should be found
def findClosest(atom1, atomList, numberCloseAtoms):
    # Create a list for closest atom(s) and distance(s)
    maxDist = 2.4
    closeAtomList = [[None, maxDist]] * numberCloseAtoms

    # Compare atom1 to every atom
    for testAtom in atomList:
        # Calculate atom distance
        testDist = calcDistance(atom1, testAtom)

        # Add testAtom to closest atom list if it is closer than at least one
        # of the atoms on the list
        for (offset, closeAtom) in enumerate(closeAtomList):
            # If atom is close, move the list up from that point and remove
            # furthest atom
            if testDist <= closeAtom[1] and testDist != 0:
                closeAtomList.insert(offset,[testAtom, testDist])
                closeAtomList.pop()
                break

    # Get rid of extra slots and display warning
    returnAtomList = [closeAtom for closeAtom in closeAtomList if maxDist != closeAtom[1]]
    if len(returnAtomList) != numberCloseAtoms:
            print("Warning: Less than requested # of closest atoms. maxDist: " + str(maxDist))

    return returnAtomList

# Determine number of periodic bonds and prints it out.
# Requires element type information
def determinePeriodicBonds(atomList):
    # Periodic box parameters
    #half_x = (dreidingAtom.xhi - dreidingAtom.xlo) / 2
    #half_y = (dreidingAtom.yhi - dreidingAtom.ylo) / 2
    #half_z = (dreidingAtom.zhi - dreidingAtom.zlo) / 2

    # Turn off any of the periodic directions if it is too small
    #if half_x < 5:
    #    half_x = 5
    #if half_y < 5:
    #    half_y = 5
    #if half_z < 5:
    #    half_z = 5

    # Store all periodic bonds in arrays
    downAtomIDSet = []
    upAtomIDSet = []
    rightAtomIDSet = []
    leftAtomIDSet = []

    # Cycle through every atom
    for atom in atomList:
        # Get ID and all the bonds associated with the atom
        currentID = atom.getID()
        bondList = atom.getBondList()

        # Cycle through every bond
        for bondAtomID in bondList:
            if bondAtomID > currentID:
                bondAtom = atomList[bondAtomID-1]

                # Calculate non-periodic distance and compare to periodic dist
                nonPeriodicDist = calcNonPeriodicDistance(atom, bondAtom)
                periodicDist = calcDistance(atom, bondAtom)
                #print(nonPeriodicDist)

                #if nonPeriodicDist > half_x or nonPeriodicDist > half_y or nonPeriodicDist > half_z:
                if nonPeriodicDist > periodicDist + 1:

                    #atom.displayAtom()
                    #bondAtom.displayAtom()

                    # Determine location of 'atom'
                    # Check if down atom
                    if abs(dreidingAtom.ylo - atom.getY()) < 5:
                        downAtomIDSet.append(currentID)
                        print('downAtomID: %d' % (currentID))
                    # Check if up atom
                    elif abs(dreidingAtom.yhi - atom.getY()) < 5:
                        upAtomIDSet.append(currentID)
                        print('upAtomID: %d' % (currentID))
                    # Check if right atom
                    elif abs(dreidingAtom.xhi - atom.getX()) < 5:
                        rightAtomIDSet.append(currentID)
                        print('rightAtomID: %d' % (currentID))
                    # Check if left atom
                    elif abs(dreidingAtom.xlo - atom.getX()) < 5:
                        leftAtomIDSet.append(currentID)
                        print('leftAtomID: %d' % (currentID))

                    # Determine location of 'bondAtom'
                    # Check if down atom
                    if abs(dreidingAtom.ylo - bondAtom.getY()) < 5:
                        downAtomIDSet.append(bondAtomID)
                        print('downAtomID: %d' % (bondAtomID))
                    # Check if up atom
                    elif abs(dreidingAtom.yhi - bondAtom.getY()) < 5:
                        upAtomIDSet.append(bondAtomID)
                        print('upAtomID: %d' % (bondAtomID))
                    # Check if right atom
                    elif abs(dreidingAtom.xhi - bondAtom.getX()) < 5:
                        rightAtomIDSet.append(bondAtomID)
                        print('rightAtomID: %d' % (bondAtomID))
                    # Check if left atom
                    elif abs(dreidingAtom.xlo - bondAtom.getX()) < 5:
                        leftAtomIDSet.append(bondAtomID)
                        print('leftAtomID: %d' % (bondAtomID))

                    #print('%d %d' % (currentID, bondAtomID))

    # Output periodic bonds to external file
    txt = ''
    txt = txt + 'down ' + ' '.join(map(str, downAtomIDSet)) + '\n'
    txt = txt + 'up ' + ' '.join(map(str, upAtomIDSet)) + '\n'
    txt = txt + 'right ' + ' '.join(map(str, rightAtomIDSet)) + '\n'
    txt = txt + 'left ' + ' '.join(map(str, leftAtomIDSet))
    return txt
                
# Determine bond order for all bonds.
# Requires element type information
def determineBondOrders(atomList):
    # Cycle through every atom
    for atom in atomList:
        # Get ID and all the bonds associated with the atom
        currentID = atom.getID()
        bondList = atom.getBondList()
        
        # Cycle through every bond
        for bondAtomID in bondList:
            if bondAtomID > currentID:
                bondAtom = atomList[bondAtomID-1]
                
                # Determine bond order from element type
                bondAtom_type = bondAtom.getElementType()
                atom_type = atom.getElementType()
                
                # If there are two resonance carbons, change to bond order 1.5
                if bondAtom_type == 'C_R' and atom_type == 'C_R':
                    changeBondOrder(bondAtom, atom, 1.5)
                    
# Determine bond type for all bonds. Compile global list
# Requires element type information
def determineBondTypes(atomList):
    # Cycle through every atom
    for atom in atomList:
        # Get ID and all the bonds associated with the atom
        currentID = atom.getID()
        bondList = atom.getBondList()
        
        # Cycle through every bond
        for bondAtomID in bondList:
            if bondAtomID > currentID:
                bondAtom = atomList[bondAtomID-1]
                
                # Determine bond type from element type pair
                bondAtom_type = bondAtom.getElementType()
                atom_type = atom.getElementType()
                
                bond_type = [atom_type, bondAtom_type]
                bond_type.sort()
                
                # Create running list of bond pair types
                for (offset, pair) in enumerate(dreidingAtom.globalBondType):
                    if pair[0] == bond_type[0] and pair[1] == bond_type[1]:
                        changeBondType(bondAtom, atom, offset+1)
                        break
                else:
                    dreidingAtom.globalBondType.append(bond_type)
                    changeBondType(bondAtom, atom, len(dreidingAtom.globalBondType))
                    
# Determine charges for all bonds
# Requires CHELPG charges from Gaussian
# Requires element type information to be accurate
# Requires bonds to be accurate
def determineCharges(atomList):
    # Cycle through every atom
    for atom in atomList:
        # Get ID and all the bonds associated with the atom
        atom.assignCharge(atomList)
        
# Count all possible dihedrals. Must run before outputting system info
def determineDihedralCount(atomList):
    # Cycle through every atom
    # Dihedral nomenclature is ijkl where j and k are middle two atoms
    for j_atom in atomList:
        # Get ID and all the bonds associated with the atom
        j_ID = j_atom.getID()
        j_bondList = j_atom.getBondList()
        
        # Cycle through every atom bonded to j to find potential k atoms
        for (offset, k_ID) in enumerate(j_bondList):
            # Forces jk pair to be sequentially unique
            if k_ID > j_ID:
                # Cycle through every atom bonded to j (again), in order to 
                # find all possible i atoms
                for i_ID in j_bondList:
                    # Check to make sure i_atom is not already k_atom
                    if i_ID != k_ID:
                        #Find list of atoms bonded to k atom
                        k_atom = atomList[k_ID-1]
                        k_bondList = k_atom.getBondList()
                        # Cycle through every atom bonded to k, in order to
                        # find all possible l atoms
                        for l_ID in k_bondList:
                            # Check to make sure l_ID is not already j_ID
                            if l_ID != j_ID:
                                # Increment dihedral count                            
                                dreidingAtom.dihedralCount = dreidingAtom.dihedralCount + 1
                                
# Count all possible impropers. Must run before outputting system info
def determineImproperCount(atomList):
    # Cycle through every atom
    # Improper nomenclature is ijkl where i is middle two atoms
    for i_atom in atomList:
        # Get atom element type.      
        i_type = i_atom.getElementType()
        
        # Get atom id and bondlist
        i_ID = i_atom.getID()
        i_bondList = i_atom.getBondList()
        
        # If this is an atom bonded to three other atoms and it forms an improper,
        # add to list
        if (i_type[2] == 'R' or i_type[2] == '2') and len(i_bondList) == 3:
            # Increment improper number
            dreidingAtom.improperCount = dreidingAtom.improperCount + 1

# Check if atom is connected in a ring by recursively checking if it finds itself
# Will only return first ring it finds
# Input: How long the ring should be
# Output: (True or False), list of ring atoms
def isRingAtom(originAtom, atom, count):
    # Get ID and all the bonds associated with the atom
    currentID = atom.getID()
    bondList = atom.getBondList()
    
    # Cycle through every bond
    for bondAtomID in bondList:
        if bondAtomID > currentID:
            bondAtom = atomList[bondAtomID-1]
            
            # Determine bond type from element type pair
            bondAtom_type = bondAtom.getElementType()
            atom_type = atom.getElementType()
    return {'x':x_a, 'y':y_a, 'z':z_a}
                                
# Prepare buckyball parameters with appropriate bond lengths and angles
# Add appropriate bond coeff and bond angle coeff
def fixBuckyParams(atomList, buckyList):
    # Add new bond type if it hasn't already been added ['C_3', 'C_3']
    bond_type = ['C_3', 'C_3']
    # Create running list of bond pair types
    for (offset, pair) in enumerate(dreidingAtom.globalBondType):
        if pair[0] == bond_type[0] and pair[1] == bond_type[1]:
            break
    else:
        dreidingAtom.globalBondType.append(bond_type)
        
    # Add new angle type if it hasn't already been added C_3
    elementType = 'C_3'
    if (elementType in dreidingAtom.angleType) == False:
        # If angle has not been added, add to the dictionary with new 
        # number
        n = len(dreidingAtom.angleType)
        dreidingAtom.angleType[elementType] = n+1
        
    # Change bond type between two atoms
    changeBondType(atom1, atom2, newBondType)
    
    # Change bond types between pentagonal bridges in buckyball
    # Cycle through every atom
    for atom in buckyList:
        # Get ID and all the bonds associated with the atom
        currentID = atom.getID()
        bondList = atom.getBondList()
        
        # Cycle through every bond
        for bondAtomID in bondList:
            if bondAtomID > currentID:
                bondAtom = atomList[bondAtomID-1]
                
                # Determine bond type from element type pair
                bondAtom_type = bondAtom.getElementType()
                atom_type = atom.getElementType()
                
                bond_type = [atom_type, bondAtom_type]
                bond_type.sort()
                
                # Create running list of bond pair types
                for (offset, pair) in enumerate(dreidingAtom.globalBondType):
                    if pair[0] == bond_type[0] and pair[1] == bond_type[1]:
                        changeBondType(bondAtom, atom, offset+1)
                        break
                else:
                    dreidingAtom.globalBondType.append(bond_type)
                    changeBondType(bondAtom, atom, len(dreidingAtom.globalBondType))
        
# Print system info as required by the intro to lammps data file
def getSystemInfo():
    # Compile number of items
    #Temporarily turn off all bonds
    #dreidingAtom.bondCount = 0
    #dreidingAtom.angleCount = 0
    #dreidingAtom.dihedralCount = 0
    text = "LAMMPS Description" + "\n\n"
    text = text + str(dreidingAtom.atomCount) + '  atoms\n'
    text = text + str(dreidingAtom.bondCount) + '  bonds\n'
    text = text + str(dreidingAtom.angleCount) + '  angles\n'
    text = text + str(dreidingAtom.dihedralCount) + '  dihedrals\n'
    text = text + str(dreidingAtom.improperCount) + '  impropers\n\n'
    
    # Compile number of each type of item
    atomTypeCount = len(getElementList())
    text = text + str(atomTypeCount) + '  atom types\n'
    bondTypeCount = len(getGlobalBondTypeList())
    text = text + str(bondTypeCount) + '  bond types\n'
    angleTypeCount = len(getAngleTypeList())
    text = text + str(angleTypeCount) + '  angle types\n'
    dihedralTypeCount = 10
    text = text + str(dihedralTypeCount) + '  dihedral types\n'
    improperTypeCount = 1
    text = text + str(improperTypeCount) + '  improper types\n\n'
    
    # Compile box information
    temp = '%.6f %.6f xlo xhi\n' % (dreidingAtom.xlo, dreidingAtom.xhi)
    text = text + temp
    temp = '%.6f %.6f ylo yhi\n' % (dreidingAtom.ylo, dreidingAtom.yhi)
    text = text + temp
    temp = '%.6f %.6f zlo zhi\n' % (dreidingAtom.zlo, dreidingAtom.zhi)
    text = text + temp
    temp = '%.6f %.6f %.6f xy xz yz\n' % (dreidingAtom.xy, dreidingAtom.xz, dreidingAtom.yz)
    text = text + temp
    
    # Return compiled text
    text = text + '\n'
    return text
    
# Print mass as required by the lammps data file
def getMassInfo():
    # Compile list of masses
    text = "Masses\n\n"
    
    # Use list of fundamental elements to make mass list    
    elementList = getElementList()
    
    # Organize atoms into lists of same elements
    # Assign mass (g/mol)
    for (offset, group) in enumerate(elementList):
        if group[0] == 'H':
            mass = 1.00794
        elif group[0] == 'B':
            mass = 10.811
        elif group[0] == 'C':
            mass = 12.0107
        elif group[0] == 'N':
            mass = 14.0067
        elif group[0] == 'O':
            mass = 15.9994
        elif group[0] == 'Ni':
            mass = 58.6934
        elif group[0] == 'CB':
            mass = 720.64
        else:
            # Will hopefully be unmovable with high mass
            warnings.warn("Unknown element type")
            mass = 100000.0
            
        temp = ' %d %.4f # %s\n' % (offset+1, mass, group[0])
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Print pair coeffs as required by the lammps data file
def getPairCoeffsInfo():
    # Compile list of pair coefficients
    text = "Pair Coeffs\n# pair_style lj/cut/coul/long nm/cut\n"
    
    # Use list of fundamental elements to make pair list    
    elementList = getElementList()
    
    # Organize atoms into lists of same elements
    # R0 = Angstroms
    # D0 = kcal/mol
    for (offset, group) in enumerate(elementList):

        # Assume pair_style hybrid
        # Default pair_style is lj/cut/coul/long
        type = 'lj/cut/coul/long'

        if group[0] == 'H':
            R0 = 3.195
            D0 = 0.0152
        elif group[0] == 'B':
            R0 = 4.02
            D0 = 0.095
        elif group[0] == 'C':
            R0 = 3.8983
            D0 = 0.0951
        elif group[0] == 'N':
            R0 = 3.6621
            D0 = 0.0774
        elif group[0] == 'O':
            R0 = 3.4046
            D0 = 0.0957
        # Use published Ni  energy (Shelnutt et al. 1991)
        elif group[0] == 'Ni':
            R0 = 2.270
            D0 = 0.055
        # C60 uses pair_style nm/cut
        # Generate custom text and move on
        elif group[0] == 'CB':
            type = 'nm/cut'
            R0 = 10.048
            D0 = 6.648

            C60_pair = '6.6480 10.0480 35.4877 8.8719'

            temp = ' %d %s %s # %s\n' % (offset+1, type, C60_pair, group[0])
            text = text + temp
            continue

        else:
            warnings.warn("Unknown element type")
            R0 = 1.0
            D0 = 0.0
            
        # Convert to lammps units: eps and sigma
        eps = D0
        sigma = R0/(2**(1/6))
            
        temp = ' %d %s %.4f %.4f # %s\n' % (offset+1, type, eps, sigma, group[0])
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Determine bond radius from element type
def getBondRadius(elementType):
    if elementType == 'H_1':
        rad = 0.33
    elif elementType == 'B_2':
        rad = 0.79
    elif elementType == 'C_3':
        rad = 0.77
    elif elementType == 'C_R':
        rad = 0.7
    elif elementType == 'C_1':
        rad = 0.602
    elif elementType == 'O_3':
        rad = 0.66
    elif elementType == 'N_3':
        rad = 0.702
    elif elementType == 'N_R':
        rad = 0.65
    elif elementType == 'N_2':
        rad = 0.615
    elif elementType == 'Ni3':
        rad = 1.215
    else:
        # Will not bond if given a huge radius unless it is specially bonded
        # using findClosest and manually bonding
        warntext = "Unknown element type for bond radius: " + elementType
        warnings.warn(warntext)
        rad = 1000
        
    return rad
    
# Determine bond energy from two element type in kcal/(mol * Angstrom^2)
def getBondEnergy(elementType1, elementType2):
    # Add bond orders to a list.  Resonance atoms get bond order 1.5
    hybrid = []   
    # Add atom 1
    if elementType1[2] == 'R':
        hybrid.append(1.5)
    else:
        hybrid.append(int(elementType1[2]))
    # Add atom 2 
    if elementType2[2] == 'R':
        hybrid.append(1.5)
    else:
        hybrid.append(int(elementType2[2]))
        
    # Check for boron and resonance exceptions
    if elementType1 == 'B_2' and elementType2 =='B_2':
        energy = 1400
    elif elementType1 == 'B_2' and elementType2 == 'C_R':
        energy = 1400
    elif elementType1 == 'C_2' and elementType2 == 'C_R':
        energy = 1400
    # Use published Ni-N energy (Shelnutt et al. 1991)
    elif elementType1 == 'Ni3' or elementType2 == 'Ni3':
        energy = 242
    else:
        # Assign energy as a function of lowest hybridization
        minHybrid = min(hybrid)
        #bondOrder = 4 - minHybrid
        bondOrder = minHybrid
        energy = 700 * bondOrder
        
    return energy
    
# Print bond coefficients as required by the lammps data file
def getBondCoeffsInfo():
    # Add header. Bond_style info not used by input file
    text = "Bond Coeffs\n# bond_style harmonic\n"
    
    # Use list of fundamental bond pair types to make bond coeff list    
    globalBondType = getGlobalBondTypeList()
    #print(globalBondType)
    
    # Cycle through list of bond pair types
    for (offset, pair) in enumerate(globalBondType):
        rad1 = getBondRadius(pair[0])
        rad2 = getBondRadius(pair[1])
        equilBond = rad1 + rad2 - 0.01
        
        energy = getBondEnergy(pair[0], pair[1])
        LAMMPSenergy = .5 * energy
            
        temp = ' %d %.4f %.6f # %s\n' % (offset+1, LAMMPSenergy, equilBond, pair)
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Print angle coefficients as required by the lammps data file
def getAngleCoeffsInfo():
    # Add header. Angle_style info not used by input file
    text = "Angle Coeffs\n# angle_style harmonic\n"
    
    # Use list of fundamental bond pair types to make angle coeff list    
    globalAngleType = getAngleTypeList()
    temp = [None] * len(globalAngleType)
    
    # Cycle through list of bond pair types
    for angleType in globalAngleType:
        # Assign energy (energy/rad^2)
        # Independent of element
        K = 100 # (kcal/mol)/rad^2
        K = K / 2 # Formatted for lammps
        
        # Assign optimal angle in degrees
        # Lammps automatically converts to radians
        if angleType == 'C_1':
            degrees = 180
        elif angleType == 'C_R' or angleType == 'C_2':
            degrees = 120
        elif angleType == 'C_3':
            degrees = 109.471
        elif angleType == 'B_2':
            degrees = 120
        elif angleType == 'O_3':
            degrees = 104.51
        elif angleType == 'N_R':
            degrees = 120
        elif angleType == 'N_2':
            degrees = 120
        elif angleType == 'N_1':
            degrees = 180
        # Use published Ni angle energy (Shelnutt et al. 1991)
        elif angleType == 'Ni3':
            degrees = 90            
            K = 36
        else:
            warntext = "Unknown angle type: " + angleType
            warnings.warn(warntext)
            degrees = 180
            
        # Print bond coeff text in appropriate row of a new list according to 
        # angle number
        index = globalAngleType[angleType]
        
        temp[index-1] = ' %d %f %.3f # %s' % (index, K, degrees, angleType)
            
    # Compile text list into one unified output
    text = text + '\n'.join(temp)
    
    # Return compiled text
    text = text + '\n\n'
    return text
    
# Print dihedral coefficients as required by the lammps data file
def getDihedralCoeffsInfo():
    # Add header. Dihedral_style info not used by input file
    text = "Dihedral Coeffs\n# dihedral_style charmm\n"
    
    # Since there will always only be 10 types, print all dihedral types in the
    # order established by Mayo et al. 2002
    
    # Formula in Mayo et al. 2002
    # E = 1/2 * V * (1 - cos(n(deg - d)))
    
    # Formula used in dihedral_style charmm
    # E = K * (1 + cos(n*deg - d))
    # K (energy)
    # n (integer >= 0)
    # d (integer value of degrees)
    
    # Cycle through all 10 dihedral types and print dihedral information
    # Type a
    V = 2
    n = 3
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (1, K, n, d_charmm, 0.0, 'Type A')
    
    # Type b
    V = 1
    n = 6
    d = 0
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (2, K, n, d_charmm, 0.0, 'Type B')
    
    # Type c
    V = 45
    n = 2
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (3, K, n, d_charmm, 0.0, 'Type C')
    
    # Type d
    V = 25
    n = 2
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (4, K, n, d_charmm, 0.0, 'Type D')
    
    # Type e
    V = 5
    n = 2
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (5, K, n, d_charmm, 0.0, 'Type E')
    
    # Type f
    V = 10
    n = 2
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (6, K, n, d_charmm, 0.0, 'Type F')
    
    # Type g
    V = 0
    n = 1
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (7, K, n, d_charmm, 0.0, 'Type G')
    
    # Type h
    V = 2
    n = 2
    d = 90
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (8, K, n, d_charmm, 0.0, 'Type H')
    
    # Type i
    V = 2
    n = 2
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (9, K, n, d_charmm, 0.0, 'Type I')
    
    # Type j
    V = 2
    n = 3
    d = 180
    K = (1/2) * V
    d_charmm = n * d + 180
    text = text + ' %d %.1f %d %d %.1f # %s\n' % (10, K, n, d_charmm, 0.0, 'Type J')

    # Return compiled text
    text = text + '\n'
    return text

# Print improper coefficients as required by the lammps data file
def getImproperCoeffsInfo():
    # Add header. improper_style info not used by input file
    text = "Improper Coeffs\n# improper_style umbrella\n"
    # Since there is only one style in dreiding for planar impropers, always include
    text = text + ' %d %.1f %.1f # %s\n' % (1, 40.0, 0.0, 'Planar, non-N column')
    
    # Return compiled text
    text = text + '\n'
    return text
    
# Systematically print all atoms in format required by lammps data file
def getAtomListInfo(atomList):
    # Add header. Atom_style info not used by input file
    text = "Atoms\n# atom_style full\n"
    
    # Get element list
    # Format: [['B', atomID1, atomID2], ['H',atomID1, atomID2]]
    elementList = getElementList()
    
    # Create element dictionary
    elementNumDict = {}
    for (offset, group) in enumerate(elementList):
        elementNumDict[group[0]] = offset + 1

    # Cycle through every atom
    for atom in atomList:
        # Get all atom information
        atomID = atom.getID()
        element = atom.getElement()
        elementType = atom.getElementType()
        elementID = elementNumDict[element]
        #Set C60 as different molecule
        if element == 'CB':
            moleculeID = 2
        else:
            moleculeID = 1
        charge = atom.getCharge()
        #For now, set all charges to 0
        #charge = 0.0
        X = atom.getX()
        Y = atom.getY()
        Z = atom.getZ()
        offsetX = 0
        offsetY = 0
        offsetZ = 0
        
        text = text + ' %d %d %d %.6f' % (atomID, moleculeID, elementID, charge)
        text = text + ' %.5f %.5f %.5f %d %d %d # %s\n' % (X, Y, Z, offsetX, offsetY, offsetZ, elementType)

    # Return compiled text
    text = text + '\n'
    return text
        
# Systematically print all bonds in format required by lammps data file
def getBondListInfo(atomList):
    # Add header.
    text = "Bonds\n\n"
    
    bondNum = 0
    # Cycle through every atom
    for atom in atomList:
        # Get ID and all the bonds associated with the atom
        currentID = atom.getID()
        bondList = atom.getBondList()
        
        # Cycle through every bond
        for (offset, bondAtomID) in enumerate(bondList):
            if bondAtomID > currentID:
                # Increment bond number
                bondNum = bondNum + 1
                
                # Get bondAtom object
                bondAtom = atomList[bondAtomID-1]
                
                # Determine bond type from element type pair
                bondAtom_type = bondAtom.getElementType()
                atom_type = atom.getElementType()
                
                bond_type = [atom_type, bondAtom_type]
                bond_type.sort()
                
                # Find bond coeff number aka bond type number
                bondTypeList = atom.getBondTypeList()
                bondTypeID = bondTypeList[offset]
                 
                text = text + ' %d %d %d %d # %s\n' % (bondNum, bondTypeID, currentID, bondAtomID, bond_type)
                
    # Return compiled text
    text = text + '\n'
    return text
                
# Systematically print all angles in format required by lammps data file
def getAngleListInfo(atomList):
    # Add header.
    text = "Angles\n\n"
    
    angleNum = 0
    # Cycle through every atom
    for atom in atomList:
        # Get all angles and print in order
        angleList = atom.getAngleList()
        
        # Cycle through every angle
        for angle in angleList:
            angleNum = angleNum + 1
            text = text + ' %d %d %d %d %d\n' % (angleNum, angle[0], angle[1], angle[2], angle[3])
            
    # Return compiled text
    text = text + '\n'
    return text

# Systemically print all dihedrals in format required by lammps data file
def getDihedralListInfo(atomList):
    # Add header.
    text = "Dihedrals\n\n "
    
    dreidingAtom.dihedralCount = 0
    # Cycle through every atom
    # Dihedral nomenclature is ijkl where j and k are middle two atoms
    for j_atom in atomList:
        # Get ID and all the bonds associated with the atom
        j_ID = j_atom.getID()
        j_bondList = j_atom.getBondList()
        j_bondOrderList = j_atom.getBondOrderList()
        
        # Cycle through every atom bonded to j to find potential k atoms
        for (offset, k_ID) in enumerate(j_bondList):
            jk_bondOrder = j_bondOrderList[offset]
            # Forces jk pair to be sequentially unique
            if k_ID > j_ID:
                # Cycle through every atom bonded to j (again), in order to 
                # find all possible i atoms
                for i_ID in j_bondList:
                    # Check to make sure i_atom is not already k_atom
                    if i_ID != k_ID:
                        #Find list of atoms bonded to k atom
                        k_atom = atomList[k_ID-1]
                        k_bondList = k_atom.getBondList()
                        # Cycle through every atom bonded to k, in order to
                        # find all possible l atoms
                        for l_ID in k_bondList:
                            # Check to make sure l_ID is not already j_ID
                            if l_ID != j_ID:
                                # Increment dihedral count                            
                                dreidingAtom.dihedralCount = dreidingAtom.dihedralCount + 1
                                
                                # Find element type for all atoms
                                i_atom = atomList[i_ID-1]
                                l_atom = atomList[l_ID-1]
                                i_type = i_atom.getElementType()
                                j_type = j_atom.getElementType()
                                k_type = k_atom.getElementType()
                                l_type = l_atom.getElementType()
                                
                                # Determine type of bond; Currently does not 
                                # check for bond order
                                
                                # Check if j,k atoms are sp3 oxygen atoms
                                if (j_type == 'O_3' and k_type == 'O_3'):
                                    # Type (h)[8] Dreiding dihedral
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 8, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                    
                                # Check if one atom is Ni3 and set to 0 energy
                                elif (j_type == 'Ni3' or k_type == 'Ni3'):
                                    # Type (g)[7] Dreiding dihedral; single bond
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 7, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                
                                # Check if one atom is sp3 oxygen and other is sp2
                                elif ((j_type == 'O_3' and k_type[2] == 'R') or
                                    (j_type == 'O_3' and k_type[2] == '2') or
                                    (j_type[2] == 'R' and k_type == 'O_3') or
                                    (j_type[2] == '2' and k_type == 'O_3')):
                                    # Type (i)[9] Dreiding dihedral
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 9, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                    
                                # Check if both j,k atoms are sp3
                                elif j_type[2] == '3' and k_type[2] == '3':
                                    # Type (a)[1] Dreiding dihedral
                                    # Implicit single bond
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 1, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                    
                                # Check if j,k atoms are sp2 and sp3
                                elif ((j_type[2] == '2' and k_type[2] == '3') or
                                      (j_type[2] == 'R' and k_type[2] == '3') or
                                      (j_type[2] == '3' and k_type[2] == '2') or
                                      (j_type[2] == '3' and k_type[2] == 'R')):
                                    # Check for the sp2-sp2-sp3-x exception
                                    if (((i_type[2] == '2' or i_type[2] == 'R') and k_type[2] == '3') or
                                        ((l_type[2] == '2' or l_type[2] == 'R') and j_type[2] == '3')):
                                        # Type (j)[10] Dreiding dihedral
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 10, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                    else:
                                        # Type (b)[2] Dreiding dihedral
                                        # Implicit single bond
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 2, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                        
                                # Check if j,k atoms are both sp2
                                elif (j_type[2] == '2' and k_type[2] == '2'):
                                    if jk_bondOrder == 2:
                                        # Type (c)[3] Dreiding dihedral; double bond
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 3, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                    elif jk_bondOrder == 1:
                                        # Type (e)[5] Dreiding dihedral; single bond
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 5, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                
                                # Check if j,k atoms are both resonance
                                elif (j_type[2] == 'R' and k_type[2] == 'R'):
                                    if jk_bondOrder == 1.5:
                                        # Type (d)[4] Dreiding dihedral;bond order 1.5
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 4, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                    elif jk_bondOrder == 1:
                                        # Type (e)[5] Dreiding dihedral; single bond
                                        dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 5, i_ID, j_ID, k_ID, l_ID]))
                                        dreidingAtom.dihedralList.append(dihedral)
                                
                                # Check if j,k atoms are a combination of sp2 and R
                                elif ((j_type[2] == '2' and k_type[2] == 'R') or
                                    (j_type[2] == 'R' and k_type[2] == '2')):
                                    # Type (e)[5] Dreiding dihedral; single bond
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 5, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                
                                # Exception exists for 2 aromatic rings
                                elif (j_type[2] == 'R' and k_type[2] == 'R'):
                                    # Type (f)[6] Dreiding dihedral; single bond
                                    # **Ignored for now
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 6, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                
                                # Check if j or k is an sp atom
                                elif (j_type[2] == '1' or k_type[2] == '1'):
                                    # Type (g)[7] Dreiding dihedral; bond order 1 or 3
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 7, i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
                                    
                                # Return warning if no dihedral found
                                else:
                                    dihedral = " ".join(map(str, [dreidingAtom.dihedralCount, 'Error', i_ID, j_ID, k_ID, l_ID]))
                                    dreidingAtom.dihedralList.append(dihedral)
    
    # Compile text list into one unified output
    text = text + '\n '.join(dreidingAtom.dihedralList)    

    # Return compiled text
    text = text + '\n\n'
    return text
    
# Systematically print all impropers in format required by lammps data file
def getImproperListInfo(atomList):
    # Add header.
    text = "Impropers\n\n"
    
    improperNum = 0
    # Cycle through every atom
    for i_atom in atomList:
        # Get atom element type.      
        i_type = i_atom.getElementType()
        
        # Get atom id and bondlist
        i_ID = i_atom.getID()
        i_bondList = i_atom.getBondList()
        
        # If this is an atom bonded to three other atoms and it forms an improper,
        # add to list
        if (i_type[2] == 'R' or i_type[2] == '2') and len(i_bondList) == 3:
            # Increment improper number
            improperNum = improperNum + 1
            
            # Find id's to other three atoms from bondList
            j_ID = i_bondList[0]
            k_ID = i_bondList[1]
            l_ID = i_bondList[2]
            
            # Get atom references
            j_atom = atomList[j_ID-1]
            k_atom = atomList[k_ID-1]
            l_atom = atomList[l_ID-1]
            
            # Get atom type for reference in final data file
            j_type = j_atom.getElementType()
            k_type = k_atom.getElementType()
            l_type = l_atom.getElementType()
            improper_type = [i_type, j_type, k_type, l_type]
            
            # There is currently only one improper type for dreiding model
            improperTypeID = 1
            
            text = text + ' %d %d %d %d %d %d # %s\n' % (improperNum, improperTypeID, i_ID, j_ID, k_ID, l_ID, improper_type)

    # Return compiled text
    text = text + '\n'
    return text
