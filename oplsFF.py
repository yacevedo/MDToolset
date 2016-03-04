"""
Author: Yaset Acevedo
Created: 4/2/14
Update Date: 4/14/14
Version: 1.3
"""

#Import modern print module
from __future__ import print_function
from __future__ import division
import warnings
import math

class oplsAtom:
    'Common base class for all atoms'
    # Initialize global class atom variables
    atomCount = 0
    buckyStart = 10^7
    buckyCount = 0
    buckyCoarse = False
    # Format: [['B', atomID1, atomID2], ['H',atomID1, atomID2]]
    elementList = []
    globalElementType = []
    
    bondCount = 0
    globalBondType = []
    
    angleCount = 0
    globalAngleType = []
    
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
        oplsAtom.atomCount = oplsAtom.atomCount + 1
        self.atomID = oplsAtom.atomCount
        
        # Create dummy variable for possible bonds, bond orders, bond types
        self.bondList = []
        self.bondTypeList = []
        # Create dummy variable for possible angles
        self.angleList = []
        self.angleTypeList = []
        
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
        else:
            # Will not bond if given a huge radius
            warnings.warn("Unknown element type")
            self.radius = 1000
            
        # Organize atoms into lists of same elements
        for (offset, group) in enumerate(oplsAtom.elementList):
            if group[0] == self.getElement():
                oplsAtom.elementList[offset].append(self.getID())
                break
        else:
            oplsAtom.elementList.append([self.getElement()])
            oplsAtom.elementList[-1].append(self.getID())
    
    def shiftXYZ(self, shiftx, shifty, shiftz):
        #Move the atom by the prescribed amount
        self.x = self.x + shiftx
        self.y = self.y + shifty
        self.z = self.z + shiftz
        
    def assignBond(self, atom2):
        # Get id of atom being bonded to
        atom2ID = atom2.getID()
        
        # Make sure that atom has not already been bonded to. Add bond if 
        # it has not already been bonded.
        if self.bondList.count(atom2ID) == 0:
            self.bondList.append(atom2ID)
            
            # Always make initial assumption of only one bond type
            self.bondTypeList.append(1)
            
            # Return boolean true if bond was successful
            return True
            
        else:
            # Return boolean false if bond was not successful
            return False
            
    def breakAllBonds(self, atomList):
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
        self.bondTypeList.pop()
            
    def assignElementType(self):
        # Current OPLS code only optimized for C60 on pentacene where all
        # carbon atoms are aromatic and all hydrogen atoms are attached to
        # aromatic carbons

        if self.elementType == 'C':
            self.elementType = self.elementType + 'A'
            
        if self.elementType == 'H':
            self.elementType = self.elementType + 'A'
                
        # Increment dihedral count depending on number of carbons since each
        # carbon is bonded to exactly three other atoms
        if self.elementType == 'CA':
            oplsAtom.improperCount = oplsAtom.improperCount + 1
            
        return self.elementType
                    
    def assignAngle(self, i_atom, k_atom, angleTypeNum):
        # Get id's of atoms in angle
        i_ID = i_atom.getID()
        j_ID = self.getID()
        k_ID = k_atom.getID()
        
        # Create sorted angle with id's
        angle = [i_ID, k_ID]
        angle.sort()
        angle = [angle[0], j_ID, angle[1]]
        
        # Make sure that angle has not already been added. Add angle if 
        # it has not already been assigned.
        for (offset, myAngle) in enumerate(self.angleList):
            if myAngle[0] == angle[0] and myAngle[1] == angle[1] and myAngle[2] == angle[2]:
                # Return boolean false if angle was not successful
                return False                    
                break
        else:
            self.angleList.append(angle)
            self.angleTypeList.append(angleTypeNum)
            # Return boolean true if bond was successful
            return True            
   
    def assignCharge(self, atomList):
        # Set atom charge depending on element type
        elementType = self.getElementType()
        
        # Base charge for each element type in pentacene
        if elementType == 'CA':
            self.charge = -0.115
        elif elementType == 'HA':
            self.charge = 0.115
        elif elementType == 'CB':
            self.charge = 0.0
        elif self.charge == 0:
            warnings.warn("Unknown element type for charge determination")
            
        # Set carbons in fused rings to 0 (See Naphthalene Fusion C)
        # Occurs when carbon aromatically linked to 3 other aromatic carbons
        elementTypeNeighborList = self.next2ElementType('CA', atomList)
        if len(elementTypeNeighborList) > 3:
            self.charge = 0.0
        
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
        if len(elementTypeNeighborList) > 1:
            elementTypeNeighborList[0] = True
        return elementTypeNeighborList
   
    def getID(self):
        return self.atomID
        
    def getBondList(self):
        # Return list of bonds [Note: these should already be presorted]
        return self.bondList
        
    def getBondTypeList(self):
        # Return bond types in the same sequence as bondList
        return self.bondTypeList
        
    def getAngleList(self):
        # Return list of angles
        return self.angleList       
        
    def getAngleTypeList(self):
        # Return list of angles
        return self.angleTypeList     
    
    def getCount(self):
        return oplsAtom.atomCount
        
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
              ', ', self.z, sep='')
              
    def __eq__(self, other):
        return calcDistance(self, other) == 0
              

# Return list of elements sorted by element
def getElementList():
    return oplsAtom.elementList
    
# Return list of element types sorted by element
def getGlobalElementTypeList():
    return oplsAtom.globalElementType
    
# Return list of bond types
def getGlobalBondTypeList():
    return oplsAtom.globalBondType
    
# Return list of angle types
def getGlobalAngleTypeList():
    return oplsAtom.globalAngleType
    
# Establish when buckyball starts and how many there are
def setBuckyStart(atomNum, buckyCount, coarse):
    oplsAtom.buckyStart = atomNum
    oplsAtom.buckyCount = buckyCount
    oplsAtom.buckyCoarse = coarse
    
# Establish layered box boundary conditions for this data file
# Primarily used for repeating COF system
def setLayerBox(periodic, a, angle, n, space):
    # Set whether periodic boundary conditions are used
    oplsAtom.periodic = periodic
    
    # For lammmps, trigonal vectors established by using xy xz yz
    # A = (xhi-xlo,0,0); 
    # B = (xy,yhi-ylo,0); 
    # C = (xz,yz,zhi-zlo)
    oplsAtom.xlo = 0
    oplsAtom.xhi = oplsAtom.xlo + a
    oplsAtom.ylo = 0
    oplsAtom.yhi = oplsAtom.ylo + a*math.sin(math.radians(angle))
    oplsAtom.zlo = -5
    oplsAtom.zhi = 5
    oplsAtom.xy = a*math.cos(math.radians(angle))
    oplsAtom.xz = 0
    oplsAtom.yz = 0
    
    # Check if this is a many layer system
    if n != -1:
        oplsAtom.zlo = 0
        oplsAtom.zhi = n * space
        
# Establish triclinic box boundary conditions for this data file
# Currently designed for C60 on pentacene system
def setTriclinicBox(periodic, a, b, c, alpha, beta, gamma):
    # Set whether periodic boundary conditions are used
    oplsAtom.periodic = periodic
    
    # For lammmps, trigonal vectors established by using xy xz yz
    # A = (xhi-xlo,0,0); 
    # B = (xy,yhi-ylo,0); 
    # C = (xz,yz,zhi-zlo)
    oplsAtom.xlo = 0
    oplsAtom.ylo = 0
    oplsAtom.zlo = 0
    
    # Formula for converting (a,b,c,alpha,beta,gamma) to (lx,ly,lz,xy,xz,yz)
    # taken from online lammps help
    oplsAtom.xhi = a
    oplsAtom.xy = b*math.cos(math.radians(gamma))
    oplsAtom.xz = c*math.cos(math.radians(beta))
    oplsAtom.yhi = math.sqrt(b**2 - oplsAtom.xy**2)
    oplsAtom.yz = (b*c*math.cos(math.radians(alpha)) - oplsAtom.xy * oplsAtom.xz)/ oplsAtom.yhi
    oplsAtom.zhi = math.sqrt(c**2 - oplsAtom.xz**2 - oplsAtom.yz**2)
  
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
    if oplsAtom.periodic:
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

# Translate x,y,z coordinates across vector 1/A, as many times as 'multiple'
def translateVector1A(pos, multiple):
    x_a = pos['x'] + (oplsAtom.xhi - oplsAtom.xlo) * multiple
    y_a = pos['y'] + 0
    z_a = pos['z'] + 0
    return {'x':x_a, 'y':y_a, 'z':z_a}

# Translate x,y,z coordinates across vector 2/B, as many times as 'multiple'
def translateVector2B(pos, multiple):
    x_b = pos['x'] + oplsAtom.xy * multiple
    y_b = pos['y'] + (oplsAtom.yhi - oplsAtom.ylo) * multiple
    z_b = pos['z'] + 0
    return {'x':x_b, 'y':y_b, 'z':z_b}
    
# Translate x,y,z coordinates across vector 3/C, as many times as 'multiple'
def translateVector3C(pos, multiple):
    x_c = pos['x'] + oplsAtom.xz * multiple
    y_c = pos['y'] + oplsAtom.yz * multiple
    z_c = pos['z'] + (oplsAtom.zhi - oplsAtom.zlo) * multiple
    return {'x':x_c, 'y':y_c, 'z':z_c}
    
def calcAVGBondLength(atom1, atom2):
    radius1 = atom1.getRadius()
    radius2 = atom2.getRadius()
    avgBondLength = radius1 + radius2
    return avgBondLength
    
# Creates a bond between to atoms by adding to their bond list
def createBond(atom1, atom2):
    if atom1.assignBond(atom2) and atom2.assignBond(atom1):
        oplsAtom.bondCount = oplsAtom.bondCount + 1
        #print(oplsAtom.bondCount, atom1.getID(), atom2.getID())
    
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
    closeAtomList = [[None, 1000]] * numberCloseAtoms
    
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
    return closeAtomList
    
# Determine element type for all atoms. Compile global list
# Requires bonds
def determineElementTypes(atomList):
    for (offset, atom1) in enumerate(atomList):
        atom1.assignElementType()
        
        # Organize atoms into lists of same element types
        for (offset, group) in enumerate(oplsAtom.globalElementType):
            if group[0] == atom1.getElementType():
                oplsAtom.globalElementType[offset].append(atom1.getID())
                break
        else:
            oplsAtom.globalElementType.append([atom1.getElementType()])
            oplsAtom.globalElementType[-1].append(atom1.getID())
                    
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
                for (offset, pair) in enumerate(oplsAtom.globalBondType):
                    if pair[0] == bond_type[0] and pair[1] == bond_type[1]:
                        changeBondType(bondAtom, atom, offset+1)
                        break
                else:
                    oplsAtom.globalBondType.append(bond_type)
                    changeBondType(bondAtom, atom, len(oplsAtom.globalBondType))
                    
# Determine angles and angle types for all atoms. Compile global list of angle types
# Requires bond and element type information
def determineAnglesAndTypes(atomList):
    # Cycle through every atom
    for j_atom in atomList:
        # Get ID and all the bonds associated with the atom
        j_ID = j_atom.getID()
        j_bondList = j_atom.getBondList()
        
        # Cycle through every atom bonded to j to find potential k atoms
        for (offset, k_ID) in enumerate(j_bondList):
            # Cycle through every atom bonded to j (again), in order to 
            # find all possible i atoms
            for i_ID in j_bondList:
                # Check to make sure i_atom is not already k_atom
                if i_ID != k_ID and k_ID > i_ID:
                    oplsAtom.angleCount = oplsAtom.angleCount + 1
                    # Make sure that global angle type has not already been added. 
                    # Add angle if it has not already been assigned.
                    i_atom = atomList[i_ID-1]
                    k_atom = atomList[k_ID-1]
                    
                    i_type = i_atom.getElementType()
                    j_type = j_atom.getElementType()
                    k_type = k_atom.getElementType()
                    
                    angle_type = [i_type, k_type]
                    angle_type.sort()
                    angle_type = [angle_type[0], j_type, angle_type[1]]
                    
                    # Create running list of bond pair types
                    for (offset, globalAngle) in enumerate(oplsAtom.globalAngleType):
                        if globalAngle[0] == angle_type[0] and globalAngle[1] == angle_type[1] and globalAngle[2] == angle_type[2]:
                            j_atom.assignAngle(i_atom, k_atom, offset+1)
                            break
                    else:
                        oplsAtom.globalAngleType.append(angle_type)
                        j_atom.assignAngle(i_atom, k_atom, len(oplsAtom.globalAngleType))
                                  
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
                                oplsAtom.dihedralCount = oplsAtom.dihedralCount + 1

# Print system info as required by the intro to lammps data file
def getSystemInfo():
    # Compile number of items
    #Temporarily turn off all bonds
    #oplsAtom.bondCount = 0
    #oplsAtom.angleCount = 0
    #oplsAtom.dihedralCount = 0
    text = "LAMMPS Description" + "\n\n"
    text = text + str(oplsAtom.atomCount) + '  atoms\n'
    text = text + str(oplsAtom.bondCount) + '  bonds\n'
    text = text + str(oplsAtom.angleCount) + '  angles\n'
    text = text + str(oplsAtom.dihedralCount) + '  dihedrals\n'
    text = text + str(oplsAtom.improperCount) + '  impropers\n\n'
    
    # Compile number of each type of item
    atomTypeCount = len(getElementList())
    text = text + str(atomTypeCount) + '  atom types\n'
    bondTypeCount = len(getGlobalBondTypeList())
    text = text + str(bondTypeCount) + '  bond types\n'
    angleTypeCount = len(getGlobalAngleTypeList())
    text = text + str(angleTypeCount) + '  angle types\n'
    dihedralTypeCount = 1
    text = text + str(dihedralTypeCount) + '  dihedral types\n'
    improperTypeCount = 1
    text = text + str(improperTypeCount) + '  improper types\n\n'
    
    # Compile box information
    temp = '%.6f %.6f xlo xhi\n' % (oplsAtom.xlo, oplsAtom.xhi)
    text = text + temp
    temp = '%.6f %.6f ylo yhi\n' % (oplsAtom.ylo, oplsAtom.yhi)
    text = text + temp
    temp = '%.6f %.6f zlo zhi\n' % (oplsAtom.zlo, oplsAtom.zhi)
    text = text + temp
    temp = '%.6f %.6f %.6f xy xz yz\n' % (oplsAtom.xy, oplsAtom.xz, oplsAtom.yz)
    text = text + temp
    
    # Return compiled text
    text = text + '\n'
    return text
    
# Print mass as required by the lammps data file
def getMassInfo():
    # Compile list of masses
    text = "Masses\n\n"
    
    # Use list of element types to make mass list    
    globalElementType = getGlobalElementTypeList()
    
    # Organize atoms into lists of same elements
    # Assign mass (g/mol)
    for (offset, group) in enumerate(globalElementType):
        if group[0] == 'HA':
            mass = 1.00794
        elif group[0] == 'CA':
            mass = 12.0107
        elif group[0] == 'CB':
            mass = 720.64
        else:
            # Will hopefully be unmovable with high mass
            warntext = "Unknown element type: " + group[0]
            warnings.warn(warntext)
            mass = 100000.0
            
        temp = ' %d %.4f # %s\n' % (offset+1, mass, group[0])
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Print bond coeffs as required by the lammps data file
# NOTE: Ni currently ignored
def getPairCoeffsInfo():
    # Compile list of masses
    text = "Pair Coeffs\n# pair_style lj/cut\n"
    
    # Use list of element types to make mass list    
    globalElementType = getGlobalElementTypeList()
    
    # Organize atoms into lists of same elements
    # R0 = Angstroms
    # D0 = kcal/mol
    for (offset, group) in enumerate(globalElementType):
        if group[0] == 'HA':
            pair_style = 'lj/cut/coul/long'
            sigma = 2.420
            eps = 0.030
        elif group[0] == 'CA':
            pair_style = 'lj/cut/coul/long'
            sigma = 3.550
            eps = 0.070
        elif group[0] == 'CB':
            pair_style = 'nm/cut'
            sigma = 10.04805
            eps = 6.648
            n = 35.4877
            m = 8.8719
        else:
            warntext = "Unknown element type: " + group[0]
            warnings.warn(warntext)
            pair_style = 'lj/cut/coul/long'
            sigma = 1.0
            eps = 0.0
            
        # Already in lammps units: eps and sigma
        #eps = D0
        #sigma = R0/(2**(1/6))
        
        if pair_style == 'lj/cut/coul/long' or pair_style == 'lj/cut':
            temp = ' %d %s %.4f %.4f # %s\n' % (offset+1, pair_style, eps, sigma, group[0])
        elif pair_style == 'nm/cut':
            temp = ' %d %s %.4f %.4f %.4f %.4f# %s\n' % (offset+1, pair_style, eps, sigma, n, m, group[0])
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Print bond coefficients as required by the lammps data file
def getBondCoeffsInfo():
    # Add header. Bond_style info not used by input file
    text = "Bond Coeffs\n# bond_style harmonic\n"
    
    # Use list of fundamental bond pair types to make bond coeff list    
    globalBondType = getGlobalBondTypeList()
    
    # Cycle through list of bond pair types
    for (offset, pair) in enumerate(globalBondType):
        # Assign energy (kcal / mol * A^2) and equilibrium bond length (A)
        if pair[0] == 'CA' and pair[1] == 'HA':
            equilBond = 1.080 
            LAMMPSenergy = 367.0
        elif pair[0] == 'CA' and pair[1] == 'CA':
            equilBond = 1.400
            LAMMPSenergy = 469.0
        else:
            warntext = "Unknown bond type: " + pair
            warnings.warn(warntext)
            equilBond = 100000
            LAMMPSenergy = 100000
        
        # Already in lammps units
        
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
    globalAngleType = getGlobalAngleTypeList()
    
    # Cycle through list of angle types
    for (offset, angleType) in enumerate(globalAngleType):
        # Assign energy (energy/rad^2) and degrees
        if angleType[0] == 'CA' and angleType[1] == 'CA' and angleType[2] == 'CA':
            K = 63.00
            degrees = 120
        elif angleType[0] == 'CA' and angleType[1] == 'CA' and angleType[2] == 'HA':
            K = 35.00
            degrees = 120
        else:
            warntext = "Unknown angle type: " + angleType
            warnings.warn(warntext)
            K = 63.00
            degrees = 180
            
        # Already in lammps units
            
        temp = ' %d %f %.3f # %s\n' % (offset+1, K, degrees, angleType)
        text = text + temp

    # Return compiled text
    text = text + '\n'
    return text
    
# Print dihedral coefficients as required by the lammps data file
def getDihedralCoeffsInfo():
    # Add header. Dihedral_style info not used by input file
    text = "Dihedral Coeffs\n# dihedral_style opls\n"
    
    # Only one dihedral type is currently present in pentacene system so only
    # this one will be coded
    
    # Formula used in dihedral_style opls
    # E = K1/2 * (1 + cos(d)) + K2/2 * (1 - cos(2 * d)) + K3/2 * (1 + cos(3 * d)) + K4/2 * (1 - cos(4 * d))
    # K1 (energy)
    # K2 (energy)
    # K3 (energy)
    # K4 (energy)    
    
    # Diehedral for all aromatic carbon centers: 'X-CA-CA-Y'
    K1 = 0.0
    K2 = 7.250
    K3 = 0.0
    K4 = 0.0
    text = text + ' %d %.3f %.3f %.3f %.3f # %s\n' % (1, K1, K2, K3, K4, 'X-CA-CA-Y')
    
    # Return compiled text
    text = text + '\n'
    return text

# Print improper coefficients as required by the lammps data file
def getImproperCoeffsInfo():
    # Add header. improper_style info not used by input file
    text = "Improper Coeffs\n# improper_style cvff\n"
    
    # There is no distinction in formula between dihedrals and impropers in opls.
    # Formula used in dihedral_style opls
    # E = K1/2 * (1 + cos(d)) + K2/2 * (1 - cos(2 * d)) + K3/2 * (1 + cos(3 * d)) + K4/2 * (1 - cos(4 * d))
    # K1 (energy)
    # K2 (energy)
    # K3 (energy)
    # K4 (energy)
    
    # OPLS defines an aromatic improper ('Z-CA-X-Y') with:
    #K1 = 0.0
    K2 = 2.200
    #K3 = 0.0
    #K4 = 0.0
    
    # Lammps provides the following function for improper_style cvff:
    # E = K(1 + J*cos(n*d))
    # Therefore, we can use:
    K = K2 / 2
    J = -1
    n = 2
    
    # Since there is only one style in opls for planar impropers, always include
    text = text + ' %d %.3f %d %d # %s\n' % (1, K, J, n, 'Z-CA-X-Y')
    
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
    buckyCarbon = 0
    for atom in atomList:
        # Get all atom information
        atomID = atom.getID()
        element = atom.getElement()
        elementType = atom.getElementType()
        elementID = elementNumDict[element]
        
        #Set each C60 as a different molecule
        if element == 'CB' or atomID >= oplsAtom.buckyStart:
            if oplsAtom.buckyCoarse:
                atomsInBucky = 1
            else:
                atomsInBucky = 60
            moleculeID = (buckyCarbon) // atomsInBucky + 2
            buckyCarbon = buckyCarbon + 1
        else:
            moleculeID = 1
        charge = atom.getCharge()
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
    
    # Acquire global list of angle types
    globalAngleType = getGlobalAngleTypeList()
    
    angleNum = 0
    # Cycle through every atom
    for atom in atomList:
        # Get all angles and print in order
        angleList = atom.getAngleList()
        angleTypeList = atom.getAngleTypeList()
        
        # Cycle through every angle
        for (offset, angle) in enumerate(angleList):
            angleNum = angleNum + 1
            
            # Collect angle type number and element type composition
            angleTypeID = angleTypeList[offset]
            angle_type = globalAngleType[angleTypeID-1]
            
            text = text + ' %d %d %d %d %d # %s\n' % (angleNum, angleTypeID, angle[0], angle[1], angle[2], angle_type)
            
    # Return compiled text
    text = text + '\n'
    return text

# Systemically print all dihedrals in format required by lammps data file
def getDihedralListInfo(atomList):
    # Add header.
    text = "Dihedrals\n\n"
    
    # Cycle through every atom
    dihedralNum = 0
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
                                dihedralNum = dihedralNum + 1
                                
                                # Find element type for all atoms
                                i_atom = atomList[i_ID-1]
                                l_atom = atomList[l_ID-1]
                                i_type = i_atom.getElementType()
                                j_type = j_atom.getElementType()
                                k_type = k_atom.getElementType()
                                l_type = l_atom.getElementType()
                                
                                # Determine type of bond; Currently does not 
                                # check for bond order
                                
                                # Check if j,k atoms are aromatic carbons
                                if (j_type == 'CA' and k_type == 'CA'):
                                    # Only opls dihedral
                                    text = text + ' %d %d %d %d %d %d\n' % (dihedralNum, 1, i_ID, j_ID, k_ID, l_ID)
                                   
                                # Return warning if no dihedral found
                                else:
                                    text = text + ' %d %d %d %d %d %d # ERROR\n' % (dihedralNum, 1, i_ID, j_ID, k_ID, l_ID)
    
    # Return compiled text
    text = text + '\n'
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
        
        # If this is an atom bonded to three other atoms and it forms an improper,
        # add to list
        if i_type == 'CA':
            # Increment improper number
            improperNum = improperNum + 1
            
            # Get atom id and bondlist
            i_ID = i_atom.getID()
            i_bondList = i_atom.getBondList()
            
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
            
            # Define improper type like a proper dihedral so that it matchs opls function
            improper_type = [j_type, i_type, k_type, l_type]
            
            # There is currently only one improper type for opls model
            improperTypeID = 1
            
            text = text + ' %d %d %d %d %d %d # %s\n' % (improperNum, improperTypeID, j_ID, i_ID, k_ID, l_ID, improper_type)

    # Return compiled text
    text = text
    return text
