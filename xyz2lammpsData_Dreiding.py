"""
Author: Yaset Acevedo
Created: 10/10/14
Update Date: 11/8/15
Version: 2.0
"""

#Import modern print module and math module
from __future__ import print_function
from __future__ import division
import math

# Import dreiding force field atom module
import dreidingFF

# *****************************************************************************
# Input Variables:
inputFile = 'hhtp_bp_x1y1z1' # .xyz automatically appended
#inputFile = 'PC_PH_periodic' # .xyz automatically appended

# Periodic bond check?
periodic = True

# Wrap atoms?
wrapAtoms = True

# Automatic charge info available?
chargeInfo = False

# External charge info available?
externalCharge = True
chargeFile = inputFile + '_charge.txt'

# Setup periodic boundary conditions
# Triclinic parameters
# For layered box, set all angles to 90 and set b=c=0
xlo = 0
xhi = 0
ylo = 0
yhi = 0
a = xhi - xlo
b = yhi - ylo
c = 3.4
alpha = 90
beta = 90
gamma = 90

#Bucky atom start? What atomic number do the bucky balls start appearing
C60_start = 10000000

# For Layer box:
# First parameter is length of side for 2d system
# Second parameter is angle between v1 and v2
#a = 21.56
#angle = 90
# Input space between layers
space = 3.4

# Bond by nearest neighbors? Otherwise, bond by approximate bond length
nearBond = True
# *****************************************************************************

# Check if this is a supercell from the file name
repeat_x = 1
repeat_y = 1
x = inputFile.find('x')
if x > -1:
    repeat_x = int(inputFile[x+1:x+2])

y = inputFile.find('y')
if y > -1:
    repeat_y = int(inputFile[y+1:y+2])

a = a*repeat_x
b = b*repeat_y

# Check if there are many layers from the file name
n = inputFile.find('z')
if n > -1:
    n = int(inputFile[n+1:n+3])

# Select input and output file
inputFile = inputFile + '.xyz'
outputFile = 'data.' + inputFile[0:-4] + '.txt'

# Import xyz atom file (e.g. from Avogadro)
f = open(inputFile, 'r')

# If external charge info available, import file (e.g. pulled from Gaussian log file)
if externalCharge == True:
    print('Will assign external charge during import')
    f_c = open(chargeFile, 'r')
    f_c.readline()
    f_c.readline()

# Record number of atoms from first line
firstLine = f.readline()
numAtoms = int(firstLine.strip())

# Record element and x,y,z data for all atoms
print('Importing from xyz file')
f.readline()
atomList = []
buckyList = []
for (offset, line) in enumerate(f):
    atomProps = line.split()
    # If Ge detected, convert this to a coarse grained C60
    if atomProps[0] == 'Ge':
        atomProps[0] = 'CB'
    newAtom = dreidingFF.dreidingAtom(atomProps[0], atomProps[1], atomProps[2],
                           atomProps[3])

    # Shift atoms so that simulations are less likely to cross the z plane.
    newAtom.shiftXYZ(3.5, 3, 0)

    # If external charge info available, assign charge to atom
    if externalCharge == True:
        chargeProps = f_c.readline().split()
        newAtom.assignManualCharge(chargeProps[1])

    #Add atom to total list of atoms
    #newAtom.displayAtom()
    atomList.append(newAtom)

    # If atom is in the buckyball, add to list of buckyball atoms
    if (offset+1 >= C60_start):
        buckyList.append(newAtom)

# Close file to save memory
f.close()

# Import boundary conditions into global class variables
if alpha == 90 and beta == 90 and gamma == 90:
    dreidingFF.setLayerBox(periodic, a, b, alpha, n, space)
else:
    dreidingFF.setTriclinicBox(periodic,a,b,c,alpha,beta,gamma)

# Cycle through all atoms and wrap atom xyz coordinates inside box
if wrapAtoms == True:
    for (offset, atom) in enumerate(atomList):
        atom.wrapXYZ()

# Bond using nearest neighbors
if nearBond == True:
    searchList = list(atomList)

    # Determine bonds for special element cases
    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all Hydrogen atoms and bond to closest 1 atoms
        if group[0] == 'H':
            print('Assigning bonds: Hydrogen')
            for (offset, atomID) in enumerate(group[1:]):
                if (offset%50 == 0):
                    print(str(offset) + ' / ' + str(len(group[1:])) + ' Hydrogens')
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # If atom has already been analyzed, skip
                if atom1 not in searchList:
                    continue

                # Gather information on current atom including element and number of current bonds
                j_bonds = atom1.getBondList()
                j_curBond = len(j_bonds)
                j_maxBond = atom1.getMaxBond()

                #Determine if all bonds have already been found for current j_atom
                if j_curBond == j_maxBond:
                    searchList.remove(atom1)
                    continue

                # Find closest 1 atoms for H
                closestAtomList = dreidingFF.findClosest(atom1, searchList, j_maxBond-j_curBond)
                # Sort list numerically for convenience
                closestAtomList.sort()

                # Bond to closest 1 atoms
                for (atom2, dist) in closestAtomList:
                    dreidingFF.createBond(atom1, atom2)

                    # Gather information on potential bonded atom including element and number of current bonds
                    k_bonds = atom2.getBondList()
                    k_curBond = len(k_bonds)
                    k_maxBond = atom2.getMaxBond()

                    #Determine if all bonds have been found for k_atom
                    if k_curBond == k_maxBond:
                        searchList.remove(atom2)
                        continue

                searchList.remove(atom1)

    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all Oxygen atoms and bond to closest 2 atoms
        if group[0] == 'O':
            print('Assigning bonds: Oxygen')
            for (offset, atomID) in enumerate(group[1:]):
                if (offset%50 == 0):
                    print(str(offset) + ' / ' + str(len(group[1:])) + ' Oxygens')
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # If atom has already been analyzed, skip
                if atom1 not in searchList:
                    continue

                # Gather information on current atom including element and number of current bonds
                j_bonds = atom1.getBondList()
                j_curBond = len(j_bonds)
                j_maxBond = atom1.getMaxBond()

                #Determine if all bonds have already been found for current j_atom
                if j_curBond == j_maxBond:
                    print(searchList.index(atom1))
                    searchList.remove(atom1)
                    continue

                # Find closest 2 atoms for O
                closestAtomList = dreidingFF.findClosest(atom1, searchList, j_maxBond-j_curBond)
                # Sort list numerically for convenience
                closestAtomList.sort()

                # Bond to closest 2 atoms
                for (atom2, dist) in closestAtomList:
                    dreidingFF.createBond(atom1, atom2)

                    # Gather information on potential bonded atom including element and number of current bonds
                    k_bonds = atom2.getBondList()
                    k_curBond = len(k_bonds)
                    k_maxBond = atom2.getMaxBond()

                    #Determine if all bonds have been found for k_atom
                    if k_curBond == k_maxBond:
                        searchList.remove(atom2)
                        continue

                searchList.remove(atom1)

    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all boron atoms and bond to closest 3 atoms
        if group[0] == 'B':
            print('Assigning bonds: Boron')
            for atomID in group[1:]:
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # If atom has already been analyzed, skip
                if atom1 not in searchList:
                    continue

                # Gather information on current atom including element and number of current bonds
                j_bonds = atom1.getBondList()
                j_curBond = len(j_bonds)
                j_maxBond = atom1.getMaxBond()

                #Determine if all bonds have already been found for current j_atom
                if j_curBond == j_maxBond:
                    searchList.remove(atom1)
                    continue

                # Find closest 3 atoms for Boron
                closestAtomList = dreidingFF.findClosest(atom1, searchList, j_maxBond-j_curBond)
                # Sort list numerically for convenience
                closestAtomList.sort()

                # Bond to closest 3 atoms
                for (atom2, dist) in closestAtomList:
                    dreidingFF.createBond(atom1, atom2)

                    # Gather information on potential bonded atom including element and number of current bonds
                    k_bonds = atom2.getBondList()
                    k_curBond = len(k_bonds)
                    k_maxBond = atom2.getMaxBond()

                    #Determine if all bonds have been found for k_atom
                    if k_curBond == k_maxBond:
                        searchList.remove(atom2)
                        continue

                searchList.remove(atom1)

    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all Ni(II) atoms and bond to closest 4 atoms
        if group[0] == 'Ni':
            print('Assigning bonds: Nickel')
            for atomID in group[1:]:
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # If atom has already been analyzed, skip
                if atom1 not in searchList:
                    continue

                # Gather information on current atom including element and number of current bonds
                j_bonds = atom1.getBondList()
                j_curBond = len(j_bonds)
                j_maxBond = atom1.getMaxBond()

                #Determine if all bonds have already been found for current j_atom
                if j_curBond == j_maxBond:
                    searchList.remove(atom1)
                    continue

                # Find closest 4 atoms for Ni(II)
                closestAtomList = dreidingFF.findClosest(atom1, searchList, j_maxBond-j_curBond)
                # Sort list numerically for convenience
                closestAtomList.sort()

                # Bond to closest 4 atoms
                for (atom2, dist) in closestAtomList:
                    dreidingFF.createBond(atom1, atom2)

                    # Gather information on potential bonded atom including element and number of current bonds
                    k_bonds = atom2.getBondList()
                    k_curBond = len(k_bonds)
                    k_maxBond = atom2.getMaxBond()

                    #Determine if all bonds have been found for k_atom
                    if k_curBond == k_maxBond:
                        searchList.remove(atom2)
                        continue

                searchList.remove(atom1)

    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all Carbon atoms and bond to closest 3 atoms
        if group[0] == 'C':
            print('Assigning bonds: Carbon')
            for (offset, atomID) in enumerate(group[1:]):
                if (offset%50 == 0):
                    print(str(offset) + ' / ' + str(len(group[1:])) + ' Carbons')
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # If atom has already been analyzed, skip
                if atom1 not in searchList:
                    continue

                # Gather information on current atom including element and number of current bonds
                j_bonds = atom1.getBondList()
                j_curBond = len(j_bonds)
                j_maxBond = atom1.getMaxBond()

                #Determine if all bonds have already been found for current j_atom
                if j_curBond == j_maxBond:
                    searchList.remove(atom1)
                    continue

                # Find closest 3 atoms for C
                closestAtomList = dreidingFF.findClosest(atom1, searchList, j_maxBond-j_curBond)
                # Sort list numerically for convenience
                closestAtomList.sort()

                # Bond to closest 3 atoms
                for (atom2, dist) in closestAtomList:
                    #print(dist)
                    dreidingFF.createBond(atom1, atom2)

                    # Gather information on potential bonded atom including element and number of current bonds
                    k_bonds = atom2.getBondList()
                    k_curBond = len(k_bonds)
                    k_maxBond = atom2.getMaxBond()

                    #Determine if all bonds have been found for k_atom
                    if k_curBond == k_maxBond:
                        searchList.remove(atom2)
                        continue

                searchList.remove(atom1)

    for (offset, group) in enumerate(dreidingFF.getElementList()):
        # Cycle through all C60 atoms and break all bonds
        if group[0] == 'CB':
            print('Breaking bonds: Buckyball')
            for (offset, atomID) in enumerate(group[1:]):
                if (offset%1 == 0):
                    print(str(offset) + ' / ' + str(len(group[1:])) + ' C60s')
                # Get atom1 reference
                atom1 = atomList[atomID-1]

                # Break all bonds with C60
                atom1.breakAllBonds()

# Bond by guessing bond distance between atoms and seeing if average bond
# distance is met
else:
    # Cycle through all pairs of atoms and determine bonds for all cases
    print('Assigning bonds: All others')
    searchList = list(atomList)
    assignedAtoms = 0
    index = 0

    # Stop while loop when all atoms have been assigned
    # Problem can occur that Nitrogen can have 2 or 3 bonds

    #while len(searchList) > 0:
    while assignedAtoms < numAtoms:
        # Produce counter to track how many atoms have been assigned bonds
        if (assignedAtoms%50 == 0):
            print(str(assignedAtoms) + ' / ' + str(numAtoms))

        # Gather information on current atom including element and number of current bonds
        j_atom = searchList[index]
        j_element = j_atom.getElement()
        j_bonds = j_atom.getBondList()
        j_curBond = len(j_bonds)

        # Guess max number of bonds from element
        if j_element == 'H':
            j_maxBond = 1
        if j_element == 'C':
            j_maxBond = 3
        if j_element == 'B':
            j_maxBond = 3
        if j_element == 'O':
            j_maxBond = 2
        if j_element == 'Ni':
            j_maxBond = 4
        if j_element == 'N':
            j_maxBond = 3
        if j_element == 'CB':
            j_maxBond = 0

        #Determine if all bonds have already been found for current j_atom
        if j_curBond == j_maxBond:
            searchList.pop(0)
            assignedAtoms = assignedAtoms + 1
            continue

        #print(j_element)
        #Loop through all atoms, including current j_atom
        for (offset, k_atom) in enumerate(searchList):
            #If on prior atom on list, check if all bonds have already been found
            #if (offset < index):
                #print(offset)

            #If on current j_atom (Self), then skip to next atom
            if (offset == index):
                continue

            # Calculate avg bond length for particular atom pair and calculate
            # distance between the atoms
            avgBondLength = dreidingFF.calcAVGBondLength(j_atom, k_atom)
            atomDist = dreidingFF.calcDistance(j_atom, k_atom)

            # If atom distance is within 25% of average bond length, create a bond
            error = math.fabs(avgBondLength - atomDist) / avgBondLength * 100

            if error < 35:
                #print(math.fabs(avgBondLength - atomDist) / avgBondLength * 100)
                dreidingFF.createBond(j_atom, k_atom)

                # Increment number of bonds for each succesful bonds
                j_curBond = j_curBond + 1

                #Determine if you can stop finding bonds for current j_atom
                if j_curBond == j_maxBond:
                    searchList.pop(0)
                    assignedAtoms = assignedAtoms + 1
                    break
        else:
            # If atom is not N, it should have found all its neighbors
            if j_element != 'N':
                print(j_element)
                print("ERROR")
            # If atom is N, assume it successfully found all of its neighbors
            else:
                assignedAtoms = assignedAtoms + 1
            index = index + 1

# Check if all bonds found
print('Checking bonds: All')
searchList = list(atomList)
assignedAtoms = 0
index = 0

while len(searchList) > 0:
    # Produce counter to track how many atoms have been assigned bonds
    if (assignedAtoms%50 == 0):
        print(str(assignedAtoms) + ' / ' + str(numAtoms))

    # Gather information on current atom including element and number of current bonds
    j_atom = searchList[index]
    j_element = j_atom.getElement()
    j_bonds = j_atom.getBondList()
    j_curBond = len(j_bonds)

    # Guess max number of bonds from element
    if j_element == 'H':
        j_maxBond = 1
    if j_element == 'C':
        j_maxBond = 3
    if j_element == 'B':
        j_maxBond = 3
    if j_element == 'O':
        j_maxBond = 2
    if j_element == 'Ni':
        j_maxBond = 4
    if j_element == 'N':
        j_maxBond = 2
    if j_element == 'CB':
        j_maxBond = 0

    #Determine if all bonds have already been found for current j_atom
    if j_element == 'N':
        if j_curBond >= j_maxBond:
            assignedAtoms = assignedAtoms + 1
        else:
            print('Error in ' + str(j_element) + '. Expected: ' + str(j_maxBond) + ' Have: ' + str(j_curBond))
    elif j_element == 'C':
        if j_curBond >= j_maxBond:
            assignedAtoms = assignedAtoms + 1
        else:
            print('Warning in ' + str(j_element) + '. Expected: ' + str(j_maxBond) + ' Have: ' + str(j_curBond)+ '. Triple bonds should be present')
    elif j_curBond == j_maxBond:
        assignedAtoms = assignedAtoms + 1
    else:
        print('Error in ' + str(j_element) + '. Expected: ' + str(j_maxBond) + ' Have: ' + str(j_curBond))

    searchList.pop(0)
print('Correct bonds:')
print(str(assignedAtoms) + ' / ' + str(numAtoms))

# Cycle through all atoms and determine periodic bonds
if periodic:
    print('Periodic bonds:')
    txt2 = dreidingFF.determinePeriodicBonds(atomList)
    outputFile2 = inputFile[0:-4] + '_periodicBonds.txt'

    f = open(outputFile2, 'w')
    print(txt2,file=f)
    f.close()

# Cycle through all atoms and determine type of Dreiding atom
# Requires bond information to be accurate
for (offset, atom1) in enumerate(atomList):
    atom1.assignElementType()

# Cycle through all atoms and determine bond orders
# Requires element type information to be accurate
dreidingFF.determineBondOrders(atomList)

# Cycle through all atoms and determine bond types
# Requires element type information to be accurate
# Also used to count number of impropers (Each C_R is an improper)
dreidingFF.determineBondTypes(atomList)

# Cycle through all atoms and determine charges
# Requires CHELPG charges from Gaussian
# Requires element type information to be accurate
# Requires bonds to be accurate
if chargeInfo:
    dreidingFF.determineCharges(atomList)

# Cycle through all atoms and determine dihedral count
# Must be run before outputting system info
# Requires bonds to be accurate
dreidingFF.determineDihedralCount(atomList)

# Cycle through all atoms and determine improper count
# Must be run before outputting system info
# Requires bonds to be accurate
dreidingFF.determineImproperCount(atomList)

# Cycle through all atoms and determine angles
# Requires bond information to be accurate
for (offset, atom1) in enumerate(atomList):
    atom1.assignAngles(atomList)
    #print(atom1.getAngleList())

# Fix bond length and angle parameters for buckyball atoms
#dreidingFF.fixBuckyParams(atomList, buckyList)

# Print intro to lammps file
# Export to data file
f = open(outputFile, 'w')
txt = dreidingFF.getSystemInfo()
txt = txt + dreidingFF.getMassInfo()
txt = txt + dreidingFF.getPairCoeffsInfo()
txt = txt + dreidingFF.getBondCoeffsInfo()
txt = txt + dreidingFF.getAngleCoeffsInfo()
txt = txt + dreidingFF.getImproperCoeffsInfo()
txt = txt + dreidingFF.getDihedralCoeffsInfo()
txt = txt + dreidingFF.getAtomListInfo(atomList)
txt = txt + dreidingFF.getBondListInfo(atomList)
txt = txt + dreidingFF.getAngleListInfo(atomList)
txt = txt + dreidingFF.getDihedralListInfo(atomList)
txt = txt + dreidingFF.getImproperListInfo(atomList)
print(txt,file=f)
f.close()
print('Done converting xyz to lammps data file')

