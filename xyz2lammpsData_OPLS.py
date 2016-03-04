"""
Author: Yaset Acevedo
Created: 4/2/14
Update Date: 4/14/14
Version: 1.3
"""

#Import modern print module and math module
from __future__ import print_function
from __future__ import division
import math

# Import opls force field atom module
import oplsFF

# *****************************************************************************
# Input Variables:
inputFile = 'Pn_supercell' # .xyz automatically appended

# Periodic bond check?
periodic = True

# Charge info available?
chargeInfo = True

# Setup periodic boundary conditions
# Triclinic parameters
# For layered box, set all angles to 90 and set b=c=0
a = 105.84
b = 130.46
c = 150.0

#bulk
#alpha = 101.9
#beta = 112.6
#gamma = 85.8

#thin
alpha = 86.7
beta = 81.4
gamma = 89.8

#Bucky atom start? What atomic number do the bucky balls start appearing
buckyStart = 10**7
buckyNum = 1
buckyCoarse = True

# For Layer box:
# First parameter is length of side for 2d system
# Second parameter is angle between v1 and v2
#a = 21.56
#angle = 90
# Input space between layers
space = 3.4
# *****************************************************************************

# Select input and output file
inputFile = inputFile + '.xyz'
outputFile = 'data.' + inputFile[0:-4] + '.txt'

# Import xyz atom file (e.g. from Avogadro)
f = open(inputFile, 'r')

# Check if there are many layers from the file name
n = inputFile.find('layers')
if n > -1:
    n = int(inputFile[n-1])

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
    newAtom = oplsFF.oplsAtom(atomProps[0], atomProps[1], atomProps[2], 
                           atomProps[3])
                           
    # Shift atoms so that no negative z positions are used.  Offset x and y
    # positions appropriately
    newAtom.shiftXYZ(1.514880, 0.577873, 10)
        
    #Add atom to total list of atoms
    #newAtom.displayAtom()
    atomList.append(newAtom)
    
    # If atom is in the buckyball, add to list of buckyball atoms
    if (offset+1 >= buckyStart):
        buckyList.append(newAtom)

# Close file to save memory    
f.close()

# Import boundary conditions into global class variables
if alpha == 90 and beta == 90 and gamma == 90:
    oplsFF.setLayerBox(periodic, a, alpha, n, space)
else:
    oplsFF.setTriclinicBox(periodic,a,b,c,alpha,beta,gamma)
    
# Tell program when bucky ball atoms begin
oplsFF.setBuckyStart(buckyStart,buckyNum,buckyCoarse)
        
# Cycle through all pairs of atoms and determine bonds for all cases
print('Assigning bonds: All others')
searchList = list(atomList)
assignedAtoms = 0

while len(searchList) > 0:
    # Produce counter to track how many atoms have been assigned bonds
    if (assignedAtoms%50 == 0):
        print(str(assignedAtoms) + ' / ' + str(numAtoms))
    
    # Gather information on current atom including element and number of current bonds
    j_atom = searchList[0]
    j_element = j_atom.getElement()
    j_bonds = j_atom.getBondList()
    j_curBond = len(j_bonds)
    
    # Guess max number of bonds from element
    if j_element == 'H':
        j_maxBond = 1
    if j_element == 'C':
        j_maxBond = 3
    if j_element == 'CB':
        j_maxBond = 0
        
    #Determine if all bonds have already been found for current j_atom
    if j_curBond == j_maxBond:
        searchList.pop(0)
        assignedAtoms = assignedAtoms + 1
        continue
        
    for (offset, k_atom) in enumerate(searchList[1:]):
        # Calculate avg bond length for particular atom pair and calculate 
        # distance between the atoms
        avgBondLength = oplsFF.calcAVGBondLength(j_atom, k_atom)
        atomDist = oplsFF.calcDistance(j_atom, k_atom)
        
        # If atom distance is within 25% of average bond length, create a bond
        error = math.fabs(avgBondLength - atomDist) / avgBondLength * 100
               
        if error < 30:
            #print((avgBondLength - atomDist) / avgBondLength * 100)
            oplsFF.createBond(j_atom, k_atom)
            
            # Increment number of bonds for each succesful bonds
            j_curBond = j_curBond + 1
            
            #Determine if you can stop finding bonds for current j_atom
            if j_curBond == j_maxBond:
                searchList.pop(0)
                assignedAtoms = assignedAtoms + 1
                break

print('Found all bonds') 
            
# Cycle through all atoms and determine type of opls atom
# Requires bond information to be accurate
oplsFF.determineElementTypes(atomList)

# Cycle through all atoms and determine bond types
# Requires element type information to be accurate
# Also used to count number of impropers (Each C_R is an improper)
oplsFF.determineBondTypes(atomList) 

# Cycle through all atoms and determine charges 
# Requires CHELPG charges from Gaussian
# Requires element type information to be accurate
# Requires bonds to be accurate
if chargeInfo:
    oplsFF.determineCharges(atomList) 

# Cycle through all atoms and determine dihedral count
# Must be run before outputting system info
# Requires bonds to be accurate
oplsFF.determineDihedralCount(atomList) 
    
# Cycle through all atoms and determine angles
# Requires bond information to be accurate
oplsFF.determineAnglesAndTypes(atomList)
#for (offset, atom1) in enumerate(atomList):
#    atom1.assignAngles()
    #print(atom1.getAngleList()) 
    
# Fix bond length and angle parameters for buckyball atoms
#oplsFF.fixBuckyParams(atomList, buckyList)
    
# Print intro to lammps file
# Export to data file
f = open(outputFile, 'w')
txt = oplsFF.getSystemInfo()
txt = txt + oplsFF.getMassInfo()
txt = txt + oplsFF.getPairCoeffsInfo()
txt = txt + oplsFF.getBondCoeffsInfo()
txt = txt + oplsFF.getAngleCoeffsInfo()
txt = txt + oplsFF.getDihedralCoeffsInfo()
txt = txt + oplsFF.getImproperCoeffsInfo()
txt = txt + oplsFF.getAtomListInfo(atomList)
txt = txt + oplsFF.getBondListInfo(atomList)
txt = txt + oplsFF.getAngleListInfo(atomList)
txt = txt + oplsFF.getDihedralListInfo(atomList)
txt = txt + oplsFF.getImproperListInfo(atomList)
print(txt,file=f)
f.close()
print('Done converting xyz to lammps data file')
 