#Generate LAMMPS input file for muscovite mica based on INTERFACE coordinates and potential parameters
#Include bond and angle data

#ARGUMENTS: in out OR nx ny nz in out

import sys
from numpy import array,zeros,random,append
#from keyboard import keyboard

##################################
## Parse command line arguments ##
##################################

args = sys.argv
del args[0]

if(len(args)not in [3,6]):
    print("Need 3 or 6 arguments: struct ff out OR nx ny nz struct ff out")
    quit()


elif (len(args)==5):
    for i in range(0,3):
        args[i]=int(args[i])

#Number of unit cells in each direction
cells=[1,1,1]
if (len(args)==6):
    print("Using command line arguments: {} {} {}".format(args[0],args[1],args[2]))
    cells[0]=int(args[0])
    cells[1]=int(args[1])
    cells[2]=int(args[2])
    
#Input & output files
inLoc=args[len(args)-3]
ffLoc=args[len(args)-2]
outLoc=args[len(args)-1]

#Read unit cell scaled coordinates from INTERFACE data file
cellFile=open(inLoc)
ffFile=open(ffLoc)

#################################
## Declare arrays & basic info ##
#################################

names=[]
elements=[]
unitCoords=[]
charges=[]
types=[]

elements=["o","a","s","k","h"]
masses={"o":39.098300,"a":26.981539,"s":28.086000,"k":39.098300,"h":1.008000}

#For inter-cell bonds
atomCoords=[]

#######################
## READ FORCE FIELD  ##
#######################

copyText=''
cellBonds=[]
cellAngles=[]

section="header"
for line in ffFile:
    if("Masses" in line):
        section="copy"
        copyText+=line
    elif("Bonds" in line):
        section="bonds"
    elif("Angles" in line):
        section="angles"
    elif ("Dihedrals" in line or "Impropers" in line or "Atoms" in line):
        section="other"
    
    elif(section=="header"):
        
        if("atoms" in line):
            nAtomsPerCell=int(line.split()[0])
        elif("bonds" in line):
            nBondsPerCell=int(line.split()[0])
        elif("angles" in line):
            nAnglesPerCell=int(line.split()[0])
    elif(section=="copy"):
        copyText+=line
    elif(section=="bonds"):
        if(line!='\n'):
            cellBonds.append([int(i) for i in line.split()])
    elif(section=="angles"):
        if(line!='\n'):
            cellAngles.append([int(i) for i in line.split()])

######################
##  Read input file ##
######################

#PBC info
pbc=cellFile.readline().split()
pbc=[float(pbc[i]) for i in range(1,4)]
a=pbc[:]
xlo=0
xhi=pbc[0]*cells[0]
ylo=0
yhi=pbc[1]*cells[1]
zlo=0
zhi=pbc[2]*cells[2]

for i in range(0,nAtomsPerCell):
    line=cellFile.readline().split()
    if(line[0] not in names):
        names.append(line[0])
    charges.append(float(line[4]))
    #Keep track of which atoms are which type
    types.append(names.index(line[0])+1)
    #Atom coordinates
    unitCoords.append([float(j) for j in [line[1],line[2],line[3]]])

######################
## Atom Coordinates ##
######################

for i in range(0,cells[0]):
    for j in range(0,cells[1]):
        for k in range(0,cells[2]):
          for atom in unitCoords:
            atomCoords.append([atom[0]+i*a[0],atom[1]+j*a[1],atom[2]+k*a[2]])
            
#######################
## Important Numbers ##
#######################

#Number of unit cells
nUnitCells=cells[0]*cells[1]*cells[2]

#Number of atoms per unit cell
if(len(unitCoords)!=nAtomsPerCell):
    print("Discrepency between structure and force field data.")
    print("Structure: {} atoms".format(len(unitCoords)))
    print("Force field: {} atoms".format(nAtomsPerCell))
    quit()

#Total number of atoms
nAtoms=nUnitCells*nAtomsPerCell

nBonds=nBondsPerCell*nUnitCells
nAngles=nAnglesPerCell*nUnitCells

nAtomTypes=len(names)

#################
## atomNumList ##
#################

atomNumList=zeros([cells[0],cells[1],cells[2],nAtomsPerCell]).astype(int)

n=1
for i in range(0,cells[0]):
    for j in range(0,cells[1]):
        for k in range(0,cells[2]):
            for atom in range(0,nAtomsPerCell):
                atomNumList[i][j][k][atom]=int(n)
                n+=1

####################
## MAKE BOND LIST ##
####################

bondList=[]

currentCellNum=0
for i in range(0,cells[0]):
    for j in range(0,cells[1]):
        for k in range(0,cells[2]):
            for count in range(0,nBondsPerCell):
                bondList.append([cellBonds[count][0]+currentCellNum*nBondsPerCell,cellBonds[count][1],cellBonds[count][2]+currentCellNum*nAtomsPerCell,cellBonds[count][3]+currentCellNum*nAtomsPerCell])
            currentCellNum+=1

##########################
## FIX CROSS-CELL BONDS ##
##########################

print("Adjusting bonds")

bondCross=[]

#Identify cross-cell bonds
for bond in bondList:
    crossDimension=[0,0,0]
    #Coordinates of two atoms involved in bond
    bondCoords=[]
    bondCoords.append(atomCoords[bond[2]-1])
    bondCoords.append(atomCoords[bond[3]-1])
    for dim in range(0,3):
        #Identify bonds which span more than half the unit cell for each dimension
        if(abs(bondCoords[0][dim]-bondCoords[1][dim])>a[dim]/2):
            crossDimension[dim]=1
    bondCross.append(crossDimension)
    
    
#Fix bonds
bondNum=0
i=[0,0,0]
for i[0] in range(0,cells[0]):
    for i[1] in range(0,cells[1]):
        for i[2] in range(0,cells[2]):
            for cellBond in range(0,nBondsPerCell):
                #Coordinates of two atoms involved in bond
                bondCoords=[]
                bondCoords.append(atomCoords[bondList[bondNum][2]-1])
                bondCoords.append(atomCoords[bondList[bondNum][3]-1])
                
                #Determine which cell to move the bond to
                newCell=i[:] #This syntax copies i so that changing newCell doesn't affect i, which is the default behavior when equating lists in Python.
                
                changedBond=False
                if(1 in bondCross[bondNum]):
                    for dim in range(0,3):
                        #If this is a cross-cell bond
                        if(bondCross[bondNum][dim]==1):
                            
                            #Determine which atom in the bond to change
                            #Want to change the one with the smaller coord in this dimension
                            #But don't want to run this part more than once
                            if(not changedBond):
                                #Remember that this bond is being changed in at least one dimension
                                changedBond=True
                                
                                #Keep the one with larger value in the first tested dimension fixed
                                if(bondCoords[1][dim]>bondCoords[0][dim]):
                                    change=0
                                    fixed=1
                                else:
                                    change=1
                                    fixed=0
                            
                            #Determine which way to move bond
                            if(bondCoords[fixed][dim]>bondCoords[change][dim]):
                                #Add
                                #if this is not the last cell in the given dimension
                                if(newCell[dim]<cells[dim]-1):
                                    newCell[dim]+=1
                                #if it is the last cell in the given dimension
                                else:
                                    #move to the beginning of that dimension
                                    newCell[dim]=0
                            else:
                                #Subtract
                                #if this is not the first cell in the given dimension
                                if(newCell[dim]>0):
                                    newCell[dim]-=1
                                #if it is the first cell in the given dimension
                                else:
                                    #move to the end of that dimension
                                    newCell[dim]=cells[dim]-1
                            
                            #Equivalent atom number in 1st unit cell
                            relAtomNum=bondList[bondNum][change+2]%nAtomsPerCell
                            
                    #Change the atom number appropriately
                    bondList[bondNum][change+2]=atomNumList[newCell[0]][newCell[1]][newCell[2]][relAtomNum-1]
                        
                bondNum+=1

######################
## MAKE ANGLES LIST ##
######################

angleList=[]

currentCellNum=0
for i in range(0,cells[0]):
    for j in range(0,cells[1]):
        for k in range(0,cells[2]):
            for count in range(0,nAnglesPerCell):
                angleList.append([cellAngles[count][0]+currentCellNum*nAnglesPerCell,cellAngles[count][1],cellAngles[count][2]+currentCellNum*nAtomsPerCell,cellAngles[count][3]+currentCellNum*nAtomsPerCell,cellAngles[count][4]+currentCellNum*nAtomsPerCell])
            currentCellNum+=1

###########################
## FIX CROSS-CELL ANGLES ## ---- This is not done yet.
########################### -- Do the same as for bonds but between 1&2 and 2&3

print("Adjusting angles")

angleCross=[]

#Identify cross-cell angles
for angle in angleList:
    crossDimension=[[],[]]
    angleCoords=[]
    angleCoords.append(atomCoords[angle[2]-1])
    angleCoords.append(atomCoords[angle[3]-1])
    angleCoords.append(atomCoords[angle[4]-1])
    for dim in range(0,3):
        #Identify angles which span more than half the unit cell for each dimension
        
        #1-2
        if(abs(angleCoords[0][dim]-angleCoords[1][dim])>a[dim]/2):
            crossDimension[0].append(1)
        else:
            crossDimension[0].append(0)
        #3-2
        if(abs(angleCoords[2][dim]-angleCoords[1][dim])>a[dim]/2):
            crossDimension[1].append(1)
        else:
            crossDimension[1].append(0)
    angleCross.append(crossDimension)
    
    
#Fix angles
angleNum=0
i=[0,0,0]
for i[0] in range(0,cells[0]):
    for i[1] in range(0,cells[1]):
        for i[2] in range(0,cells[2]):
            for cellAngle in range(0,nAnglesPerCell):
                #First or second pair of atoms in angle
                pairNum=0
                
                #Keep middle atom fixed, move others if necessary
                for change in [0,2]:
                    #Coordinates of pair involved in angle
                    angleCoords=[]
                    angleCoords.append(atomCoords[angleList[angleNum][3]-1]) #Fixed - Middle
                    angleCoords.append(atomCoords[angleList[angleNum][2+change]-1]) #Change
                    
                    if(1 in angleCross[angleNum][pairNum]):
                        #Determine which cell to move the angle to
                        newCell=i[:]
                        for dim in range(0,3):
                            #If this is a cross-cell angle
                            if(angleCross[angleNum][pairNum][dim]==1):
                                #Want to keep the middle atom fixed
                                
                                #Determine which way to move angle
                                if(angleCoords[0][dim]>angleCoords[1][dim]):
                                    #Add
                                    #if this is not the last cell in the given dimension
                                    if(newCell[dim]<cells[dim]-1):
                                        newCell[dim]+=1
                                    #if it is the last cell in the given dimension
                                    else:
                                        #move to the beginning of that dimension
                                        newCell[dim]=0
                                else:
                                    #Subtract
                                    #if this is not the first cell in the given dimension
                                    if(newCell[dim]>0):
                                        newCell[dim]-=1
                                    #if it is the first cell in the given dimension
                                    else:
                                        #move to the end of that dimension
                                        newCell[dim]=cells[dim]-1
                                
                                #Equivalent atom number in 1st unit cell
                                relAtomNum=angleList[angleNum][2+change]%nAtomsPerCell
                            
                        #Change the atom number appropriately
                        angleList[angleNum][2+change]=atomNumList[newCell[0]][newCell[1]][newCell[2]][relAtomNum-1]
                        
                    pairNum+=1
                        
                angleNum+=1

################
## WRITE FILE ##
################e

#File object
outFile=open(outLoc,'w')

#Initial comment
outFile.write('#Muscovite Mica - Heinz INTERFACE Force Field\n') 

#Counters
atomNum=1
atomType=1
molNum=1

#Headers section
print("Writing headers")
outFile.write("{} atoms\n".format(nAtoms))
outFile.write("{} bonds\n".format(nBonds))
outFile.write("{} angles\n".format(nAngles))
outFile.write("\n")
outFile.write("{} atom types\n".format(nAtomTypes))
outFile.write("34 bond types\n")
outFile.write("101 angle types\n")
outFile.write("\n")

#PBC
outFile.write("{} {} xlo xhi\n".format(xlo,xhi))
outFile.write("{} {} ylo yhi\n".format(ylo,yhi))
outFile.write("{} {} zlo zhi\n".format(zlo,zhi))
outFile.write("\n")

#Atom section
print("Writing Atoms section")
outFile.write("Atoms #full\n\n")
#id mol type q x y z cells[0] cells[1] cells[2]

#Loop through all unit cells
for i in range(0,cells[0]):
    for j in range(0,cells[1]):
        for k in range(0,cells[2]):

            count=0
            cellList=[]
            for atom in unitCoords:
                q=charges[count]
                x,y,z=atomCoords[atomNum-1]
                atomType=types[count]
                element=elements.index(names[types[count]-1][0])+1
                
                #Append atomCoords
                atomCoords.append([x,y,z])
                
                #Write file
                outFile.write("{} {} {} {} {} {} {} 0 0 0 # {}\n".format(atomNum,element,atomType,q,x,y,z,names[types[count]-1]))
                
                #Increment atom counters
                atomNum+=1
                count+=1
            
            #Increment molecule counter
            molNum+=1
            


#Bond section
print("Writing Bonds section")
outFile.write("\nBonds\n\n")
for bond in bondList:
    outFile.write(' '.join([str(num) for num in bond]))
    outFile.write("\n")
outFile.write("\n")

#Angles section
print("Writing Angles section")
outFile.write("\nAngles\n\n")
for angle in angleList:
    outFile.write(' '.join([str(num) for num in angle]))
    outFile.write("\n")
outFile.write("\n")

#Paste copied data
outFile.write(copyText)


print("Done! :)")

