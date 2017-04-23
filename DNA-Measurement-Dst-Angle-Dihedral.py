# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 21:02:01 2017

@author: Dimitri Cabaud
This program calculate the pseudo-bond length, bending angles 
and dihedral angles for C-alpha carbons of residues in any PDB file. 

"""
import numpy as np
from numpy.linalg import norm
import math

def import_data_from_pdb(nameofthefile):
    """
    Import the data of C-alpha atoms from a PDB file
    Atom/Atom serial number/Atom Name/Alternate location indicator/Residue Name/Chain Identifier/Residue Sequence number/Code for insertions of residues/x orthogonal angstrom coordinate/y orthogonal angstrom coordinate/z orthogonal angstrom coordinate
    :param string name of the file
    :return a list of lists containing all the informations 
    """
    data = []
    
    with open(nameofthefile, "r") as f:
        for line in f.readlines():
            if '  CA  ' in line and 'ATOM' in line:
                data.append(line.strip().split())
    return (data)
    

def distance(data, index1, index2):
    """
    Calculates the distance between two C-alpha atoms in the backbone of a polypeptide/protein
    :param list of list of C-alpha molecules informations:
           Atom/Atom serial number/Atom Name/Alternate location indicator/Residue Name/Chain Identifier/Residue Sequence number/Code for insertions of residues/x orthogonal angstrom coordinate/y orthogonal angstrom coordinate/z orthogonal angstrom coordinate
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :return distance in angstroms
    """
    coordinates = []
    line1, line2 = 0,0
    for i in range(len(data)):
        if str(data[i][1]) == str(index1):
            line1 = i
        if str(data[i][1]) == str(index2):
            line2 = i

    coordinates.append(data[line1][6:9])
    coordinates.append(data[line2][6:9])
    return np.linalg.norm(coordinates)

def dihedral(data, index1, index2, index3, index4):
    """
    Calculates the dihedral angle of a set of 4 C-alpha atoms in a pdb
    :param list of list of C-alpha molecules informations:
        Atom/Atom serial number/Atom Name/Alternate location indicator/Residue Name/Chain Identifier/Residue Sequence number/Code for insertions of residues/x orthogonal angstrom coordinate/y orthogonal angstrom coordinate/z orthogonal angstrom coordinate
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :return dihedral angle
    """
    
    coordinates=np.zeros([4,3])   
    lines = []
    for i in range(len(data)):
        if str(data[i][1]) == str(index1):
            lines.append(i)
               
        if str(data[i][1]) == str(index2):
            lines.append(i)

        if str(data[i][1]) == str(index3):
            lines.append(i)

        if str(data[i][1]) == str(index4):
            lines.append(i)
            
    coordinates[0] = data[lines[0]][6:9]
    coordinates[1] = data[lines[1]][6:9]
    coordinates[2] = data[lines[2]][6:9]
    coordinates[3] = data[lines[3]][6:9]
      
    matrix_vector=coordinates[1:,:]-coordinates[0:-1,:]      
    
    cosphi=np.dot(np.cross(matrix_vector[0,:],-matrix_vector[1,:]),np.cross(matrix_vector[1,:],-matrix_vector[2,:]))/(np.linalg.norm(np.cross(matrix_vector[0,:],matrix_vector[1,:]))*np.linalg.norm(np.cross(matrix_vector[1,:],matrix_vector[2,:])))
    sinphi=np.dot(matrix_vector[1,:],np.cross(np.cross(matrix_vector[0,:],-matrix_vector[1,:]),np.cross(matrix_vector[1,:],-matrix_vector[2,:])))/(np.linalg.norm(matrix_vector[1,:])*np.linalg.norm(np.cross(matrix_vector[0,:],matrix_vector[1,:]))*np.linalg.norm(np.cross(matrix_vector[1,:],matrix_vector[2,:])))

    if cosphi>0:
        phi=np.arcsin(sinphi)*180/np.pi
    elif sinphi<0:
        phi=-np.arccos(cosphi)*180/np.pi
    else:
        phi=np.arccos(cosphi)*180/np.pi
    return phi
 
 
 
def angle(data, index1, index2, index3): 
    """
    Calculates the angle of a set of 3 C-alpha atoms in a pdb
    :param list of list of C-alpha molecules informations:
        Atom/Atom serial number/Atom Name/Alternate location indicator/Residue Name/Chain Identifier/Residue Sequence number/Code for insertions of residues/x orthogonal angstrom coordinate/y orthogonal angstrom coordinate/z orthogonal angstrom coordinate
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :param integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
    :return angle
    """
    lines = []
    for i in range(len(data)):
        if str(data[i][1]) == str(index1):
            lines.append(i)
               
        if str(data[i][1]) == str(index2):
            lines.append(i)

        if str(data[i][1]) == str(index3):
            lines.append(i)
   
    coordinates1 = np.array([float(data[lines[0]][6]),float(data[lines[0]][7]),float(data[lines[0]][8])])
    coordinates2 = np.array([float(data[lines[1]][6]),float(data[lines[1]][7]),float(data[lines[1]][8])])
    coordinates3 = np.array([float(data[lines[2]][6]),float(data[lines[2]][7]),float(data[lines[2]][8])])
            
    x = coordinates2-coordinates1 
    y = coordinates2-coordinates3 
    Vec1_2 = norm(x)
    Vec2_3 = norm(y)
    Norm1_2 = x / Vec1_2;
    Norm2_3 = y / Vec2_3;
    res = Norm1_2[0] * Norm2_3[0] + Norm1_2[1] * Norm2_3[1] + Norm1_2[2] * Norm2_3[2];
    angle = math.acos(res)*180.0/ math.pi
    return angle 

def interface():
    print("Hello, please type the name of the pdb file you want to use and put it in the same folder of this program")
    print("For example: 1ubq.pdb")
    name_pdb = str(input())
    print('\n')
    data = import_data_from_pdb(name_pdb)
    print("Ok so now what type of measurement do you want to do ?")
    print("Type 1 for distance between two C-alpha atoms")
    print("Type 2 for angle between three C-alpha atoms")
    print("Type 3 for dihedral angle between four C-alpha atoms")
    answer = int(input())
    print('\n')
    if answer == 1:
        print("Please type the index of the first C-alpha atom")
        calpha1 = int(input())
        print('\n')
        print("Please type the index of the second C-alpha atom")
        calpha2 = int(input())
        print('\n')
        dst = distance(data,calpha1,calpha2)
        print('Distance between Alpha Carbon: ',calpha1,'and Alpha Carbon: ',calpha2,' = ',dst,' Angstroms')
    if answer == 2:
        print("Please type the index of the first C-alpha atom")
        calpha1 = int(input())
        print('\n')
        print("Please type the index of the second C-alpha atom")
        calpha2 = int(input())
        print('\n')
        print("Please type the index of the third C-alpha atom")
        calpha3 = int(input())
        print('\n')
        angleCalpha = angle(data,calpha1,calpha2,calpha3)
        print('Angle between Alpha Carbon: ',calpha1,', Alpha Carbon: ',calpha2,' and Alpha Carbon: ',calpha3,' = ',angleCalpha,' °')
    if answer == 3:
        print("Please type the index of the first C-alpha atom")
        calpha1 = int(input())
        print('\n')
        print("Please type the index of the second C-alpha atom")
        calpha2 = int(input())
        print('\n')
        print("Please type the index of the third C-alpha atom")
        calpha3 = int(input())
        print('\n')
        print("Please type the index of the fourth C-alpha atom")
        calpha4 = int(input())
        print('\n')
        dihedralCalpha = dihedral(data,calpha1,calpha2,calpha3,calpha4)
        print('Dihedral angle between Alpha Carbon: ',calpha1,', Alpha Carbon: ',calpha2,', Alpha Carbon: ',calpha3,' and Alpha Carbon: ',calpha4,' = ',dihedralCalpha,' °')


"""
Main
""" 
interface()
    