#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 15 22:40:41 2017

@author: Aretas
"""
###Accepts serine protease PDB file as an input
import sys
import math

class PDBatom:
    def __init__(self, string):
        self.atom_number = string[4:11].strip()
        self.name = string[12:17].strip()
        self.residue = string[17:21].strip()
        self.chain = string[21:22].strip()
        self.residue_number = string[22:26].strip()
        self.x_cord = float(string[30:38].strip())
        self.y_cord = float(string[38:46].strip())
        self.z_cord = float(string[46:54].strip())
        self.factor = string[60:67].strip()
        self.fullline = string.strip("\n")

#function to measure distance between two atoms and sort by the max lenght of distance
def distance (list1, list2, distance_dict, max_distance):
    for element1 in list1:
        for element2 in list2:
            x_dist = (element1.x_cord - element2.x_cord)**2
            y_dist = (element1.y_cord - element2.y_cord)**2
            z_dist = (element1.z_cord - element2.z_cord)**2
            distance = math.sqrt (x_dist + y_dist + z_dist)
            if distance < max_distance:
                distance_dict[element1, element2] = distance
                #print (distance_dict)
            #stores the distance in the dictionary atom1, atom2 : distance

#a function to find the angle when all three sides are known;
#a is the HS; b is HD; c is DS
def gamma_angle (a, b, c):
    gamma  = math.acos ((a**2 + b**2 - c**2) / (2*a*b))
    return gamma

#this part deals with the angles of between triad residues to determine the triad
#finding the catalytic triad
def find_triad1 (SH_dict, SD_dict, HD_dict):
    for key, value in SH_dict.items():
        serine_atom = key[0]
        his_atom = key[1]
        for key1, value1 in SD_dict.items():
            serine_atom2 = key1[0]
            asp_atom2 = key1[1]
            if serine_atom == serine_atom2:
                for key2, value2 in HD_dict.items():
                    his_atom3 = key2[0]
                    asp_atom3 = key2[1]
                    if his_atom == his_atom3 and asp_atom2 == asp_atom3 and \
                        serine_atom == serine_atom2 and \
                        value1 > value2 and value1 > value:
                        site_angle = gamma_angle (SH_dict[key], HD_dict[key2], SD_dict[key1])
                        if site_angle > 2 and site_angle < 3.1:
                            radians = site_angle
                            global triad_atoms
                            triad_atoms = serine_atom.atom_number + "-" + \
                                his_atom3.atom_number + "-" + asp_atom3.atom_number
                            triad_residues = serine_atom.residue_number + "-" + \
                                his_atom3.residue_number + "-" + asp_atom3.residue_number
                            print ("The triad atoms of {} chain {} are:".format\
                                (sys.argv[1], his_atom3.chain),\
                                triad_atoms, "The triad residue numbers are: {}"\
                                .format(triad_residues))
                            #print (value, value1, value2)
                            #print (radians)

if __name__ == "__main__":

    #Setting dictionaries and variables to use
    ser_OG_list = []  #atom number : xyz coordinates as list [x, y, z]
    asp_OD2_list = []
    his_ND1_list = []

    try:
        with open(sys.argv[1]) as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    atom = PDBatom(line)
                    if atom.residue == "SER" and atom.name == "OG":
                        ser_OG_list.append(atom)
                    if atom.residue == "ASP" and atom.name == "CG":
                        asp_OD2_list.append(atom)
                    if atom.residue == "HIS" and atom.name == "CE1":
                        his_ND1_list.append(atom)
    except IOError:
        print ('Could not open the file!', sys.argv[1])
        exit()

    #set the lowest distance threshold allowed
    max_distance_SerHis = 4
    max_distance_SerAsp = 9
    max_distance_HisAsp = 5
    distance_dictSH = {} #distance between two atoms; two atom numbers : distance
    distance_dictSD = {}
    distance_dictHD = {}

    radians = 0
    triad_atoms = None

    #finds distance between all S, H, D residues within the set maximum distance
    distance (ser_OG_list, his_ND1_list, distance_dictSH, max_distance_SerHis)
    #print ("------")
    distance (ser_OG_list, asp_OD2_list, distance_dictSD, max_distance_SerAsp)
    #print ("------")
    distance (his_ND1_list, asp_OD2_list, distance_dictHD, max_distance_HisAsp)
    #print ("------")

    find_triad1 (distance_dictSH, distance_dictSD, distance_dictHD)

    if triad_atoms == None:
        print ("Failed to detect the catalytic triad for {}!".format(sys.argv[1]))
