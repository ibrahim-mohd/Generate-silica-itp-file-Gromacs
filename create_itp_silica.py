# Written By Mohd Ibrahim
# email: ibrahim.mohd@biophys.mpg.de

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda
import argparse


parser = argparse.ArgumentParser(description='Creates itp file for silica')
parser.add_argument('-f', dest='gro_file', type=str, default='unitcell.gro',help='unitcell gro file')
parser.add_argument('-o', dest='itp_file', type=str, default='silica.itp',help='output itp file')
parser.add_argument('-d', dest='defaults', type=bool, default=True, help='if True creates as standalone itp file')

## Input file ##################33
args        = parser.parse_args()
input_file  = args.gro_file
output_file = args.itp_file
Defaults    = args.defaults

u            = mda.Universe (input_file)
crystal      = u.select_atoms ("all")
bond_lengths = 1.65 # longest bond, the hydrogen involving bonds are shorter


#  The parameters are taken from the paper:
#  Force Field and a Surface Model Database for Silica to Simulate Interfacial Properties in Atomic Resolution  ####
#  (Fateme S. Emami, Valeria Puddu, Rajiv J. Berry, Vikas Varshney, Siddharth V. Patwardhan, Carole C. Perry, and Hendrik Heinz)
#  Chemistry of Materials Journal ACS



# Atomnames are as:
# si--- Silicon: osib----> oxygen connected to silicon in bulk: osis----> oxygen silicon in silinol: hsi----> silinol hydrogen


############################################################################################################################################
## Forcefield parameters:

''' Bond and Angle values'''
'''Bond'''
# eq bond distances
si_o_r0   = round ( 1.65/10, 5) 
o_h_r0   = round (0.945/10, 5)

# spring constants
si_o_keq   =  round (285*(2*4.184*100), 2) # convert to gromacs uints, terms in bracket are conversion factor to Gromacs units
o_h_keq    =  round (495*(2*4.184*100), 2)

# Angle info
'''Angles'''
o_si_o_theta0  = 109.5
si_o_si_theta0 = 149.0
si_o_h_theta0  = 115.0

# spring constants
o_si_o_keq  = round ( 100*(2*4.184), 2)
si_o_si_keq = round (100*(2*4.184), 2)
si_o_h_keq  = round( 50*(2*4.184),2)

##########################################################
## the dictionar keys are "s" intead of "si" just for convienece later

bond_info = dict (s_o =  dict (keq = si_o_keq, r0 = si_o_r0),
                  o_s =  dict (keq = si_o_keq, r0 = si_o_r0),
                  o_h =  dict (keq = o_h_keq, r0 = o_h_r0),
                  h_o =  dict (keq = o_h_keq, r0 = o_h_r0))

angle_info  =  dict (o_s_o = dict (keq = o_si_o_keq, theta0  =  o_si_o_theta0), 
                     s_o_s = dict (keq = si_o_si_keq, theta0 =  si_o_si_theta0),
                     s_o_h = dict (keq = si_o_h_keq, theta0  =  si_o_h_theta0),
                     h_o_s = dict (keq = si_o_h_keq, theta0  =  si_o_h_theta0)
                      )


lj_params   = dict ( si  = dict (sigma = 4.15/10, epsilon=0.093*4.184, at_no=14),
                     sis = dict (sigma = 4.15/10, epsilon=0.093*4.184, at_no=14),
                     sio = dict (sigma = 4.15/10, epsilon=0.093*4.184, at_no=14),
                     osib= dict (sigma = 3.47/10, epsilon=0.054*4.184, at_no=8),
                     osi = dict (sigma = 3.47/10, epsilon=0.122*4.184, at_no=8),
                     ohs = dict (sigma = 3.47/10, epsilon=0.122*4.184, at_no=8),
                     hso = dict (sigma = 1.085/10,epsilon=0.015*4.184, at_no=1))


silica_charges = dict (si   = 1.10,  sis= 1.10, sio=0.725,
                       osib = -0.55, ohs = -0.675, osi = -0.9,
                       hso  = 0.4)


silica_masses  = dict (si = 28.08600,sis=28.08600, sio=28.08600,
                       osib = 15.99940, ohs = 15.99940, osi = 15.99940,
                       hso = 1.008)




## Find all the pairs within the 'bond_length' cutoff

pairs, dists  =  mda.lib.distances.capped_distance(crystal.positions, crystal.positions, min_cutoff=1e-4,
                                                   max_cutoff = bond_lengths, box=u.dimensions)
## List the bonds according to above pairs
# gromacs indices starts from 1 instead of 0 

Bonds = []

for p in pairs:
    
    b1 = (p[0], p[1])
    b2 = (p[1], p[0])
    
    if crystal.positions [p[0],2] - crystal.positions [p[1],2] > bond_lengths: continue # For periodic molecule in xy-dimension

    if b1 not in Bonds and b2 not in Bonds: Bonds.append (b1)

## Get the Angles based on the bond information above

Angles = []

for index, b in enumerate (Bonds):
    
    for a in Bonds [index+1:]:
        
        if b[0] == a [0]:
            angle = (min ([b[1], a[1]]),b[0], max ([b[1], a[1]]))
            Angles.append (angle)
        
        if b[0] == a [1]:
            angle = (min ([b[1], a[0]]),b[0], max ([b[1], a[0]]))
            Angles.append (angle)
        
        
        if b[1] == a [0]:
            angle = (min ([b[0], a[1]]),b[1], max ([b[0], a[1]]))
            Angles.append (angle)
        
        if b[1] == a [1]:
            angle = (min ([b[0], a[0]]),b[1], max ([b[0], a[0]]))
            Angles.append (angle)
            
    
    
    
## Write an itp file based on above information

itp_file  = open(output_file, "w+")
Defaults = False # Just for now, can use as an input later, set to True if want to create a standalone silica itp
if Defaults:
    
    itp_file.write ("[ defaults ]  \n;; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
    itp_file.write ("1               2               yes             0.5     0.83\n")


## Atom types 
itp_file.write ("[ atomtypes ] \n;name   atomic no.      mass     chBrge   ptype   sigma         epsilon\n")

for atom_name in list (lj_params.keys ()):
    
    sigma, epsilon  = lj_params [atom_name]["sigma"],  lj_params [atom_name]["epsilon"], 
    mass, atomic_no = silica_masses [atom_name], lj_params [atom_name]["at_no"]
    
    itp_file.write ("%6s %8d %12.5f  %8.5f %5s %10.5f %12.5f \n" % (atom_name, atomic_no, mass, 0.000, "A",
                                                                    sigma, epsilon))

    
#################Molecualr type
itp_file.write ("\n[ moleculetype ] \n; name \t rexcl \n silica \t 3\n")
itp_file.write ("\n[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   ;bond_type\n")


## Write [ atoms ]

qtot = 0

for name, index in zip (crystal.names, crystal.indices):
    
    nr, type_, resid,res,atom  = index+1, name,index+1, name, name
    cgnr, charge, mass        = index+1, silica_charges [name], silica_masses[name]
    
    qtot += charge
    
    itp_file.write ("%5d %5s %5d %5s %5s %5d %10.4f %10.5f \t ;qtot = %6.3f\n" % (nr, type_, resid,res, atom,cgnr, charge, mass, qtot ) )
    
    
## Write [bonds]    
itp_file.write ("\n[ bonds ] \n; ;   ai     aj  funct   r             k\n")
#itp_file.write ("[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   ;bond_type\n")

for index, (i,j) in enumerate (Bonds) :
    
    atom1, atom2 = crystal.names [i], crystal.names [j]
    
    key = atom1 [0] + "_" + atom2[0]
    #print (bond_info [key])
    r0, keq = bond_info [key]["r0"],  bond_info [key]["keq"]
    
    itp_file.write ("%8d %6d %5d %8.5f  %8.5f \t ;%5s-%5s\n" % (i+1, j+1,1, r0, keq, atom1, atom2) )
    
# Write [ Angles ]
itp_file.write ("\n[ angles ] \n;   ai     aj     ak    funct   theta         cth\n")
 
for index, (i,j,k) in enumerate (Angles) :
    
    atom1, atom2, atom3 = crystal.names [i], crystal.names [j], crystal.names [k]
    
    key = atom1 [0] + "_" + atom2[0] + "_" +  atom3[0] 
    #print (bond_info [key])
    theta0, keq = angle_info [key]["theta0"],  angle_info [key]["keq"]
    
    itp_file.write ("%6d %6d %6d %6d %12.5f %12.5f \t ;%5s-%5s-%5s\n" % (i+1, j+1, k+1,1, theta0, keq, atom1, atom2, atom3) )
    

itp_file.close ()

