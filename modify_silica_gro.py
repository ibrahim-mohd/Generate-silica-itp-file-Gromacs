# Written By Mohd Ibrahim
# email: ibrahim.mohd@biophys.mpg.de

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda
import argparse


def create_universe (n_atoms, name, resname, positions):

    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate


    u_new.add_TopologyAttr('name',    [name]*n_atoms)
    u_new.add_TopologyAttr('resname', [resname]*n_atoms)
    u_new.atoms.positions = positions
    
    return u_new


def get_index (dict_si):

    index = ''

    for  i in list (dict_si.keys()):

        if not dict_si [i]:  index = index + ' ' + i

    return index



parser = argparse.ArgumentParser(description='Creates itp file for silica')
parser.add_argument('-f', dest='gro_file', type=str, default='conf.gro',help='input gro file')
parser.add_argument('-o', dest='ouput_file', type=str, default='out.gro',help='output gro file')
parser.add_argument('-hf', dest='prot_percentage', type=float, default=90, help='percentage of surface oxygen protonated')

## Input file ##################33
args        = parser.parse_args()
input_file  = args.gro_file
output_file = args.ouput_file
prot_percentage    = args.prot_percentage
bond_lengths       = 1.65
u            = mda.Universe (input_file)

all_si      =   u.select_atoms ("name si")

z_min_si    =  min (all_si.positions [:,2]) + 1

z_max_si    =  max (all_si.positions [:,2]) - 1
###########

crystal           =  u.select_atoms ("all")


## Find the top oxygen connected to each surface silicon atom  for top surface
top_surface_si       =  u.select_atoms ("name si and prop z>%f"%(z_max_si))
bottom_surface_si    =  u.select_atoms ("name si and prop z< %f"%(z_min_si))

pairs, dists      =  mda.lib.distances.capped_distance(top_surface_si.positions, crystal.positions, min_cutoff=1e-4,
                                                   max_cutoff = bond_lengths, box=None)

dict_si_top = {}

for key in top_surface_si.indices: dict_si_top [str(key)] = []

for i, j in pairs:
    
    index_si = top_surface_si.indices [i]
    
    if crystal.positions [j, 2] > top_surface_si.positions [i, 2]:  dict_si_top [str (index_si)].append (j)
    
## Find the two top oxygen connected to each surface silicon atom  for bottom surface

pairs, dists         =  mda.lib.distances.capped_distance(bottom_surface_si.positions, crystal.positions, min_cutoff=1e-4,
                                                   max_cutoff = bond_lengths, box=None)

dict_si_bottom = {}

for key in bottom_surface_si.indices: dict_si_bottom [str(key)] = []

for i, j in pairs:
    
    index_si = bottom_surface_si.indices [i]
    
    if crystal.positions [j, 2] < bottom_surface_si.positions [i, 2]:  dict_si_bottom [str (index_si)].append (j)
    
    


############# Add oxygen to make both surface symmetirc 

# Top oxygen

index = get_index (dict_si_top)
A     = u.select_atoms ("index %s"% (index))
positions = A.positions
positions [:,2] += 1.6
top_ox_add = create_universe (n_atoms=len (A),  name="osib", resname="osib", positions=positions)


## Bottom oxygen
index = get_index (dict_si_bottom)
B     = u.select_atoms ("index %s"% (index))
positions =B.positions
positions [:,2] -= 1.6
bottom_ox_add = create_universe (n_atoms=len (B),  name="osib", resname="osib", positions=positions)



#################  Create Universe with added oxygens ##################################33

u_all = mda.Merge (crystal, bottom_ox_add.select_atoms ("all"), top_ox_add.select_atoms ("all"))


###########################################################
all_si      =   u_all.select_atoms ("name si")
z_min_si =  min (all_si.positions [:,2]) + 1

z_max_si =  max (all_si.positions [:,2]) - 1
###########

crystal           =  u_all.select_atoms ("all")
## Find the top oxygen connected to each surface silicon atom  for top surface
top_surface_si       =  u_all.select_atoms ("name si and prop z>%f"%(z_max_si))

pairs, dists      =  mda.lib.distances.capped_distance(top_surface_si.positions, crystal.positions, min_cutoff=1e-4,
                                                   max_cutoff = bond_lengths, box=None)

dict_si_top = {}

for key in top_surface_si.indices: dict_si_top [str(key)] = 0

for i, j in pairs:
    
    index_si = top_surface_si.indices [i]
    
    if crystal.positions [j, 2] > top_surface_si.positions [i, 2]:
        if dict_si_top [str (index_si)]:
            
            j_existing = dict_si_top [str (index_si)]
            
            if  crystal.positions [j, 2] >  crystal.positions [j_existing, 2]: dict_si_top [str (index_si)] = j
            
        else: dict_si_top [str (index_si)] = j
       
       
        #dict_si_top [str (index_si)].append (j)
    
## Find the  top oxygen connected to each surface silicon atom  for bottom surface
bottom_surface_si    =  u_all.select_atoms ("name si and prop z< %f"%(z_min_si))

pairs, dists         =  mda.lib.distances.capped_distance(bottom_surface_si.positions, crystal.positions, min_cutoff=1e-4,
                                                   max_cutoff = bond_lengths, box=None)

dict_si_bottom = {}

for key in bottom_surface_si.indices: dict_si_bottom [str(key)] = 0

for i, j in pairs:
    
    index_si = bottom_surface_si.indices [i]
    
    if crystal.positions [j, 2] < bottom_surface_si.positions [i, 2]:
        
        if dict_si_bottom [str (index_si)]:
            
            j_existing = dict_si_bottom [str (index_si)]
            
            if  crystal.positions [j, 2] <  crystal.positions [j_existing, 2]: dict_si_bottom [str (index_si)] = j
            
        else: dict_si_bottom [str (index_si)] = j
        
        
############################## Choose random surface oxygens to protonate ###############################################

hydrogen_bond_length = 1.2

prot_percentage /= 100

n_hydrogen =  int (prot_percentage *len (dict_si_top))


Top_groups_si_prot     =  np.random.choice (list (dict_si_top.keys ()), n_hydrogen, replace=False) # keys to protonate
Bottom_groups_si_prot  =  np.random.choice (list (dict_si_bottom.keys ()), n_hydrogen,replace=False) # keys to pro

# Top hydrogens
oxygen_indices_top =  ""
silicon_indices_top = ""

for index in Top_groups_si_prot: 

    oxygen_indices_top =  oxygen_indices_top + ' ' + str (dict_si_top [index])
    
    silicon_indices_top =  silicon_indices_top + ' ' + index

top_ox      = u_all.select_atoms ("index %s" %(oxygen_indices_top))

hydrogen_pos = top_ox.positions
hydrogen_pos [:,2] += hydrogen_bond_length
name = "hso"
top_hydrogen = create_universe (n_atoms=n_hydrogen, name = name, resname=name, positions=hydrogen_pos)



# bottom hydrogens
oxygen_indices_bottom =  "" # protonated oxygen
silicon_indices_bottom   = "" # silicon connected to protonated oxygen
for index in Bottom_groups_si_prot: 

    oxygen_indices_bottom  =  oxygen_indices_bottom + ' ' + str (dict_si_bottom [index])
    silicon_indices_bottom =  silicon_indices_bottom + ' ' + index 

bottom_ox      = u_all.select_atoms ("index %s" %(oxygen_indices_bottom))

hydrogen_pos = bottom_ox.positions
hydrogen_pos [:,2] -= hydrogen_bond_length
name = "hso"
bottom_hydrogen  = create_universe (n_atoms=n_hydrogen, name = name, resname=name, positions=hydrogen_pos)

################# Find silicon and oxygens not protonated ##########################3

############### Top Surface  #####################

no_prot_oxygen_indices_top  =  "" # protonated oxygen
no_prot_silicon_indices_top =  "" # silicon connected to protonated oxygen

for index in list (dict_si_top.keys ()):
    
    if index not in Top_groups_si_prot:
        no_prot_oxygen_indices_top  =  no_prot_oxygen_indices_top + ' ' + str (dict_si_top [index])
        no_prot_silicon_indices_top =  no_prot_silicon_indices_top + ' ' + index 

############### Bottom Surface  #####################
no_prot_oxygen_indices_bottom =  ""
no_prot_silicon_indices_bottom  =  ""

for index in list (dict_si_bottom.keys ()):
    
    if index not in Bottom_groups_si_prot:
        no_prot_oxygen_indices_bottom  =  no_prot_oxygen_indices_bottom + ' ' + str (dict_si_bottom [index])
        no_prot_silicon_indices_bottom =  no_prot_silicon_indices_bottom + ' ' + index    
    
#################### Create universe for each of above cases #################################

protonated_si    =  u_all.select_atoms ("index %s or index %s" %(silicon_indices_bottom,
                                                                 silicon_indices_top))

protonated_ox    = u_all.select_atoms ("index %s or index %s" % (oxygen_indices_bottom,
                                                                 oxygen_indices_top) )

################ Non protonated #################################3

no_protonated_si =  u_all.select_atoms ("index %s or index %s" %(no_prot_silicon_indices_bottom,
                                                             no_prot_silicon_indices_top) )

no_protonated_ox = u_all.select_atoms ("index %s or index %s" % (no_prot_oxygen_indices_bottom,
                                                                 no_prot_oxygen_indices_top))


### Rest of the group: i.e excluding above surface atoms
surface_index = ''

for i,j in zip ( list (dict_si_top.keys ()), list (dict_si_top.values ())):
    surface_index = surface_index + ' ' + i + ' ' +str (j)  + ' '

for i,j in zip ( list (dict_si_bottom.keys ()), list (dict_si_bottom.values ())):
    surface_index = surface_index + ' ' + i + ' ' +str (j) + ' '
    
rest = u_all.select_atoms ("all and not index %s" % (surface_index))


#################### Create universe for each of above cases #################################
# 1: Si connected to protonated oxygen (sis): silinol silicon 

name = "sis" # silinol si

n_atoms = len (protonated_si)
positions = protonated_si.positions

protonated_si = create_universe (n_atoms=n_atoms, name=name, resname=name, positions=positions)


# 2: Si connected to non-protonated oxygen: Siloxide group

name = "sio" # silinol si

n_atoms = len (no_protonated_si)
positions = no_protonated_si.positions

no_protonated_si = create_universe (n_atoms=n_atoms, name=name, resname=name,
                                positions=positions)

# 3: Protonated oxygens: oh

name = "ohs" # silinol si

n_atoms = len (protonated_ox)
positions = protonated_ox.positions

protonated_ox = create_universe (n_atoms=n_atoms, name=name, resname=name, positions=positions)


# 4: Non-protonated oxygens: oh

name = "osi" # silinol si

n_atoms = len (no_protonated_ox)
positions = no_protonated_ox.positions

no_protonated_ox = create_universe (n_atoms=n_atoms, name=name, resname=name, positions=positions)




a = mda.Merge (rest,
               no_protonated_si.select_atoms ("all"),
               protonated_si.select_atoms ("all"),
               protonated_ox.select_atoms ("all"),
               no_protonated_ox.select_atoms ("all"),
               top_hydrogen.select_atoms ("all"),
               bottom_hydrogen.select_atoms ("all"))
               
a.dimensions =  u.dimensions

a.select_atoms ("all").write (output_file)
