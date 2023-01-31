#%%

# SHEBANG for python interpreter
#!"C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\environments\conda_prot_properties\python"

#### IMPORT PACKAGES ####
# from Bio.PDB import *
from string import ascii_uppercase
from collections import defaultdict
import sys
import fnmatch
import os
# print(sys.path)
# print(sys.executable)
from pymol import cmd
from pymol import stored
import ipywidgets
import freesasa
import itertools
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.mol as mol
import numpy as np
import re
import shutil
import pandas as pd

#### FUNCTIONS ####

def create_SAS_dict(prot_structure):
    '''Input: protein structure from freesasa.
       Output: dictionary of relative SAS values where dictionary key = chain,value = list of rel. SAS values
    ''' 
    
    # ### How to access one residue SAS values ###
    # residue_SAS = rAs["A"]["2"]
    # df = pd.DataFrame({})
    # print(residue_SAS.total)
    # print(residue_SAS.sideChain)
    # print(residue_SAS.relativeTotal)
    # print(residue_SAS.relativeSideChain)

    # yahO.bonds = struc.connect_via_residue_names(yahO) # Based on residues that appear in Chemical Component Dictionary
    # resd1 = yahO[:, (yahO.res_id == 1)]
    # print("Bonds (atoms names):")
    # print(resd1.atom_name[resd1.bonds.as_array()[:, :2]])

    result = freesasa.calc(prot_structure, freesasa.Parameters({'algorithm' : freesasa.LeeRichards, 'n-slices': 960}))
    area_classes = freesasa.classifyResults(result, prot_structure)
    rAs = result.residueAreas()
    
    ### Creates a less confusing where dictionary key = chain, value = list of relative SAS values ### 
    key_list = [k for k, v in rAs.items()] # because there are 2 keys in this dictionary
    rel_SAS_dict = {}  

    for key_1 in key_list[0]:
        rel_SAS = []
        for key_2, value in rAs[key_1].items():
            rel_SAS.append(value.relativeTotal)
        rel_SAS_dict[key_1] = rel_SAS
    
    return rel_SAS_dict

def create_bonds_dict(prot_structure):
    prot_structure.bonds = struc.connect_via_distances(prot_structure[0]) # connect_via_residues was not giving correct bonds
    
    bonds_dict = defaultdict(dict)
    chains = list(ascii_uppercase)
    tot_residues = struc.get_residue_count(prot_structure)

    # Create a bonds dictionary with 2 keys: Key 1 : Chain ID, Key 2: Residue ID
    prot_length = np.shape(prot_structure.res_id)[0]
    for i in range(0,prot_length):
        bonds_dict[prot_structure.chain_id[i]][prot_structure.res_id[i]] = []

    # Create a bonds dictionary and assign the values to the keys: Value = all bonds of a residue
    key_list = [k for k, v in bonds_dict.items()]
    #print(key_list)
        
    for key_1 in key_list:
        for key_2, value in bonds_dict[key_1].items():
            residue = prot_structure[:, (prot_structure.res_id == key_2)]
            bonds_dict[key_1][key_2] = residue.atom_name[residue.bonds.as_array()[:, :2]]   
    return key_list, bonds_dict

def count_labile_protons(key_list, bonds_dict, path_to_pdb):    
    ### Using Biotite to report atom bonds###
    # https://www.biotite-python.org/index.html
 
    ### Count labile protons in the group
    counter = 0
    for key_1 in key_list: # For group
        for key_2, value in bonds_dict[key_1].items():
            for i in range(0,np.shape(value)[0]):
                if (fnmatch.fnmatch(value[i,0], 'N*') or fnmatch.fnmatch(value[i,0], 'O*')) and fnmatch.fnmatch(value[i,1], 'H*'):
                    counter += 1
                elif fnmatch.fnmatch(value[i,0], 'H*') and (fnmatch.fnmatch(value[i,1], 'N*') or fnmatch.fnmatch(value[i,1], 'O*')):  
                    counter += 1

    print("this is the counter before accounting for charged AA: {}". format(counter))
    cmd.load(path_to_pdb)
    aaseq = cmd.get_fastastr(selection="(all)")
    pattern = re.compile("")

    ### Amino acids that is + charge (has extra proton)
    count_K, count_H, count_R = aaseq.count('K'), aaseq.count('H'), aaseq.count('R')

    ### Amino acids that is - charge (missing a proton)
    count_E, count_D = aaseq.count('E'), aaseq.count('D')
   
    # 
    counter = counter - count_K - count_H - count_R + count_E + count_D
    print("this is the counter after accounting for charged AA: {}". format(counter))

    # "Close" the file from pymol
    cmd.remove(selection = "(all)")
    return counter
    
def list_hb(selection, selection2=None, cutoff=3.2,
        angle=45, mode=1, hb_list_name='hbonds',print_distances=1,write_distances_file=None):
    
    # based on a script found at pastebin.com/5mkwWJdd,
    # which is in turn based on my script list_hbonds.py
    # Copyright (c) 2010 Robert L. Campbell

  """
  USAGE

  list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                     [angle (default=45)], [mode (default=1)],
                     [hb_list_name (default='hbonds')], [print_distances (default=1)]
                     [write_distances (default=None)]

  The script automatically adds a requirement that atoms in the
  selection (and selection2 if used) must be either of the elements N or
  O.

  If mode is set to 0 instead of the default value 1, then no angle
  cutoff is used, otherwise the angle cutoff is used and defaults to 45
  degrees.

  e.g.
  To get a list of all H-bonds within chain A of an object and to save it to a file
    list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds,write_distances_file=abc-hbonds.dat

  To get a list of H-bonds between chain B and everything else:
    list_hb 1tl9 & c. b, 1tl9 &! c. b
  """
  if write_distances_file:
    hb_data = open(write_distances_file, 'w')

  cutoff = float(cutoff)
  angle = float(angle)
  mode = int(mode)
  print_distances=int(print_distances)

  # ensure only N and O atoms are in the selection
  selection = selection + " & e. n+o"
  if not selection2:
    hb = cmd.find_pairs(selection , selection, mode=mode, cutoff=cutoff, angle=angle)
  else:
    selection2 = selection2 + " & e. n+o"
    hb = cmd.find_pairs(selection, selection2, mode=mode, cutoff=cutoff, angle=angle)

# convert hb list to set to remove duplicates
  hb_set = set()
  for atoms in hb:
    a = [atoms[0],atoms[1]]
    a.sort()
    hb_set.add(tuple(a))

# convert set back to list and sort for easier reading
  hb = list(hb_set)
  hb.sort(key=lambda x: x[0][1])

  stored.listA = []
  stored.listB = []
  stored.listC = []

  for pairs in hb:
    cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]),
            'stored.listA.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]),
            'stored.listB.append( "%1s/%3s`%s/%s/%i " % (chain, resn, resi, name, ID),)')

    stored.listC.append(cmd.distance(hb_list_name, "%s and index %s" %
                (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0],
                    pairs[1][1])))

  for line in enumerate(stored.listA):
    # if print_distances:
    #   print("%s   %s   %.2f" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
    if write_distances_file:
      hb_data.write("%s   %s   %.2f\n" % (stored.listA[line[0]], stored.listB[line[0]], stored.listC[line[0]]))
  if write_distances_file:
    hb_data.close()

  return stored.listA, stored.listB, stored.listC

def prune_same_residue(A,B,C):

  index_remove = []
  for i in range(0,len(A)):
    res1, res2 = A[i], B[i]
    if fnmatch.fnmatch(res1[:7],res2[:7]):
      index_remove.append(i)
    
  A = np.delete(A, index_remove).tolist()
  B = np.delete(B, index_remove).tolist()
  C = np.delete(C, index_remove).tolist()

  return A, B, C



def to_excel_indv(rel_SAS_dict, counter,A,B,C):
    # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_excel.html#:~:text=To%20write%20to%20multiple%20sheets,necessary%20to%20save%20the%20changes.
    table_setup = [A,B,C]
    table_setup = np.asarray(table_setup).T
    
    writer = pd.ExcelWriter(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\test.xlsx')
    df1 = pd.DataFrame(table_setup, columns = ['atom 1', 'atom 2', 'distance in angstrom'])
    df1.to_excel(writer, sheet_name = 'internal hydrogen bonds', index=False)
    print(counter)
    df2 = pd.DataFrame([counter], columns = ['labile proton count'])
    df2.to_excel(writer, sheet_name = 'labile proton count', index=False)

    for key1 in [*rel_SAS_dict]:
      rel_SAS_list = [*rel_SAS_dict[key1]]
    table_setup =  np.asarray([range(1,len(rel_SAS_list)+1), rel_SAS_list]).T

    df3 = pd.DataFrame(table_setup, columns = ['residue number', 'relative SAS values'])
    df3.to_excel(writer, sheet_name = 'relative SAS', index=False)
    
    writer.save()

def calculate_metrics (rel_SAS_dict, counter, A,B,C):
    table_setup = np.asarray([A,B,C]).T
    num_hyd_bonds = np.shape(table_setup)[0]

    ### Total relative SASA (SUM(relative SASA / number of residues))
    for key1 in [*rel_SAS_dict]:
      rel_SAS_list = [*rel_SAS_dict[key1]]
    
    total_rel_SASA = np.sum(np.asarray(rel_SAS_list))

    num_residues = len(rel_SAS_list)

    return num_hyd_bonds, counter, total_rel_SASA, num_residues

def to_excel_metrics(num_hyd_bonds_list, total_rel_SASA_list, labileprotons_list, num_residues_tot_list, array):
    # set up results table
    table_setup = np.asarray([array[:,0].tolist(), array[:,1].tolist(), num_residues_tot_list, num_hyd_bonds_list, labileprotons_list, total_rel_SASA_list]).T
    print(np.shape(table_setup))
    writer = pd.ExcelWriter(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Results_compiled.xlsx')
    df1 = pd.DataFrame(table_setup, columns = ['identifier', 'description', 'number of residues', 'Number of intramolecular hydrogen bonds', 'Number of labile protons', 'total relative SASA'])
    df1.to_excel(writer, sheet_name = 'results_compiled', index=False)
    
    writer.save()


def main():
    ### Read excel file of proteins to run script on
    df = pd.read_excel(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Proteins_to_run.xlsx')
    array = df.to_numpy()
    num_hyd_bonds_list = []
    total_rel_SASA_list = []
    labileprotons_list = []
    num_residues_tot_list = []   

    for i in array[:,0]: # identifier
    #for i in ['A']:
      print("yes")
      path_to_folder = r'C:\\Users\\Jihyun.Park\\OneDrive - USDA\\Documents\\Alphafold predictions\\Protein List\\output\\' + i + r'_output\\'+ i
      path_to_pdb = path_to_folder + r'\\relaxed_model_1_pred_0'

      orig_prot = path_to_pdb+ r'.pdb'
      protonated_prot = path_to_folder+ r'\\proteinprotonated.pdb'
      # Using FreeSASAy
      # https://freesasa.github.io/doxygen/CLI.html

      structure = freesasa.Structure(orig_prot)
      ### Create dictionary of relative SAS ###
      rel_SAS_dict = create_SAS_dict(structure)

      ### Protonate the structure in the pdb file and save using pymol ###
      cmd.reinitialize()
      cmd.load(orig_prot)
      cmd.h_add(selection = "(all)" )
      cmd.remove(selection = 'solvent')     

      cmd.save(protonated_prot,selection="(all)", state=int(-1) , format='pdb')
      cmd.remove(selection = "(all)")


      ### Create dictionary of atom bond reports ###
      # Path to protonated file

      pdb_file = pdb.PDBFile.read(protonated_prot)
      protein = pdb_file.get_structure()
      key_list, bonds_dict = create_bonds_dict(protein)


      # Count number of labile protons
      counter = count_labile_protons(key_list, bonds_dict, protonated_prot)

      ### Try to print out hydrogen bonds

      cmd.reinitialize()
      cmd.load(protonated_prot)
      #A,B,C = list_hb("(all)", cutoff=3.2, angle=180, write_distances_file=r"C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\test.txt")
      
      # For side chain atoms only:
      A,B,C = list_hb("sidechain", cutoff=3.2, angle=180, write_distances_file=r"C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\test.txt")
      
      # To prune same residue hydrogen bonds
      A,B,C = prune_same_residue(A, B, C)

      # Append the metrics to list
      num_hyd_bonds, counter, total_rel_SASA, num_residues = calculate_metrics(rel_SAS_dict, counter, A,B,C)
      num_hyd_bonds_list.append(num_hyd_bonds)
      total_rel_SASA_list.append(total_rel_SASA)
      labileprotons_list.append(counter)
      num_residues_tot_list.append(num_residues)
      
      # If wanting to run through pymol cmd:
      #cmd.do(r"run C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\list_hb.py")
      #print('second')
      #cmd.do(r"list_hb (all), cutoff=3.2, angle=180, write_distances_file=C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\test.txt")

      ### Write individual results to excel file
      to_excel_indv(rel_SAS_dict, counter,A,B,C)

    
      ### Write customized report of each protein to excel file
    to_excel_metrics(num_hyd_bonds_list, total_rel_SASA_list, labileprotons_list, num_residues_tot_list, array)

if __name__ == "__main__":
    main()


#%%
