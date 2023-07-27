from pymol import cmd
from pymol import stored
import numpy as np
import pandas as pd
import freesasa
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

class SS:
    '''
    Reports secondary structure
    Usage:
    object = rp.SS('path to folder as string', '1 letter residue as string')
    '''
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num
    
    def set_filepaths(self):
        path_to_pdb = self.path_to_folder + r'\\relaxed_model_'+ str(self.model_num) + '_pred_0'

        orig_prot = path_to_pdb+ r'.pdb'
        protonated_prot = self.path_to_folder+ r'\\proteinprotonated.pdb'

        return orig_prot, protonated_prot

    def protonate(self):
        orig_prot, protonated_prot = self.set_filepaths()
        '''
        Protonate the structure in the pdb file and save using pymol
        '''
        cmd.reinitialize()
        cmd.load(orig_prot)
        cmd.h_add(selection = "(all)" )
        cmd.remove(selection = 'solvent')     

        cmd.save(protonated_prot,selection="(all)", state=int(-1) , format='pdb')
        cmd.remove(selection = "(all)")
    
    def find_residue(self):
        '''
        Locate instances of specific amino acid residue
        '''
        self.protonate()
        orig_prot, protonated_prot = self.set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)
        
        aaseq = cmd.get_fastastr(selection="(all)")
        self.list_aa = ''.join(aaseq.split()[1:]) # [1:] removes the header

        index_match = []
        for i in range(0,len(self.list_aa)):
            if self.list_aa[i] == self.residue and self.residue in ['D','E','N'] and i < len(self.list_aa) - 1: # Don't care about aspartic/glutamic acids at the end or beginning
                index_match.append(i)
            elif self.list_aa[i] == self.residue and self.residue == 'P' and i > 0:
                index_match.append(i)
        
        return index_match
    
    def identify_secondary(self):
        '''
        Identify secondary structure of every residue in selected protein. Returns list of secondary structure for each residue. uses DSSP for 8 structure types.
        '''

        # DSSP:
        #Kabsch W & Sander C (1983) DSSP: definition of secondary structure of proteins given a set of 3D coordinates. Biopolymers 22: 2577–2637
        
        orig_prot, protonated_prot = self.set_filepaths()
        p = PDBParser()
        structure = p.get_structure('protein',protonated_prot)
        model = structure[0]
        dssp = DSSP(model, protonated_prot)
        ss = []
        res_pos = []

        for i in range(0, len(self.list_aa)):
            ss.append( dssp[('A',(' ', i+1, ' '))][2] ) # Key set up
            res_pos.append( i+1 )

        return res_pos, ss
    
    def get_res_pos(self):
        '''
        Returns list of indexed position of matched residue
        '''
        index_match = self.find_residue()
        res_pos, ss = self.identify_secondary()
        
        resd_pos = [res_pos[i] for i in index_match ]
        return resd_pos

    def get_res_ss(self):
        '''
        Returns corresponding list of secondary structure of matched residue
        '''
        index_match = self.find_residue()
        res_pos, ss = self.identify_secondary()
        
        res_ss = [ss[i] for i in index_match]
        return res_ss

class SAS:
    '''
    Reports SAS
    Protein and its strcutural properties
    Usage:
    object = rp.SAS('path to folder as string', '1 letter residue as string')
    '''
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num


    def rel_SAS(self):
        '''
        Output: list of relative SAS values for each residue. Cannot work for proteins with multiple chains.
        ''' 

        # Reuse from previous class
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()

        prot_structure = freesasa.Structure(orig_prot)

        result = freesasa.calc(prot_structure, freesasa.Parameters({'algorithm' : freesasa.LeeRichards, 'n-slices': 960}))
        area_classes = freesasa.classifyResults(result, prot_structure)
        rAs = result.residueAreas()
        
        ### Creates a less confusing where dictionary key = chain, value = list of relative SAS values ### 
        key_list = [k for k, v in rAs.items()] # because there are 2 keys in this dictionary 
        rel_SAS = []

        for key_1 in key_list[0]:
            for key_2, value in rAs[key_1].items():
                rel_SAS.append(value.relativeTotal)
        res_pos =  range(1, len(rel_SAS)+1)
        return res_pos, rel_SAS
    
    def get_res_pos(self):
        '''
        Returns list of indexed position of matched residue
        '''
        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()
        res_pos, rel_SAS = self.rel_SAS()
        
        resd_pos = [res_pos[i] for i in index_match ]
        return resd_pos
    
    def get_res_SAS(self):
        '''
        Returns corresponding list of secondary structure of matched residue
        '''
        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()
        res_pos, rel_SAS = self.rel_SAS()
        
        res_SAS = [rel_SAS[i] for i in index_match]
        return res_SAS

class adjacent:
    '''
    Reports the C-terminal side residue adjacent to the residue of interest
    Usage:
    object = rp.adjacent('path to folder as string', '1 letter residue as string')
    '''
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num
    
    def adj(self):
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)

        aaseq = cmd.get_fastastr(selection="(all)")
        list_aa = ''.join(aaseq.split()[1:])

        N_termadj, C_termadj = [], []
        for i in range(0, len(list_aa) - 1): # Because last residue of C-terminal side won't have an adjacent residue
            C_termadj.append(list_aa[i+1])
        for i in range(1, len(list_aa)): # Because first residue of N-terminal side won't have an adjacent residue
            N_termadj.append(list_aa[i-1])
        
        return C_termadj, N_termadj
    
    def get_adj(self):
        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()
        C_termadj, N_termadj = self.adj()

        res_Ctermadj, res_Ntermadj = [],[]
        for i in index_match:
                res_Ctermadj.append(C_termadj[i])
                res_Ntermadj.append(N_termadj[i-1])
                
        return res_Ctermadj, res_Ntermadj


class hyd_bonds:
    '''
    Reports hbonds
    Protein and its structural properties
    Usage:
    object = rp.hyd_bonds('path to folder as string', '1 letter residue as string')
    '''
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num
    
    def int_hbonds(self, selection, min_cutoff, max_cutoff, angle):
        '''
        returns a list of tuples.
        Tuple format: (Atom1, Atom2, Residue 1, Residue 2, Distance between pair)'''
        # Reuse from previous class
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)

        hyd_bonds = cmd.find_pairs(selection + "& symbol n+o", 'all' + "& symbol n+o", mode=1, cutoff=float(max_cutoff), angle=float(angle))
        pairs = []

        for i in hyd_bonds:

            if (i[1], i[0]) in hyd_bonds: #unique only
                hyd_bonds.remove((i[1], i[0])) 
        
            stored.tuple_pair = () ## Tuple of (Atom1, Residue 1, Atom2, Residue 2, Distance between pair)
            
            # write string for selections
            sel_atom1 = i[0][0] + ' and index ' + str(i[0][1])
            sel_atom2 = i[1][0] + ' and index ' + str(i[1][1])### Calculate distance
            
            dist=cmd.distance('temp', sel_atom1, sel_atom2, cutoff=float(max_cutoff), mode = 0)            
            cmd.delete('temp') # To remove pre-drawn labels of bonds

            if dist > min_cutoff:
                ### First atom of the pair                
                cmd.iterate(sel_atom1,
                            'stored.tuple_pair += (name, resi)' )
                            #'print("Element: {} Residue: {} atom index: {}".format(name, resi, index))')
                
                ### Second atom of the pair
                cmd.iterate(sel_atom2,
                            'stored.tuple_pair += (name, resi)' )
                
                ### Record distance
                dist=cmd.distance('internal_hbonds', sel_atom1, sel_atom2, cutoff=float(max_cutoff), mode = 0)
                stored.tuple_pair += (round(dist,2),)
                pairs.append(stored.tuple_pair)
        
        return pairs   
    
    def count_int_hbonds(self, pairs):
        
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)

        ## Get total residue number of protein       
        aaseq = cmd.get_fastastr(selection="(all)")
        list_aa = ''.join(aaseq.split()[1:])
        
        hbonds_count = []
        for i in range(0,len(list_aa)):
            counter = 0
            for j in pairs:                
                if (i+1 == int(j[1])) or (i+1 == int(j[3])):
                    counter += 1
            hbonds_count.append(counter)
        
        return hbonds_count
    
    def get_count_hbonds(self, selection, min_cutoff, max_cutoff, angle):
        '''
        Returns corresponding list of count of intramolecular hydrogen bonds for residue of interest
        '''

        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()

        pairs = self.int_hbonds(selection, min_cutoff, max_cutoff, angle)
        hbonds_count = self.count_int_hbonds(pairs)
        
        res_hbonds_count = [hbonds_count[i] for i in index_match]
        return res_hbonds_count

class fragment_length:
    '''
    Considers the relative position of the fragment 0-1 from n-terminus to c-terminus
    '''
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num
    
    def frag_lengths(self):
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        
        aaseq = cmd.get_fastastr(selection="(all)")
        list_aa = ''.join(aaseq.split()[1:])
    

        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()
        length_y_ion = [] # of residues for each y-ion fragment
        length_b_ion = [] 
        total_prot_length = len(list_aa)

        # y-ion and b-ion
        for i in index_match:
            if self.residue in ['D','E','N'] and i < len(list_aa) - 1: # Because D,E,N cannot be at the C-term end and produce a fragment ion
                length_y_ion.append(i+1)
                length_b_ion.append(total_prot_length - (i+1))
            elif self.residue == 'P' and i > 0:
                length_y_ion.append(i+1)                        # Because P cannot be at the N-term beginning and produce a fragment ion
                length_b_ion.append(total_prot_length - (i+1))
    
        return length_y_ion, length_b_ion, total_prot_length
    
    def get_rel_position(self):
        length_y_ion, length_b_ion, total_prot_length = self.frag_lengths()

        rel_y_ion, rel_b_ion = [],[]
        for i in range(0, len(length_y_ion)):
            
            rel = length_y_ion[i]/total_prot_length
            rel_y_ion.append(rel)
            rel = length_b_ion[i]/total_prot_length
            rel_b_ion.append(rel)

        
        return rel_y_ion, rel_b_ion


class saltbridge:
    def __init__(self, path_to_folder, residue, model_num):
        self.path_to_folder = path_to_folder
        self.residue = residue
        self.model_num = model_num
    
    def int_saltbridge(self, selection, min_cutoff, max_cutoff, angle):    
    # Reuse from previous class
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)
        '''
        Input: selection as string, minimum cutoff distance in angstrom, maximum cutoff distance in angstrom, angle
        Output: Prints and returns a list of tuples.
        Tuple format: (Atom1, Residue 1, Atom2, Residue 2, Distance between pair)
        Only considers LYS, ARG as possible cations and ASP, GLU as possible anions'''

        selection_string = "(" + selection + " & resn Lys & sidechain & n. NZ)" + " | " + "(" + selection + " & resn Arg & sidechain & n. NH1)" + " | " + "(" + selection + " & resn Asp & sidechain & n. OD2)" + " | " + "(" + selection + " & resn Glu & sidechain & n. OE2)" 
        saltbridge_bonds = cmd.find_pairs(selection_string, selection_string, mode=1, cutoff=float(max_cutoff), angle=float(angle))
        pairs = []

        saltbridge_bonds.sort()

        for i in saltbridge_bonds:

            if (i[1], i[0]) in saltbridge_bonds: #unique only
                saltbridge_bonds.remove((i[1], i[0])) 
                
            stored.tuple_pair = () ## Tuple of (Atom 1, Residue number 1, 3 letter residue Atom 1, Atom 2, Residue number 2, 3 letter residue Atom 2, Distance between pair)
            
            # write string for selections
            sel_atom1 = i[0][0] + ' and index ' + str(i[0][1])
            sel_atom2 = i[1][0] + ' and index ' + str(i[1][1])### Calculate distance
            
            dist=cmd.distance('temp', sel_atom1, sel_atom2, cutoff=float(max_cutoff), mode = 0)            
            cmd.delete('temp') # To remove pre-drawn labels of bonds
            
            if dist > min_cutoff:

                ### First atom of the pair                
                cmd.iterate(sel_atom1, 'stored.tuple_pair += (name, resi, resn)' )
                
                ### Second atom of the pair
                cmd.iterate(sel_atom2,'stored.tuple_pair += (name, resi, resn)' )

                # Check if only cationic (Lys,Arg) and anionic (Asp, Glu) pairs:
                if (stored.tuple_pair[2] in ['LYS','ARG'] and stored.tuple_pair[5] in ['ASP','GLU']) or (stored.tuple_pair[5] in ['LYS','ARG'] and stored.tuple_pair[2] in ['ASP','GLU']):

                    ### Record distance
                    dist=cmd.distance('internal_hbonds', sel_atom1, sel_atom2, cutoff=float(max_cutoff), mode = 0)
                    stored.tuple_pair += (round(dist,2),)
                    pairs.append(stored.tuple_pair)
                
        return pairs  

    def count_int_saltbridge(self, pairs):
        
        orig_prot, protonated_prot = SS(self.path_to_folder, self.residue, self.model_num).set_filepaths()
        cmd.reinitialize()
        cmd.load(protonated_prot)

        ## Get total residue number of protein       
        aaseq = cmd.get_fastastr(selection="(all)")
        list_aa = ''.join(aaseq.split()[1:])
        
        saltbridge_count = []
        for i in range(0,len(list_aa)):
            counter = 0
            for j in pairs:                
                if (i+1 == int(j[1])) or (i+1 == int(j[4])): # if matches residue number
                    counter += 1
            saltbridge_count.append(counter)
        
        return saltbridge_count
    
    def get_count_saltbridge(self, selection, min_cutoff, max_cutoff, angle):
        '''
        Returns corresponding list of count of intramolecular hydrogen bonds for residue of interest
        '''

        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()

        pairs = self.int_saltbridge(selection, min_cutoff, max_cutoff, angle)
        saltbridge_count = self.count_int_saltbridge(pairs)
        
        res_saltbridge_count = [saltbridge_count[i] for i in index_match]
        return res_saltbridge_count

class centrality:
    '''
    Reads centrality measurements from NAPS output:
    https://bioinf.iiit.ac.in/NAPS/index.php
    '''
    def __init__(self, path_to_folder, path_to_NAPS, residue, model_num):
        self.path_to_folder = path_to_folder
        self.path_to_NAPS = path_to_NAPS
        self.residue = residue
        self.model_num = model_num
    
    def set_filepaths(self):
        path_to_parameters= self.path_to_NAPS + '_parameters.txt'
        path_to_centrality = self.path_to_NAPS + '_centrality.txt'
        path_to_edgecentrality = self.path_to_NAPS + '_edgecentrality.txt'
        
        return path_to_parameters, path_to_centrality, path_to_edgecentrality

    def read_centrality(self):
        path_to_parameters, path_to_centrality, path_to_edgecentrality = self.set_filepaths()
        df1 = pd.read_csv(path_to_centrality, sep = '\t+|\s+\t+|\t+\s+', engine = 'python', index_col = False)
        df2 = pd.read_csv(path_to_parameters, sep = '\t+|:', engine = 'python', header = None, index_col = 0)
        centrality = df1.to_numpy()
        parameters = df2.T.to_numpy()
        # ADJUST centrality

        centrality[:,6] = centrality[:,6]/parameters[0,3]
        return centrality
    
    
    def get_centralityspecs(self):
        
        index_match = SS(self.path_to_folder, self.residue, self.model_num).find_residue()

        centrality_array = self.read_centrality()

        # AA degree, CC, closeness, betweeness, eigenvector centrality, eccentricity, ANDegree
        centrality = centrality_array[np.asarray(index_match),:]

        return centrality[:,1:]
        