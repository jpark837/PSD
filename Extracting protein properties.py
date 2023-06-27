import pandas as pd
from openpyxl import load_workbook
import numpy as np
import os
import json

# Modules specific to me:
import residueproperties as rp
import SNRCalculation as sc
import analysis as ana
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def to_excel_writer(i, filepath, counter):
    df = pd.DataFrame(np.asarray([i]))
    if counter == 0:
        writer = pd.ExcelWriter(filepath)
        df.to_excel(writer, sheet_name = 'test', index=False, header = False)
    else:
        reader = pd.read_excel(filepath)
        writer = pd.ExcelWriter(filepath, mode = "a", engine = "openpyxl", if_sheet_exists = "overlay")
        df.to_excel(writer, sheet_name = 'test', index=False, header = False, startrow = len(reader)+1)  
    writer.save()

def extract_strucfeat(residue, path_to_folder, path_to_spectrum, path_to_predicted, path_to_NAPS, precursormz, model_num):
    ### Get list of position and ss from residues
    A_SS = rp.SS(path_to_folder, residue, model_num)
    pos = A_SS.get_res_pos() # without parantheses calls pointer instead
    ss = A_SS.get_res_ss()
    
    ### Get relative SASA values from residues
    A_SAS = rp.SAS(path_to_folder, residue, model_num)
    SAS = A_SAS.get_res_SAS()
    
    ### Get hydrogen bonds from residues
    A_hb = rp.hyd_bonds(path_to_folder, residue, model_num)
    hbcount_side = A_hb.get_count_hbonds('sidechain', 2.5,  3.2, 180)
    hbcount_back = A_hb.get_count_hbonds('backbone', 2.5,  3.2, 180)

    ### Get C-terminal residue adjacent from aspartic acid 'D' residues
    A_adj = rp.adjacent(path_to_folder, residue, model_num)
    Ctermadj, Ntermadj = A_adj.get_adj()

    ### Get length of y- and b-ion
    A_length = rp.fragment_length(path_to_folder, residue, model_num)
    rel_y_ion, rel_b_ion = A_length.get_rel_position()

    ### Get salt bridge count from residue
    A_sb = rp.saltbridge(path_to_folder, residue, model_num)
    saltbridge = A_sb.get_count_saltbridge('all', 0, 4, 180)

    ### Get centrality measurements from residues
    A_cent = rp.centrality(path_to_folder, path_to_NAPS, residue, model_num)
    centrality= A_cent.get_centralityspecs()

    ### Get SNR calculation for each...
    if os.path.exists(path_to_spectrum):          
        SNR_PP = sc.fragmentSNR(path_to_spectrum, path_to_predicted, precursormz, residue)
        snr_b, snr_y = SNR_PP.get_SNRI()
        
    else: 
        snr_b, snr_y = [], []
    
    return pos, ss, SAS, hbcount_side, hbcount_back, snr_b, snr_y, Ctermadj, Ntermadj, rel_y_ion, rel_b_ion, saltbridge, centrality
        

def main():
    df = pd.read_excel(r'C:\Users\Jihyun.Park\Documents\Python Scripts\Proteins_to_run_RA.xlsx')
    array = df.to_numpy()

    counter = 0
    sslist, SASlist, hbcount_backlist, hbcount_sidelist, snr_blist, snr_ylist, Ctermadjlist, Ntermadjlist = [], [], [], [], [], [], [], []
    rel_y_ionlist, rel_b_ionlist, saltbridgelist = [], [], []
    centralityarr = np.empty([1,8], dtype = int)
    for i in range(0, np.shape(array)[0]):
        
        path_to_folder = r'C:\Users\Jihyun.Park\Documents\\Alphafold predictions\\Protein List\\output\\' + array[i,0] + r'_output\\'+ array[i,0]
        path_to_spectrum = r'C:\Users\Jihyun.Park\Documents\Python Scripts\Spectras for protein list\Compiled List\\' + array[i,0] + r'.txt'
        path_to_predicted = r'C:\Users\Jihyun.Park\Documents\Python Scripts\Predicted fragments\\DENP\\' + array[i,0] + r'.csv'
        path_to_NAPS = r'C:\Users\Jihyun.Park\Documents\Python Scripts\NAPS\Weighted\\' + array[i,0]
        precursormz = array[i,2]
        residue = 'D'
        model_num = array[i,4]
        A_SS = rp.SS(path_to_folder, residue, model_num)
        pos = A_SS.get_res_pos() # check if there is matching residue
        if len(pos) != 0:

            ### concatenate results:
            if os.path.exists(path_to_spectrum):  

                pos, ss, SAS, hbcount_side, hbcount_back, snr_b, snr_y, Ctermadj, Ntermadj, rel_y_ion, rel_b_ion, saltbridge, centrality = extract_strucfeat(residue, path_to_folder, path_to_spectrum, path_to_predicted, path_to_NAPS, precursormz, model_num)
                sslist += ss
                SASlist += SAS
                hbcount_sidelist += hbcount_side
                hbcount_backlist += hbcount_back
                snr_blist += snr_b
                snr_ylist += snr_y
                Ctermadjlist += Ctermadj
                Ntermadjlist += Ntermadj
                rel_y_ionlist += rel_y_ion
                rel_b_ionlist += rel_b_ion
                saltbridgelist += saltbridge
                centralityarr = np.append(centralityarr, centrality, axis = 0)

                if len(snr_b) != len(rel_b_ion):
                    print('yes')
            
    
    centralityarr = centralityarr[1:,:] # remove first value...
    # plot and analyze data
    p = ana.analysis(snr_blist, snr_ylist, sslist, SASlist, Ctermadjlist, Ntermadjlist, hbcount_backlist, hbcount_sidelist, rel_y_ionlist, rel_b_ionlist, saltbridgelist, centralityarr)
    p.MLR()

    p = ana.plotting(snr_blist, snr_ylist, rel_y_ionlist, rel_b_ionlist)
    snr_max, rel_ionlength = p.snr_handle()
    p.plot_columnscatter(saltbridgelist, 'salt bridge count')
    p.plot_columnscatter(sslist, 'secondary structure')
    p.plot_columnscatter(Ctermadjlist, 'C terminal adjacent residue')
    p.plot_columnscatter(Ntermadjlist, 'N terminal adjacent residue')

        
        


if __name__ == "__main__":
    main()
