import pandas as pd
from openpyxl import load_workbook
import numpy as np
import os


# Modules specific to me:
import residueproperties as rp
import SNRCalculation as sc
import plotting as pltg
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

def extract_strucfeat(residue, path_to_folder, path_to_spectrum, path_to_predicted, path_to_NAPS):
    ### Get list of position and ss from aspartic acid 'D' residues
    A_SS = rp.SS(path_to_folder, residue)
    pos = A_SS.get_res_pos() # without parantheses calls pointer instead
    ss = A_SS.get_res_ss()
    
    ### Get relative SASA values from aspartic acid 'D' residues
    A_SAS = rp.SAS(path_to_folder, residue)
    SAS = A_SAS.get_res_SAS()
    
    ### Get hydrogen bonds from aspartic acid 'D' residues
    A_hb = rp.hyd_bonds(path_to_folder, residue)
    hbcount = A_hb.get_count_hbonds('(all)', 2.5,  3.2, 180)

    ### Get C-terminal residue adjacent from aspartic acid 'D' residues
    A_adj = rp.adjacent(path_to_folder, residue)
    Ctermadj = A_adj.get_adj()

    ### Get centrality measurements from aspartic acid 'D' residues
    A_cent = rp.centrality(path_to_folder, path_to_NAPS, residue)
    degree = A_cent.get_centralityspecs()

    ### Get SNR calculation for each...
    if os.path.exists(path_to_spectrum):          
        SNR_PP = sc.fragmentSNR(path_to_spectrum, path_to_predicted, residue)
        snr_b, snr_y = SNR_PP.get_SNRI()
        
    else: 
        snr_b, snr_y = [], []
    
    return pos, ss, SAS, hbcount, snr_b, snr_y, Ctermadj, degree
        

def main():
    df = pd.read_excel(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Proteins_to_run.xlsx')
    array = df.to_numpy()
    #array = np.asarray([['M','blah']])

    counter = 0
    sslist, SASlist, hbcountlist, snr_blist, snr_ylist, Ctermadjlist, degreelist = [], [], [], [], [], [], []
    for i in range(0, np.shape(array)[0]):
        path_to_folder = r'C:\\Users\\Jihyun.Park\\OneDrive - USDA\\Documents\\Alphafold predictions\\Protein List\\output\\' + array[i,0] + r'_output\\'+ array[i,0]
        path_to_spectrum = r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Spectras for protein list\Compiled List\\' + array[i,0] + r'.txt'
        path_to_predicted = r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Predicted fragments\\' + array[i,0] + r'.csv'
        path_to_NAPS = r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\NAPS\Unweighted\\' + array[i,0]

        residue = 'D'
        pos, ss, SAS, hbcount, snr_b, snr_y, Ctermadj, degree = extract_strucfeat(residue, path_to_folder, path_to_spectrum, path_to_predicted, path_to_NAPS)

        ### concatenate results:
        if os.path.exists(path_to_spectrum):   
            sslist += ss
            SASlist += SAS
            hbcountlist += hbcount
            snr_blist += snr_b
            snr_ylist += snr_y
            Ctermadjlist += Ctermadj
            degreelist += degree
        
        for j in ["Identifier: " + array[i,0] + ". " + array[i,1], 'Residue position of ' + residue ,pos, 
                  'secondary structure of residue position of '  + residue, ss, 'relative SASA of residue position of ' + residue, SAS,
                  'Number of hydrogen bonds of residue position of ' + residue, hbcount, 'C-terminal residue adjacent to residue position of '  + residue, Ctermadj, 
                  'centrality degree of residue position of ' + residue, degree,
                  'SNR score of b-fragment of '  + residue, snr_b, 'SNR score of y-fragment of ' + residue, snr_y, ' ']:
            to_excel_writer(j,r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\prop_res.xlsx', counter)
            counter += 1
    
    # plot and analyze data
    p = pltg.plotting(snr_blist, snr_ylist)
    # p.plot_scatter(SASlist, 'relative SASA')
    # p.plot_scatter(hbcountlist, 'number of hydrogen bonds')
    # p.plot_scatter(degreelist, 'degree')
    p.plot_columnscatter(degreelist, 'degree')
    # p.plot_columnscatter(sslist, 'secondary structure')
    # p.plot_columnscatter(Ctermadjlist, 'C terminal adjacent residue')

        
        


if __name__ == "__main__":
    main()
