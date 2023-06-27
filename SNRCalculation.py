import numpy as np
import scipy.io
import pandas as pd
import fnmatch

class fragmentSNR:
    
    def __init__(self, path_to_spectrum, path_to_predicted, precursormz, residue):
        self.path_to_spectrum = path_to_spectrum
        self.path_to_predicted = path_to_predicted
        self.precursormz = precursormz
        self.residue = residue

    def extract_data(self, path_to_spectrum, path_to_predicted):
        '''
        For extracting data from csv files
        '''
        df1 = pd.read_csv(path_to_spectrum, sep = '\s+|\t+|\s+\t+|\t+\s+', skiprows = 1, engine = 'python', index_col = False)        
        df2 = pd.read_csv(path_to_predicted,sep = '\s+|\t+|\s+\t+|\t+\s+',  engine = 'python', index_col = False)
        
        array_spectrum, array_fragments = df1.to_numpy(), df2.to_numpy() #4th column = b-ion, 5th column = y-ion

        return array_spectrum, array_fragments

    def match_fragment(self, array_fragments, residue):
        '''
        For pruning results, returns index of matched residue for array of predicted fragment ions
        '''

        index_match = []
        for i in range(0, np.shape(array_fragments)[0]):
            if array_fragments[i, 1] == residue:
                index_match.append(i)

        return np.asarray(index_match)
    
    def window_regexmatcher(self, array_spectrum, fragment_mz, window_size):
        '''
        Takes spectrum array, returns index of matching regex within a window size (max signal)
        '''
        fragment_mz = int(float(fragment_mz))
        index_match = np.where(np.logical_and(array_spectrum[:,0] > fragment_mz-window_size, array_spectrum[:,0] < fragment_mz+window_size) )[0] # a tuple, the first value is an array of indexes
        values_index_match = array_spectrum[index_match,1]

        if np.any(values_index_match) == False: # If there is no matching signal, then automatically assume as 0
            idx_max = None 
        else:
            idx_max = index_match[np.argmax(values_index_match)]

        return idx_max
    
    def prune_signal(self, array_spectrum, array_fragments):
        '''
        Takes spectrum array, finds index of matching fragment, within a m/z window of n (returns the max signal)
        and prunes the array, removing all signal matching the fragment ions provided
        '''
        b_ion, y_ion = array_fragments[:,3], array_fragments[:,4]

        index_b, index_y, index_max, index_bwater, index_ywater, index_yammonia = np.empty([1,0], dtype = int), np.empty([1,0], dtype = int), np.empty([1,0], dtype = int), np.empty([1,0], dtype = int), np.empty([1,0], dtype = int), np.empty([1,0], dtype = int) # Initialize empty index array
        index_preammonia, index_prewater = np.empty([1,0], dtype = int), np.empty([1,0], dtype = int)
        ### Prune b-ion fragments
        for fragment_mz in b_ion:    
            x = self.window_regexmatcher(array_spectrum.astype(int), fragment_mz, 5)
            index_b  = np.append(index_b, x)
            
            # Remove water (b-18)
            if x != None:
                y = self.window_regexmatcher(array_spectrum.astype(int), array_spectrum[x, 0]-18, 2)
                index_bwater = np.append(index_bwater, y)
        
        ### Prune y-ion fragments
        for fragment_mz in y_ion:    
            x = self.window_regexmatcher(array_spectrum.astype(int), fragment_mz, 5)
            index_y = np.append(index_y, x)
            
            # Remove water and ammonia(y-17, y-18)
            if x != None:
                y = self.window_regexmatcher(array_spectrum.astype(int), array_spectrum[x, 0]-17, 2)
                index_yammonia = np.append(index_yammonia, y)
                y = self.window_regexmatcher(array_spectrum.astype(int), array_spectrum[x, 0]-18, 2)
                index_ywater = np.append(index_ywater, y)

        ### Prune precursor
        x = self.window_regexmatcher(array_spectrum.astype(int), self.precursormz, 5)

        # Remove water and ammonia(y-17, y-18)
        if x != None:
            y = self.window_regexmatcher(array_spectrum.astype(int), int(self.precursormz)-17, 2)
            index_preammonia = np.append(index_preammonia, y)
            y = self.window_regexmatcher(array_spectrum.astype(int), int(self.precursormz)-18, 2)
            index_prewater = np.append(index_prewater, y)
        
        ## Add the indexes for fragment ions (y + b) and precursor
        index_max = np.append(index_b, index_y)
        index_deleteall = np.append(index_max, x)
        index_delete = index_deleteall[index_deleteall!=np.array(None)].astype(int)

    
        for i in [index_bwater, index_yammonia, index_ywater, index_preammonia, index_prewater]:
            index_delete = np.append(index_delete, i[i!=np.array(None)].astype(int))
        
        index_delete = np.unique(index_delete) # prune duplicates

        pruned_array_spectrum = np.delete(array_spectrum, index_delete, axis = 0) # pruned array, removing the signal matching the fragment

        return index_b, index_y, index_max, pruned_array_spectrum
    
    def get_SNRI(self):
        array_spectrum, array_fragments = self.extract_data(self.path_to_spectrum, self.path_to_predicted)
        index_b, index_y, index_max, pruned_array_spectrum = self.prune_signal(array_spectrum, array_fragments)

        ### calculate standard deviation
        sd = pruned_array_spectrum[:,1].std(axis=0, ddof=0)

        m_b, m_y = [], []
        for i in range(0, np.shape(index_b)[0]):
            if index_b[i] == None: # If there is no signal of the predicted fragment ion, then automatically assign as 0
                m_b.append(0)
            else:
                m_b.append(array_spectrum[index_b[i],1])
            
            if index_y[i] == None:
                m_y.append(0)
            else:
                m_y.append(array_spectrum[index_y[i],1])
        
        snr_b, snr_y = np.where(sd == 0, 0, np.divide(m_b,sd)), np.where(sd == 0, 0, np.divide(m_y,sd)) # removes issue of 0/0

        idx_matchres = self.match_fragment(array_fragments, self.residue)

        if idx_matchres.size != 0:
            snr_b, snr_y = snr_b[idx_matchres], snr_y[idx_matchres]
            snr_b, snr_y = np.around(snr_b, decimals = 2).tolist(), np.around(snr_y, decimals = 2).tolist()
        else:
            snr_b, snr_y = [], []

        return snr_b, snr_y