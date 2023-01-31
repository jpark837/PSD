import numpy as np
import scipy.io
import pandas as pd
import fnmatch

# https://github.com/hrtlacek/SNR/blob/main/SNR.ipynb
# https://github.com/hrtlacek/SNRconda install -c anaconda scipy

#df = pd.read_excel(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\TestMSMSSpectra.xlsx')
#To read from text/csv file:
df = pd.read_csv(r'C:\Users\Jihyun.Park\OneDrive - USDA\Documents\Python Scripts\Raw ASCII files\230126_YS_CKF_D15_PM7819_RM6607W.txt', sep = '\s+|\t+|\s+\t+|\t+\s+', skiprows = 1, engine = 'python', index_col = False)
array = df.to_numpy()
print(df)
print(array)

#print(type(array[0,0]))
def regexmatcher(array, regex_pattern):
    '''
    Takes array, finds regex pattern and  returns list of matched index
    '''
    index_match = []
    for i in range(0,len(array)):
        if fnmatch.fnmatch(array[i,0],regex_pattern):
            index_match.append(i)
    return index_match

def windowfinder(array, window_size):


def signaltonoise(array, axis=0, ddof=0):

    str_b = array.astype(str)
    
    index_match = regexmatcher(str_b, '5142.*')
    b = np.delete(array, index_match, axis = 0) # pruned array

    # Convert back to int
    b = b.astype(int)
    #calculate standard deviation
    # sd = b[:,1].std(axis=axis, ddof=ddof)
    sd = array[:,1].std(axis=axis, ddof=ddof)

    index_match = regexmatcher(str_b, '5142.*')
    m = array[index_match[0],1]
    print('this is m')
    print(m)
    print('this is sd')
    print(sd)
    return np.where(sd == 0, 0, np.divide(m,sd))

snr = signaltonoise(array)
print('this is snr')
print(snr)