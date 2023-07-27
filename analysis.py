import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels
import statsmodels.formula.api as smf

class plotting:
    def __init__(self, snr_blist, snr_ylist, rel_y_ionlist, rel_b_ionlist):
        self.snr_blist = snr_blist
        self.snr_ylist = snr_ylist
        self.rel_y_ionlist = rel_y_ionlist
        self.rel_b_ionlist = rel_b_ionlist
    
    def snr_handle(self):
        '''
        For finding max snr between b and y fragments
        '''
        snr_max = []
        rel_ionlength = []

        for i in range(0,len(self.snr_blist)):
            if self.snr_blist[i] > self.snr_ylist[i]:
                snr_max.append(self.snr_blist[i])
                rel_ionlength.append(self.rel_b_ionlist[i])
            else:
                snr_max.append(self.snr_ylist[i])
                rel_ionlength.append(self.rel_y_ionlist[i])
        return snr_max, rel_ionlength
    
    def plot_scatter(self, x, property):
        '''
        Matplotlib scatter plot of x vs snr_max
        '''
        snr_max, rel_ionlength = self.snr_handle()
        plt.scatter(np.asarray(x), np.asarray(snr_max))
        plt.xlabel(property)
        plt.ylabel('max SNR')
        plt.show()

    def plot_columnscatter(self,x, property):
        '''
        seaborn based column scatter
        '''
        snr_max, rel_ionlength = self.snr_handle()
        d = {property: x, 'max SNR': snr_max}
        df = pd.DataFrame(data = d)

        sns.set_style("whitegrid")
        sns.boxplot(data = df, x = property, y = 'max SNR')
        plt.show()

class analysis:
    def __init__(self, snr_blist, snr_ylist, sslist, SASlist, Ctermadjlist, Ntermadjlist, hbcount_backlist, hbcount_sidelist, rel_y_ionlist, rel_b_ionlist, saltbridgelist, centralityarr):
        self.snr_blist = snr_blist
        self.snr_ylist = snr_ylist
        self.sslist = sslist
        self.SASlist = SASlist
        self.Ctermadjlist = Ctermadjlist
        self.Ntermadjlist = Ntermadjlist
        
        #degree, CC, closeness, betweeness, eigenvector centrality, eccentricity, ANDegree
        self.centralityarr = centralityarr
        self.hbcount_backlist = hbcount_backlist
        self.hbcount_sidelist = hbcount_sidelist
        self.rel_y_ionlist = rel_y_ionlist
        self.rel_b_ionlist = rel_b_ionlist
        self.saltbridgelist = saltbridgelist
    
    def MLR(self):
        self.snr_max, self.rel_ionlength = plotting(self.snr_blist, self.snr_ylist, self.rel_b_ionlist, self.rel_y_ionlist).snr_handle()
        self.centralityarr = self.centralityarr.astype('float64')
        d = {'snr':self.snr_max, 'backbone': np.asarray(self.hbcount_backlist), 'sidechain': np.asarray(self.hbcount_sidelist), 
             'secondarystructure': np.asarray(self.sslist), 'SAS':np.asarray(self.SASlist), 'Ctermadj': np.asarray(self.Ctermadjlist), 'Ntermadj': np.asarray(self.Ntermadjlist),
             'relionlength': np.asarray(self.rel_ionlength), 'saltbridgecount': np.asarray(self.saltbridgelist),
             'degree':self.centralityarr[:,0], 'CC':self.centralityarr[:,1], 'closeness': self.centralityarr[:,2], 'betweeness':self.centralityarr[:,3],
             'eigenvec_centrality': self.centralityarr[:,4], 'eccentricity': self.centralityarr[:,5],
             'AND':self.centralityarr[:,6], 'strength':self.centralityarr[:,7]}

        df = pd.DataFrame(data=d)
        df.to_excel(excel_writer = r"C:\Users\Jihyun.Park\Documents\Python Scripts\df_all.xlsx")

