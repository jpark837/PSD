import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class plotting:
    def __init__(self, snr_b, snr_y):
        self.snr_b = snr_b
        self.snr_y = snr_y
    
    def snr_handle(self):
        '''
        For finding max snr between b and y fragments
        '''
        snr_max = np.maximum(np.asarray(self.snr_b), np.asarray(self.snr_y))
        return snr_max
    
    def plot_scatter(self, x, property):
        '''
        Matplotlib scatter plot of x vs snr_max
        '''
        snr_max = self.snr_handle()
        plt.scatter(np.asarray(x), np.asarray(snr_max))
        plt.xlabel(property)
        plt.ylabel('max SNR')
        plt.show()

    def plot_columnscatter(self,x, property):
        '''
        seaborn based column scatter
        '''
        snr_max = self.snr_handle()
        d = {property: x, 'max SNR': snr_max}
        df = pd.DataFrame(data = d)
        sns.set_style("whitegrid")
        sns.boxplot(data = df, x = property, y = 'max SNR')
        plt.show()
