import pandas as pd
from scipy.stats import mannwhitneyu

'''
The "neighborhood testing" script performes the statistical analysis to the "neighborhood visualization" script.
The script takes the .csv files that were previously generated with the "neighborhood visualization" script.
'''

#folder directories
folder_dir_high = r'F:\ATF6 FAK\All\High\CD8\CD8_immune_neighborhood.csv'
folder_dir_low = r'F:\ATF6 FAK\All\Low\CD8\CD8_immune_neighborhood.csv'

#data import
df_high = pd.read_csv(folder_dir_high)
df_low = pd.read_csv(folder_dir_low)

types = ['CD4','CD8','CD20','CD15','CD68','CD11c'] #cell phenotypes for analysis

for type in types: #goes through phenotypes

    print(type)

    #performs Mann Whitney u test
    U1, p = mannwhitneyu(df_high[type].to_list(), df_low[type].to_list())
    nx, ny = len(df_high[type].to_list()), len(df_low[type].to_list())
    U2 = nx*ny - U1

    #gives you the Bonferroni corrected p value
    print(p)
    if p < (0.05/18):
        print('found one')

    print('------------')
