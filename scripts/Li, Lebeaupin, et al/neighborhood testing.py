import pandas as pd
from scipy.stats import mannwhitneyu

'''
The "neighborhood testing" script performes the statistical analysis to the "neighborhood visualization" script.
The script takes the .csv files that were previously generated with the "neighborhood visualization" script.
'''

#file directories
file_dir_1 = r''
file_dir_2 = r''

#data import
df_1 = pd.read_csv(file_dir_1)
df_2 = pd.read_csv(file_dir_2)

types = ['CD4','CD8','CD20','CD15','CD68','CD11c'] #cell phenotypes for analysis

for type in types: #goes through phenotypes

    print(type)

    #performs Mann Whitney u test
    U1, p = mannwhitneyu(df_1[type].to_list(), df_2[type].to_list())
    nx, ny = len(df_1[type].to_list()), len(df_2[type].to_list())
    U2 = nx*ny - U1

    #gives you the Bonferroni corrected p value
    print(p)
    if p < (0.05/18):
        print('found one')

    print('------------')
