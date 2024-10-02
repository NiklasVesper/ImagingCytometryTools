import pandas as pd
from scipy.stats import mannwhitneyu

folder_dir_high = r'F:\ATF6 FAK\All\High\CD8\CD8_immune_neighborhood.csv'
folder_dir_low = r'F:\ATF6 FAK\All\Low\CD8\CD8_immune_neighborhood.csv'

df_high = pd.read_csv(folder_dir_high)
df_low = pd.read_csv(folder_dir_low)

types = ['CD4','CD8','CD20','CD15','CD68','CD11c']

for type in types:

    print(type)

    U1, p = mannwhitneyu(df_high[type].to_list(), df_low[type].to_list())
    
    nx, ny = len(df_high[type].to_list()), len(df_low[type].to_list())
    U2 = nx*ny - U1

    print(p)

    if p < (0.05/18):
        print('found one')

    print('------------')
