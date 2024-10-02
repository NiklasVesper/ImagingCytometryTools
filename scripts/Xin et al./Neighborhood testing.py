import random
import scandir as sd
import os
import pandas as pd
import re
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

#Filedirlist = []


'''
for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory


    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith("CD8_neighborhood.csv"):


            filename = os.path.basename(file)
            filename_string = str(filename)
            filedir_string = str(filedir)[:-len(filename)]

            if 'Tumor_1' in filename_string:
                Filedirlist.append(filedir)

pdlist = []

for el in Filedirlist:
    pdlist. append(pd.read_csv(el))

Big_data = pd.concat([t for t in pdlist], ignore_index=True)

Big_interaction = []
Total_Cell_count = []

for index, cell in Big_data.iterrows():

    if cell['Cell_types'] == "['Immune cell', 'CD8 T cell']":
        Total_Cell_count.append('+1')

        Neigboorhood = re.split(r'], ', cell['CD8 neigboors'])

        Cell_list = []
        for x in Neigboorhood:
            Cell_list.append(x.replace('[', '').replace(']', ''))
            Big_interaction.append(x.replace('[', '').replace(']', ''))

print(len(Total_Cell_count))
print(Big_interaction)
'''
