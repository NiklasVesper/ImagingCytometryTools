import pandas as pd
import scandir as sd
import os
import numpy as np
import matplotlib.pyplot as plt
import re

from ImagingCytometryToolsIMC.neigboorhood import neigboorhood_cell_type
from ImagingCytometryToolsIMC.neigboorhood import neigboorhood_cell_type_FAK
from ImagingCytometryToolsIMC.neigboorhood import neigboorhood_cell_type_H


folder_dir = r'D:\Test data\20240731_NV_9299_HCC\20240716_NV_9299_HCC__9299_Tumor_1'

for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith("subcell_neighborhood.csv"): # returns all files that end with txt

            filename = os.path.basename(file) # gives you the file name
            filename_string = str(filename)
            filedir_string = str(filedir)[:-len(filename)] # gives you the file directory

            neighborhood_analysis = pd.read_csv(filedir) #opens the found directory

            cell_types = [] # empty list for all cell types
            for index, cell in neighborhood_analysis.iterrows():
                cell_type = []


                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:#
                    cell_type.append('Immune cell')
                    if cell['MeanIntensity_CD3_Nucleus'] > 0.2:
                        if cell['MeanIntensity_CD4_Nucleus'] > 0.2:
                            cell_type.append('CD4 T cell')
                        if cell['MeanIntensity_CD8_Nucleus'] > 1:
                            cell_type.append('CD8 T cell')
                    if cell['MeanIntensity_CD20_Nucleus'] > 0.3:
                        cell_type.append('B cell')
                    if cell['MeanIntensity_CD15_Nucleus'] > 1:
                        cell_type.append('Granulocyte')
                    if cell['MeanIntensity_CD68_Nucleus'] > 0.5:
                        cell_type.append('Myeloid cell')
                    if cell['MeanIntensity_CD11c_Nucleus'] > 0.3:
                        cell_type.append('DC cell')

                elif cell['MeanIntensity_CD45_Nucleus'] < 0.2:
                    cell_type.append('Tissue cell')

                cell_types.append(cell_type)
            neighborhood_analysis['Cell_types'] = cell_types

            immune_type = []
            for index, cell in neighborhood_analysis.iterrows():

                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:
                    if cell['MeanIntensity_CD3_Nucleus'] > 0.2:
                        if cell['MeanIntensity_CD8_Nucleus'] > 1:
                            immune_type.append('CD8 T cell')
                        else:
                            immune_type.append('no')
                    else:
                        immune_type.append('no')
                else:
                    immune_type.append('no')

            neighborhood_analysis['immune_type'] = immune_type

            neighborhood_analysis['CD8 neigboors'] = neigboorhood_cell_type_H(['Immune cell', 'CD8 T cell'],'CD8 T cell',neighborhood_analysis)

            dir = filedir_string + filename_string[:-4] + "_CD8_neighborhood" + ".csv"

            print(dir)
            neighborhood_analysis.to_csv(dir)

            #for CD20 positive cells
            '''
                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:
                    if cell['MeanIntensity_CD20_Nucleus'] > 0.3:
                        immune_type.append('B cell')
                    else:
                        immune_type.append('no')
                else:
                    immune_type.append('no')

            neighborhood_analysis['immune_type'] = immune_type

            neighborhood_analysis['CD20 neigboors'] = neigboorhood_cell_type_H(['Immune cell', 'B cell'], 'B cell',neighborhood_analysis)

            dir = filedir_string + filename_string[:-4] + "_CD20_neighborhood" + ".csv"
            '''

            # for CD11c positive cells
            '''
                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:
                    if cell['MeanIntensity_CD11c_Nucleus'] > 0.3:
                        immune_type.append('DC cell')
                    else:
                        immune_type.append('no')
                else:
                    immune_type.append('no')

            neighborhood_analysis['immune_type'] = immune_type

            neighborhood_analysis['CD11c neigboors'] = neigboorhood_cell_type_H(['Immune cell', 'DC cell'],'DC cell',neighborhood_analysis)

            dir = filedir_string + filename_string[:-4] + "_CD11c_neighborhood" + ".csv"
            '''

        '''

        if cell['MeanIntensity_CD45_Nucleus'] >= 0.3:
            cell_type.append('Immune cell')
            if cell['MeanIntensity_CD3_Nucleus'] > 0.3:
                if cell['MeanIntensity_CD4_Nucleus'] > 1:
                    cell_type.append('CD4 T cell')
                if cell['MeanIntensity_CD8_Nucleus'] > 1.4:
                    cell_type.append('CD8 T cell')
                if cell['MeanIntensity_CD20_Nucleus'] > 0.3:
                    cell_type.append('B cell')
            elif cell['MeanIntensity_CD3_Nucleus'] < 0.3:
                if cell['MeanIntensity_CD20_Nucleus'] > 0.3:
                    cell_type.append('B cell')
                if cell['MeanIntensity_CD15_Nucleus'] > 1.1:
                    cell_type.append('Granulocyte')
                if cell['MeanIntensity_CD68_Nucleus'] > 1.4:
                    cell_type.append('Myeloid cell')
        elif cell['MeanIntensity_CD45_Nucleus'] < 0.3:
            cell_type.append('Tissue cell')
            if cell['MeanIntensity_Ecadherin_Cell'] > 0.8:
                cell_type.append('Hepatocyte')
            if cell['MeanIntensity_Collagen_Cell'] > 5:
                cell_type.append('Fibroblast')



        if cell['MeanIntensity_CD45_Nucleus'] >= 0.3:
            if cell['MeanIntensity_CD3_Nucleus'] > 0.3:


                if cell['MeanIntensity_PD1_Cell'] > 0.18 and cell['MeanIntensity_TCF1_Nucleus'] > 0.18: # tcf1 Pd1 = th 0.18
                    TCF1_PD1.append('TCF-1 pos PD1 pos')

                if cell['MeanIntensity_PD1_Cell'] > 0.18 and cell['MeanIntensity_TCF1_Nucleus'] < 0.18:
                    TCF1_PD1.append('TCF-1 neg PD1 pos')

                if cell['MeanIntensity_PD1_Cell'] < 0.18 and cell['MeanIntensity_TCF1_Nucleus'] > 0.18:
                    TCF1_PD1.append('TCF-1 pos PD1 neg')

                if cell['MeanIntensity_PD1_Cell'] < 0.18 and cell['MeanIntensity_TCF1_Nucleus'] < 0.18:
                    TCF1_PD1.append('TCF-1 neg PD1 neg')



    file['CD8 TCF1+ PD1+ neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD8 T cell'],'TCF-1 pos PD1 pos',file)
    file['CD8 TCF1+ PD1- neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD8 T cell'], 'TCF-1 pos PD1 neg',file)
    file['CD8 TCF1- PD1+ neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD8 T cell'], 'TCF-1 neg PD1 pos',file)
    file['CD8 TCF1- PD1- neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD8 T cell'], 'TCF-1 neg PD1 neg',file)

    file['CD4 TCF1+ PD1+ neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD4 T cell'], 'TCF-1 pos PD1 pos',file)
    file['CD4 TCF1+ PD1- neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD4 T cell'], 'TCF-1 pos PD1 neg',file)
    file['CD4 TCF1- PD1+ neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD4 T cell'], 'TCF-1 neg PD1 pos',file)
    file['CD4 TCF1- PD1- neigboors'] = neigboorhood_cell_type(['Immune cell', 'CD4 T cell'], 'TCF-1 neg PD1 neg',file)




        if cell['MeanIntensity_CD4_Nucleus'] > 0.2:
            cell_type.append('CD4 T cell')
        if cell['MeanIntensity_CD8_Nucleus'] > 0.2:
            cell_type.append('CD8 T cell')
        if cell['MeanIntensity_CD20_Nucleus'] > 0.2:
            cell_type.append('B cell')
        if cell['MeanIntensity_CD15_Nucleus'] > 0.2:
            cell_type.append('Granulocyte')
        if cell['MeanIntensity_CD68_Nucleus'] > 0.4:
            cell_type.append('Myeloid cell')

        if cell['MeanIntensity_CD4_Nucleus'] < 0.2:
            if cell['MeanIntensity_CD8_Nucleus'] < 0.2:
                if cell['MeanIntensity_CD20_Nucleus'] < 0.2:
                    if cell['MeanIntensity_CD15_Nucleus'] < 0.2:
                        if cell['MeanIntensity_CD68_Nucleus'] < 0.4:
                            cell_type.append('Tumor')


        if cell['MeanIntensity_CD4_Nucleus'] < 0.2:
            if cell['MeanIntensity_CD8_Nucleus'] < 0.2:
                if cell['MeanIntensity_CD20_Nucleus'] < 0.2:
                    if cell['MeanIntensity_CD15_Nucleus'] < 0.2:
                        if cell['MeanIntensity_CD68_Nucleus'] < 0.4:

                            if cell['MeanIntensity_FAK_Cytoplasm'] > 0.5:
                                FAK.append('FAK cyto high')
                            else:
                                FAK.append('FAK cyto low')

                        else:
                            FAK.append('no')
                    else:
                        FAK.append('no')
                else:
                    FAK.append('no')
            else:
                FAK.append('no')
        else:
            FAK.append('no')


    file['FAK high'] = neigboorhood_cell_type_FAK(['Tumor'], 'FAK cyto high', file)
    file['FAK low'] = neigboorhood_cell_type_FAK(['Tumor'], 'FAK cyto low', file)
'''