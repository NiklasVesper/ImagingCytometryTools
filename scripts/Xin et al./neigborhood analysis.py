import pandas as pd
import scandir as sd
import os

from ImagingCytometryTools.neigboorhood import neigboorhood_cell_type

folder_dir = r'D:\ATF6'# folder directory

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
            
