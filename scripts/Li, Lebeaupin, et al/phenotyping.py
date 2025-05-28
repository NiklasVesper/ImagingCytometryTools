import pandas as pd
import scandir as sd
import os

import warnings
import re

'''
The "phenotyping" script assigns cells their type and after that assigns each cell that could be phenotyped they neighborhood.
The script takes files that were previously generated with the "neighborhood analysis" script.

Continue now with the "neighborhood visualization" script.
'''

warnings.simplefilter(action='ignore', category=FutureWarning)

#assigns cells their neighboring cell types
def neigboorhood_cell_type(Cell_type,Subtype,file):

    Neigboorhood_fin= [] #a list for the whole neighborhood phenotypes
    for index, cell in file.iterrows(): #loops through all cells

        if cell['Cell_types'] == Cell_type and cell['cell_type_neighborhood'] == Subtype: #checks through conditions of the desired cell types

            Neighborhood = re.split(r',', cell['Neighborhood'].replace('[', '').replace(']', '')) #splits the sting found in the neighborhood column

            Neighborhood_list = [] #list for the neighboring cell numbers
            for x in Neighborhood: #goes through the neighborhood of a cell
                try:
                    Neighborhood_list.append(int(float(x))) #appends the neighboring cell numbers into a list
                except ValueError: #goes on if there is nothing in the neighborhood
                    pass

            Neighborhood_cell = [] #a list for single cell neighborhood phenotypes
            for index, cell in file.iterrows(): #loops through all cells

                if cell['Cell_number'] in Neighborhood_list: #checks if the number of a cell is in the neighborhood list
                    Neighborhood_cell.append(cell['Cell_types']) #appends the phenotype into a new list

            Neigboorhood_fin.append(Neighborhood_cell) #appends all surrounding phenotypes into a new list

        else:
            Neigboorhood_fin.append('currently not of interest')

    return(Neigboorhood_fin)

folder_dir = r'' #folder directory

#Example for neighborhood analysis of CD8 T cells. Any other cell lineage or subset can be analyzed the same way.
for paths, dirs, files in sd.walk(folder_dir): #goes through all files and folders in given directory

    for file in os.listdir(paths): #goes through all files in a folder
        filedir = os.path.join(paths, file) #returns full file directory

        if filedir.endswith("subcell_neighborhood.csv"): #checks if the file has 'subcell_neighborhood.csv' in its name

            filename = os.path.basename(file) #gives you the file name
            filename_string = str(filename) #turns the filename into a string
            filedir_string = str(filedir)[:-len(filename)] #gives you the file directory as a string

            neighborhood_analysis = pd.read_csv(filedir) #opens the found directory

            cell_types = [] #empty list for all cell types
            for index, cell in neighborhood_analysis.iterrows():
                cell_type = [] #empty list for signular cell type

                '''
                The cell types are assigned based on thresholds that were determined by gating on the images in CellProfiler.
                '''

                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:
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
            neighborhood_analysis['Cell_types'] = cell_types #adds a new column with the cell types to the pandas data frame

            cell_type_neighborhood = []
            for index, cell in neighborhood_analysis.iterrows():

                '''
                Defines the cell type that one is interested in.
                '''

                if cell['MeanIntensity_CD45_Nucleus'] >= 0.2:
                    if cell['MeanIntensity_CD3_Nucleus'] > 0.2:
                        if cell['MeanIntensity_CD8_Nucleus'] > 1:
                            cell_type_neighborhood.append('CD8 T cell')
                        else:
                            cell_type_neighborhood.append('no')
                    else:
                        cell_type_neighborhood.append('no')
                else:
                    cell_type_neighborhood.append('no')

            neighborhood_analysis['cell_type_neighborhood'] = cell_type_neighborhood #adds a new column for the cell type of choice to the pandas data frame

            neighborhood_analysis['CD8 neighbors'] = neigboorhood_cell_type(['Immune cell', 'CD8 T cell'],'CD8 T cell',neighborhood_analysis) #assigns cell type of choice their neighboring cell types

            dir = filedir_string + filename_string[:-4] + "_CD8_neighborhood" + ".csv" #create a new file name and directory

            print(dir)
            neighborhood_analysis.to_csv(dir) #exports the file
