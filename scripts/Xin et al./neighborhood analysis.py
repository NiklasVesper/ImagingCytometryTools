import pandas as pd
import os
import scandir as sd

import numpy as np
import shapely.geometry
import warnings

'''
The "neighborhood analysis" script examines the spatial proximity of all cells within the tissue.
The script takes files that were previously generated with the "match cells with their nuclei" script.
'''

warnings.simplefilter(action='ignore', category=FutureWarning)

#spatial proximity of all cells within the tissue
def neighborhood(Cells):
    print('I am looking for my friends :)')

    #creates lists of all unique images
    UniqueNames_Full_cell = Cells.ImageNumber.unique()

    #creates a pandas data frame dictionary of all unique images
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  # creates a dictionary of all unique images

    #resets indexes for identified objects in the image blocks
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    Cell_number = []
    Neighborhood = []

    #takes connected images and they objects form from the dictionary
    for key, value in DataFrameDict_Full_cell.items():

        #gets all x and y positions of cells on an image and creates an array of them
        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        count_cell = -1 #counter
        for cell in position_cells: #goes throw all cell positions
            cell_position = shapely.geometry.Point(cell) #creates a point at the center of the cell
            neighborhood_radius = DataFrameDict_Full_cell['AreaShape_MaxFeretDiameter'].mean() #gives you the mean MaxFeretDiameter of all cells
            cell_neighborhood = shapely.geometry.Point(cell).buffer(neighborhood_radius * 1.5) #defines the radius for the spatial proximity
            count_cell = count_cell + 1 #ticks up the counter
            Cell_number.append(count_cell)

            Whats_in_the_hood = [] #list for the individual neighborhood of a cell
            count_other_cell = -1 #counter
            for other_cell in position_cells: #goes throw all other cell positions
                other_cell_position = shapely.geometry.Point(other_cell) #creates a point at the center of the cell
                count_other_cell = count_other_cell + 1 #ticks up the counter

                if other_cell_position == cell_position: #skipps the own cell position
                    continue
                elif other_cell_position.within(cell_neighborhood) == True: #if a cell is within the defined radius it will be considered to be within the neighborhood
                    Whats_in_the_hood.append(count_other_cell)

            Neighborhood.append(Whats_in_the_hood)

    Cells['Neighborhood'] = Neighborhood #creates a new columb for the neighborhood
    Cells['Cell_number'] = Cell_number #creates a new columb for the cell number

    return (Cells)


folder_dir = r''  #folder directory

for paths, dirs, files in sd.walk(folder_dir):  #goes throw all files and folders in given directory

    for file in os.listdir(paths):  #goes throw all files in a folder
        filedir = os.path.join(paths, file)  #returns full file directory

        if filedir.endswith("single_cells_and_organells.csv"):  #checks if the file has 'single_cells_and_organells.csv' in its name

            filename = os.path.basename(file) #gives you the file name
            filename_string = str(filename) #turns the filename into a string

            filedir_string = str(filedir)[:-len(filename) - 5]  #gives you the file directory as a string

            filedir_images_string = filedir_string + r'\neigboorhood'  #creates a new directory

            if os.path.isdir(filedir_images_string) == True: #checks if the new directory allready exist
                pass
            else: #creates a new directory if necessary
                os.makedirs(filedir_images_string)

            neigboorhood_all = neighborhood(pd.read_csv(filedir))  #calculates the spatial proximity of all cells within the tissue
            neigboorhood_all.to_csv(filedir_images_string + '/' + filename_string[:-18] + '_neighborhood.csv')  #exports the file
