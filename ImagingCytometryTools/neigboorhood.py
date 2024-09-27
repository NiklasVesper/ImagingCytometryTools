import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import re

warnings.simplefilter(action='ignore', category=FutureWarning)

def neigboorhood(Cells):
    print('I am looking for my friends :)')

    # create unique list of names of images
    UniqueNames_Full_cell = Cells.ImageNumber.unique()  # creates a list of all unique images

    # create a data frame dictionary to store your data frames
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  # creates a dictionary of all unique images

    # Resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    Neigboorhood = []
    Cell_number = []

    # Takes connectet images and they identefied objects form from the Dictionary
    for key, value in DataFrameDict_Full_cell.items():

        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        count_cell = -1
        for cell in position_cells:
            cell_position = shapely.geometry.Point(cell)
            neighborhood_radius = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].mean()
            cell_neighborhood = shapely.geometry.Point(cell).buffer(neighborhood_radius * 1.5)
            count_cell = count_cell + 1
            Cell_number.append(count_cell)

            Whats_in_the_hood = []
            count_other_cell = -1
            for other_cell in position_cells:
                other_cell_position = shapely.geometry.Point(other_cell)
                count_other_cell = count_other_cell + 1

                if other_cell_position == cell_position:
                    continue
                elif other_cell_position.within(cell_neighborhood) == True:
                    Whats_in_the_hood.append(count_other_cell)

            Neigboorhood.append(Whats_in_the_hood)

    Cells['Neigboorhood'] = Neigboorhood
    Cells['Cell_number'] = Cell_number

    return(Cells)
