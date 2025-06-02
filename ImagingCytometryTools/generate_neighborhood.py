import pandas as pd
import numpy as np
import shapely.geometry
import warnings

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
            neighborhood_radius = Cells['AreaShape_MaxFeretDiameter'].mean() #gives you the mean MaxFeretDiameter of all cells
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
