import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import os
import scandir as sd
import re

warnings.simplefilter(action='ignore', category=FutureWarning)

from ImagingCytometryToolsIMC.get_markers import get_markers_from_segmentation

# for mostly round cells
def cell_to_organell_basic(Cells, Cytoplasm, Nucleus, Nucleus_count):

    print('I am working :)')

    #create unique list of names of images
    UniqueNames_Full_cell = Cells.ImageNumber.unique() #creates a list of all unique images
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique() #creates a list of all unique images
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique() #creates a list of all unique images

    #create a data frame dictionary to store your data frames
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell} #creates a dictionary of all unique images
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm} #creates a dictionary of all unique images
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus} #creates a dictionary of all unique images

    #creates a data frame to store matched subcellular locations
    columns = ['ImageNumber','Location_Center_X','Location_Center_Y','AreaShape_MinFeretDiameter','AreaShape_MaxFeretDiameter',] #list for essential values and markers
    for markers in get_markers_from_segmentation(Cells): #get all markers
        columns.append(markers + '_Cell') #creates a column for the full cell
        columns.append(markers + '_Cytoplasm') #creates a column for the cytoplasm
        columns.append(markers + '_Nucleus') #creates a column for the nucleus

    single_nuclear_cells = pd.DataFrame(columns=columns) #empty dataframe for all markers and subcellular locations

    #Resets indexes for image blocks and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()

    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    #Takes connectet images and they identefied objects form from the Dictionary
    for key, value in DataFrameDict_Full_cell.items():

        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list() #gets all x positions of cells on a image
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list() #gets all y positions of cells on a image
        position_cells = np.array(list(zip(x_cell_t, y_cell_t))) #creates an array of all the x and y positions of the cells

        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list() #gets all x positions of cytoplasms on a image
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list() #gets all y positions of cytoplasms on a image
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t))) #creates an array of all the x and y positions of the cytoplasms

        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list() #gets all x positions of nuclei on a image
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list() #gets all y positions of nuclei on a image
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t))) #creates an array of all the x and y positions of the nuclei
        
        count_cell = -1 #counter
        for cell in position_cells: #goes throw all cell positions
            cell_diameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell] #creates the smallest possible cell diameter
            cell_shape = shapely.geometry.Point(cell).buffer(cell_diameter/2) #creates a circle representing the cell area
            count_cell = count_cell + 1 #adds a counter to the counter

            count_cytoplasm = -1 #counter
            for cytoplasm in position_cytoplasm: #goes throw all cytoplasm positions
                cytoplasm_position = shapely.geometry.Point(cytoplasm) #creates a point at the center of the cytoplasm object
                count_cytoplasm = count_cytoplasm + 1 #adds a counter to the counter

                if cytoplasm_position.within(cell_shape) == True:

                    count_nucleus = -1 #counter
                    nucleus_list = [] #list of all found nuclei
                    for nucleus in position_nucleus: #goes throw all nuclei positions
                        nucleus_position = shapely.geometry.Point(nucleus)
                        count_nucleus = count_nucleus + 1 #adds a counter to the counter

                        if nucleus_position.within(cell_shape) == True:
                            nucleus_list.append(count_nucleus)

                    if len(nucleus_list) == Nucleus_count:

                        new_row = {}  # Creates a new Row for all the localistions

                        ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                        Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                        Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                        AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                        AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]

                        new_row.update({'ImageNumber': ImageNumber,
                                        'Location_Center_X': Location_Center_X,
                                        'Location_Center_Y': Location_Center_Y,
                                        'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                        'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter})

                        for marker in get_markers_from_segmentation(Cells):
                            MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                            MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                            MeanIntensity_Nucleus = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]

                            new_row.update({marker + '_Cell': MeanIntensity_Cell,
                                       marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                       marker + '_Nucleus': MeanIntensity_Nucleus})

                        single_nuclear_cells.loc[len(single_nuclear_cells)] = new_row
                    break
                    
    return(single_nuclear_cells)
