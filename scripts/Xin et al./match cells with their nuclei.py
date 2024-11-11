import pandas as pd
import os
import scandir as sd

import numpy as np
import shapely.geometry
import warnings
import re

'''
The "match cells with their nuclei" script links objects that were independently identified in a CellProfiler pipeline.

The CellProfiler pipeline can be found under:
https://github.com/MrTheLukasBoom/ImagingCytometryTools/blob/main/Cellprofiler%20pipelines/Xin%20et%20al.%20Segmentation%20Deeplearning%20Imaging%20Mass%20Cytometry.cpproj

For this CellProfiler pipeline Cellpose 2 is required as a plugin. 
Instructions for the plugin installation can be found under:
https://github.com/CellProfiler/CellProfiler-plugins/tree/master

CellProfiler: https://doi.org/10.1186/s12859-021-04344-9
Cellpose 2.0: https://doi.org/10.1038/s41592-022-01663-4
'''

warnings.simplefilter(action='ignore', category=FutureWarning)

#gets protein markers from a CellProfiler file
def get_markers_from_segmentation(Cells):
    first_column = list(Cells.columns.values.tolist()) #gets you the keys from a pandas data frame as list
    clean_markers = []  #empty list for markers

    for element in first_column:  #gets you each element in the key list
        marker_cond = re.findall('MeanIntensity_', element)  #condition to find the markers

        if bool(marker_cond) == True:  #if there is an element that fits the description do the following
            marker = element[element.find('MeanIntensity_'):]  #find all the symbols after mean intensity
            marker_str = str(marker)  #turns marker into strings for the list
            clean_markers.append(marker_str)  #append to marker list

    return (clean_markers)  #returns marker list


#match cells with their nuclei for mostly round cells
def cell_to_organell_basic(Cells, Cytoplasm, Nucleus):
    print('I am working :)')

    #creates lists of all unique images for the different localizations
    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique()
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique()

    #creates a pandas data frame dictionary of all unique images for the different localizations
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm}
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus}

    #creates column name list for the pandas data frame that stores the relevant information of matched subcellular locations
    columns = ['ImageNumber', 'Location_Center_X', 'Location_Center_Y', 'AreaShape_MinFeretDiameter', 'AreaShape_MaxFeretDiameter', ]
    for markers in get_markers_from_segmentation(Cells):
        columns.append(markers + '_Cell')
        columns.append(markers + '_Cytoplasm')
        columns.append(markers + '_Nucleus')

    single_nuclear_cells = pd.DataFrame(columns=columns)  #empty pandas data frame that stores the relevant information of matched subcellular locations

    #resets indexes for identified objects in the image blocks
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()

    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    #takes connected images and they identified objects form from the dictionary
    for key, value in DataFrameDict_Full_cell.items():

        #gets all x and y positions of cells on an image and creates an array of them
        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        #gets all x and y positions of cytoplasm on an image and creates an array of them
        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list()
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list()
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t)))

        #gets all x and y positions of nuclei on an image and creates an array of them
        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list()
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list()
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t)))

        count_cell = -1  #counter
        for cell in position_cells:  #goes throw all cell positions
            cell_diameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]  #creates the smallest possible cell diameter
            cell_shape = shapely.geometry.Point(cell).buffer(cell_diameter / 2)  #creates a circle representing the cell area
            count_cell = count_cell + 1  #ticks up the counter

            count_cytoplasm = -1  #counter
            for cytoplasm in position_cytoplasm:  #goes throw all cytoplasm positions
                cytoplasm_position = shapely.geometry.Point(cytoplasm)  #creates a point at the center of the cytoplasm object
                count_cytoplasm = count_cytoplasm + 1  #ticks up the counter

                if cytoplasm_position.within(cell_shape) == True:  #checks if the center of the cytoplasm is within the cell

                    count_nucleus = -1  #counter
                    nucleus_list = []  #list for all found nuclei
                    for nucleus in position_nucleus:  #goes throw all nuclei positions
                        nucleus_position = shapely.geometry.Point(nucleus) #creates a point at the center of the nucleus object
                        count_nucleus = count_nucleus + 1  #ticks up the counter

                        if nucleus_position.within(cell_shape) == True:  #checks if the center of the nucleus is within the cell
                            nucleus_list.append(count_nucleus) #appends the nucleus count to the list for all found nuclei

                    if len(nucleus_list) == 1:  #checks the length of the nucleus list

                        new_row = {}  #a new row for the pandas data frame

                        ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                        Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                        Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                        AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                        AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]

                        #adds relevant information to the new row
                        new_row.update({'ImageNumber': ImageNumber,
                                        'Location_Center_X': Location_Center_X,
                                        'Location_Center_Y': Location_Center_Y,
                                        'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                        'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter})

                        #adds the mean intensities of the full cell and subcellular localizations to the new row
                        for marker in get_markers_from_segmentation(Cells):
                            MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                            MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                            MeanIntensity_Nucleus = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]

                            new_row.update({marker + '_Cell': MeanIntensity_Cell,
                                            marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                            marker + '_Nucleus': MeanIntensity_Nucleus})

                        single_nuclear_cells.loc[len(single_nuclear_cells)] = new_row  #adds the new row to the pandas data frame
                    break

    return (single_nuclear_cells)


folder_dir = r''  #folder directory

for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file)  #returns full file directory

        if filedir.endswith(".csv"):  #returns all files that end with ".csv"

            filename = os.path.basename(file) #gives you the file name
            filename_string = str(filename) #turns the filename into a string

            if 'RunCellpose_C' in filename_string: #checks if 'RunCellpose_C' is in the string
                dir_list = []  #a list to hold the other directories

                dir_list.append(filedir)
                filedir_string = str(filedir)[:-len(filename) - 5]  #gives you the file directory as a string

                filedir_images_string = filedir_string + r'\subcellular localisation'  #creates a new directory

                if os.path.isdir(filedir_images_string) == True:  #checks if the new directory allready exists
                    pass
                else: #creates a new directory if necessary
                    os.makedirs(filedir_images_string)

                for paths, dirs, files in sd.walk(filedir_string + '/csv'):  #goes throw all files and folders in given directory

                    for file in os.listdir(paths):  #goes throw all files in a folder
                        filedir = os.path.join(paths, file) #returns full file directory
                        filedir_string = str(filedir)  # creates a file directory string

                        if 'RunCellpose_N' in filedir_string:  #checks if 'RunCellpose_N' is in the string
                            dir_list.append(filedir)
                        if 'Cytoplasm' in filedir_string:  #checks if 'Cytoplasm' is in the string
                            dir_list.append(filedir)

                        if len(dir_list) == 3:  #if three files are found the matching is performed

                            print(dir_list)

                            #reads in the found files
                            Full_cell = pd.read_csv(dir_list[0])  
                            Cytoplasm = pd.read_csv(dir_list[1])
                            Nucleus = pd.read_csv(dir_list[2])

                            single_cells_and_organells = cell_to_organell_basic(Full_cell, Cytoplasm, Nucleus,)  #links objects
                            single_cells_and_organells.to_csv(filedir_images_string + '/' + filename_string[:-18] + '_single_cells_and_organells.csv') #exports the file
