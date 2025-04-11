import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import os
import scandir as sd
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=FutureWarning)

from ImagingCytometryTools.get_stuff import get_markers_from_segmentation
from ImagingCytometryTools.get_stuff import get_pixel_values

from ImagingCytometryTools.object_pixel_extraction import get_connected_black_pixels
from ImagingCytometryTools.object_pixel_extraction import find_edge_pixels
from ImagingCytometryTools.object_pixel_extraction import pixel_sort

'''
Functions that match different subcellular localizations together.

Here two options are available.
The first function extracts the full cellular area and check whether the center of the nucleus is within that area.
The second function checks whether the full nuclear area is within the cell.
'''

#extracts the full cellular area and check whether the center of the nucleus is within that area
def cell_to_organell_basic(Cells, Cytoplasm, Nucleus, FileDirectory, ImageName, outline_cell_image_dir, outline_nucleus_image_dir,Nucleus_count):

    #create unique list of names of images
    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique()
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique()

    #create a data frame dictionary to store your data frames
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm}
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus}

    #creates a list to store matched subcellular locations and all the related information
    if Nucleus_count == 1:
        #list for essential values and markers
        columns = ['FileDirectory',
                   'ImageName',
                   'ImageNumber',
                   'Location_Center_X',
                   'Location_Center_Y',
                   'AreaShape_MinFeretDiameter',
                   'AreaShape_MaxFeretDiameter',
                   'Cell_outline',
                   'Cell_pixel_area',
                   'Nuclear_outline',
                   'Nuclear_pixel_area']
        for markers in get_markers_from_segmentation(Cells):  #gets protein markers from the file
            columns.append('PathName_' + markers.replace('MeanIntensity_', '')) #path name of the image
            columns.append('FileName_' + markers.replace('MeanIntensity_', '')) #file name of the image
            columns.append('Intensity_' + markers + '_Cell')  #creates a column for the full cell
            columns.append('Intensity_' + markers + '_Cytoplasm')  #creates a column for the cytoplasm
            columns.append('Intensity_' + markers + '_Nucleus')  #creates a column for the nucleus

    #creates a list to store matched subcellular locations and all the related information of dividing cells
    elif Nucleus_count == 2:
        #list for essential values and markers
        columns = ['FileDirectory',
                   'ImageName',
                   'ImageNumber',
                   'Location_Center_X',
                   'Location_Center_Y',
                   'AreaShape_MinFeretDiameter',
                   'AreaShape_MaxFeretDiameter',
                   'Cell_outline',
                   'Cell_pixel_area',
                   'Nuclear_outline_1',
                   'Nuclear_pixel_area_1',
                   'Nuclear_outline_2',
                   'Nuclear_pixel_area_2']
        for markers in get_markers_from_segmentation(Cells):  #gets protein markers from the file
            columns.append('PathName_' + markers.replace('MeanIntensity_', '')) #path name of the image
            columns.append('FileName_' + markers.replace('MeanIntensity_', '')) #file name of the image
            columns.append('Intensity_' + markers + '_Cell')  #creates a column for the full cell
            columns.append('Intensity_' + markers + '_Cytoplasm')  #creates a column for the cytoplasm
            columns.append('Intensity_' + markers + '_Nucleus_1')  #creates a column for the first nucleus
            columns.append('Intensity_' + markers + '_Nucleus_2') #creates a column for the second nucleus

    cells_with_sub_cell_loc = pd.DataFrame(columns=columns)  #empty dataframe for subcellular locations and all the related information

    #resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()
    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()
    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    file_dir_list_cell = [] #list for the cell image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  #goes throw all files and folders in given directory

        for file in os.listdir(paths):  #goes throw all files in a folder
            filedir = os.path.join(paths, file)  #returns full file directory

            if filedir.endswith(".tiff"):  #returns all files that end with .tiff
                file_dir_list_cell.append(filedir) #appends the found .tiff files into the list

    file_dir_list_nucleus = [] #list for the nucleus image outline directories
    for paths, dirs, files in sd.walk(outline_nucleus_image_dir):  #goes throw all files and folders in given directory

        for file in os.listdir(paths):  #goes throw all files in a folder
            filedir = os.path.join(paths, file)  #returns full file directory

            if filedir.endswith(".tiff"):  #returns all files that end with .tiff
                file_dir_list_nucleus.append(filedir) #appends the found .tiff files into the list

    for key, value in tqdm(DataFrameDict_Full_cell.items()): #takes connected images and they identified objects form from the dictionary

        #extracts x and y position for each individual cell
        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        #extracts x and y position for each individual cytoplasm
        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list()
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list()
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t)))

        #extracts x and y position for each individual nucleus
        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list()
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list()
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1]) #gets the cell outline image associated with the image key
        outline_image_nucleus = get_pixel_values(file_dir_list_nucleus[key - 1]) #gets the nucleus outline image associated with the image key

        count_cell = -1  #cell counter
        for cell in position_cells:
            count_cell = count_cell + 1  # add 1 to the cell counter
            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0]))

            #if the start pixel is white it could happen that we jump outside the cell boundaries so we need to exclude these cells
            if connected_black_pixels_cell == None:
                continue

            edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method = 'edge')

            # excludes the cell if the area is below 4 because you can not construct a polygon with it
            if len(edge_pixels_cell) < 4:
                continue

            elif len(edge_pixels_cell) >=4:

                pixel_sorted_cell = pixel_sort(edge_pixels_cell, method = 'frog_jump')
                cell_shape = shapely.geometry.Polygon(pixel_sorted_cell)

                count_cytoplasm = -1  #cytoplasm counter
                for cytoplasm in position_cytoplasm:
                    count_cytoplasm = count_cytoplasm + 1  # add 1 to the cytoplasm counter
                    cytoplasm_position = shapely.geometry.Point(cytoplasm)  #creates a point at the center of the cytoplasm object

                    if cytoplasm_position.within(cell_shape) == True:  #checks if the center of the cytoplasm is within the cell area

                        count_nucleus = -1  #nucleus counter
                        nucleus_list = []  #list for all found nuclei
                        for nucleus in position_nucleus:
                            count_nucleus = count_nucleus + 1  # add 1 to the nucleus counter
                            nucleus_position = shapely.geometry.Point(nucleus)

                            if nucleus_position.within(cell_shape) == True:  #checks if the center of the nucleus is within the cell area
                                nucleus_list.append(count_nucleus)

                        if len(nucleus_list) == 1 and Nucleus_count == 1:

                            connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))

                            if connected_black_pixels_nucleus == None:
                                continue

                            else:
                                edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method = 'edge')
                                pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method= 'frog_jump')

                                new_row = {}  #creates a new row for all the localizations

                                #adds all the relevant information
                                ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                                Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                                Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                                AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                                AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]
                                new_row.update({'FileDirectory': FileDirectory,
                                                'ImageName': ImageName,
                                                'ImageNumber': ImageNumber,
                                                'Location_Center_X': Location_Center_X,
                                                'Location_Center_Y': Location_Center_Y,
                                                'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                                'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter,
                                                'Cell_outline': pixel_sorted_cell,
                                                'Cell_pixel_area': connected_black_pixels_cell,
                                                'Nuclear_outline': pixel_sorted_nucleus,
                                                'Nuclear_pixel_area': connected_black_pixels_nucleus},
                                               )

                                #adds the mean intensities to the new row
                                for marker in get_markers_from_segmentation(Cells):
                                    PathName = DataFrameDict_Full_cell[key]['PathName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                    FileName = DataFrameDict_Full_cell[key]['FileName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                    MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                                    MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                                    MeanIntensity_Nucleus = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]
                                    new_row.update({
                                        'PathName_' + marker.replace('MeanIntensity_', ''): PathName,
                                        'FileName_' + marker.replace('MeanIntensity_', ''): FileName,
                                        'Intensity_' + marker + '_Cell': MeanIntensity_Cell,
                                        'Intensity_' + marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                        'Intensity_' + marker + '_Nucleus': MeanIntensity_Nucleus})

                                cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row  #adds the new row to the data frame

                        break

                        if len(nucleus_list) == 2 and Nucleus_count == 2:  # checks the lenth of the nucleus list

                            connected_black_pixels_nucleus_1 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                            edge_pixels_nucleus_1 = find_edge_pixels(connected_black_pixels_nucleus_1, method = 'edge')
                            pixel_sorted_nucleus_1 = pixel_sort(edge_pixels_nucleus_1, method= 'frog_jump')

                            connected_black_pixels_nucleus_2 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[1]][1]),round(position_nucleus[nucleus_list[1]][0]))
                            edge_pixels_nucleus_2 = find_edge_pixels(connected_black_pixels_nucleus_2, method = 'edge')
                            pixel_sorted_nucleus_2 = pixel_sort(edge_pixels_nucleus_2, method= 'frog_jump')

                            new_row = {}  #creates a new row for all the localizations

                            #adds all the relevant information
                            ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                            Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                            Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                            AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                            AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]
                            new_row.update({'FileDirectory': FileDirectory,
                                            'ImageName': ImageName,
                                            'ImageNumber': ImageNumber,
                                            'Location_Center_X': Location_Center_X,
                                            'Location_Center_Y': Location_Center_Y,
                                            'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                            'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter,
                                            'Cell_outline': pixel_sorted_cell,
                                            'Cell_pixel_area': connected_black_pixels_cell,
                                            'Nuclear_outline_1': pixel_sorted_nucleus_1,
                                            'Nuclear_pixel_area_1': connected_black_pixels_nucleus_1,
                                            'Nuclear_outline_2': pixel_sorted_nucleus_2,
                                            'Nuclear_pixel_area_2': connected_black_pixels_nucleus_2},
                                           )

                            #adds the mean intensities to the new row
                            for marker in get_markers_from_segmentation(Cells):
                                PathName = DataFrameDict_Full_cell[key]['PathName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                FileName = DataFrameDict_Full_cell[key]['FileName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                                MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                                MeanIntensity_Nucleus_1 = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]
                                MeanIntensity_Nucleus_2 = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[1]]
                                new_row.update({
                                    marker.replace('MeanIntensity_', '') + '_PathName': PathName,
                                    marker.replace('MeanIntensity_', '') + '_FileName': FileName,
                                    marker + '_Cell': MeanIntensity_Cell,
                                    marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                    marker + '_Nucleus_1': MeanIntensity_Nucleus_1,
                                    marker + '_Nucleus_2': MeanIntensity_Nucleus_2})

                            cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row  #adds the new row to the data frame

                        break

    return(cells_with_sub_cell_loc)

#extracts the full cellular area and check whether the area of the nucleus is within the cellular area
def cell_to_organell_advanced(Cells, Cytoplasm, Nucleus, FileDirectory, ImageName, outline_cell_image_dir, outline_nucleus_image_dir,Nucleus_count):

    #create unique list of names of images
    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique()
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique()

    #create a data frame dictionary to store your data frames
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm}
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus}

    #creates a list to store matched subcellular locations and all the related information
    # creates a list to store matched subcellular locations and all the related information
    if Nucleus_count == 1:
        columns = ['FileDirectory', 'ImageName', 'ImageNumber', 'Location_Center_X', 'Location_Center_Y', 'AreaShape_MinFeretDiameter','AreaShape_MaxFeretDiameter', 'Cell_outline', 'Cell_pixel_area', 'Nuclear_outline','Nuclear_pixel_area']  # list for essential values and markers
        for markers in get_markers_from_segmentation(Cells):  # gets protein markers from the file
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))  # path name of the image
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))  # file name of the image
            columns.append('Intensity_' + markers + '_Cell')  # creates a column for the full cell
            columns.append('Intensity_' + markers + '_Cytoplasm')  # creates a column for the cytoplasm
            columns.append('Intensity_' + markers + '_Nucleus')  # creates a column for the nucleus

    # creates a list to store matched subcellular locations and all the related information of dividing cells
    elif Nucleus_count == 2:
        columns = ['FileDirectory','ImageName', 'ImageNumber', 'Location_Center_X', 'Location_Center_Y', 'AreaShape_MinFeretDiameter','AreaShape_MaxFeretDiameter', 'Cell_outline', 'Cell_pixel_area', 'Nuclear_outline_1','Nuclear_pixel_area_1', 'Nuclear_outline_2','Nuclear_pixel_area_2']  # list for essential values and markers
        for markers in get_markers_from_segmentation(Cells):  # gets protein markers from the file
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))  # path name of the image
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))  # file name of the image
            columns.append('Intensity_' + markers + '_Cell')  # creates a column for the full cell
            columns.append('Intensity_' + markers + '_Cytoplasm')  # creates a column for the cytoplasm
            columns.append('Intensity_' + markers + '_Nucleus_1')  # creates a column for the first nucleus
            columns.append('Intensity_' + markers + '_Nucleus_2')  # creates a column for the second nucleus

    cells_with_sub_cell_loc = pd.DataFrame(columns=columns)  # empty dataframe for subcellular locations and all the related information

    #resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()
    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()
    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    file_dir_list_cell = [] #list for the cell image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  #goes throw all files and folders in given directory

        for file in os.listdir(paths):  #goes throw all files in a folder
            filedir = os.path.join(paths, file)  #returns full file directory

            if filedir.endswith(".tiff"):  #returns all files that end with .tiff
                file_dir_list_cell.append(filedir) #appends the found .tiff files into the list

    file_dir_list_nucleus = [] #list for the nucleus image outline directories
    for paths, dirs, files in sd.walk(outline_nucleus_image_dir):  #goes throw all files and folders in given directory

        for file in os.listdir(paths):  #goes throw all files in a folder
            filedir = os.path.join(paths, file)  #returns full file directory

            if filedir.endswith(".tiff"):  #returns all files that end with .tiff
                file_dir_list_nucleus.append(filedir) #appends the found .tiff files into the list

    for key, value in tqdm(DataFrameDict_Full_cell.items(), colour='white'): #takes connected images and they identified objects form from the dictionary

        #extracts x and y position for each individual cell
        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        #extracts x and y position for each individual cytoplasm
        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list()
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list()
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t)))

        #extracts x and y position for each individual nucleus
        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list()
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list()
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1]) #gets the cell outline image associated with the image key
        outline_image_nucleus = get_pixel_values(file_dir_list_nucleus[key - 1]) #gets the nucleus outline image associated with the image key

        count_cell = -1  #cell counter
        for cell in position_cells:
            count_cell = count_cell + 1  # add 1 to the cell counter

            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0]))
            edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method = 'edge')

            # excludes the cell if the area is below 4 because you can not construct a polygon with it
            if len(edge_pixels_cell) < 4:
                continue

            elif len(edge_pixels_cell) >=4:

                pixel_sorted_cell = pixel_sort(edge_pixels_cell, method= 'frog_jump')
                cell_shape = shapely.geometry.Polygon(pixel_sorted_cell)

                count_cytoplasm = -1  #cytoplasm counter
                for cytoplasm in position_cytoplasm:
                    count_cytoplasm = count_cytoplasm + 1  #add 1 to the cytoplasm counter
                    cytoplasm_position = shapely.geometry.Point(cytoplasm)  #creates a point at the center of the cytoplasm object

                    if cytoplasm_position.within(cell_shape) == True:  #checks if the center of the cytoplasm is within the cell area

                        count_nucleus = -1  #nucleus counter
                        nucleus_list = []  #list for all found nuclei
                        for nucleus in position_nucleus:
                            count_nucleus = count_nucleus + 1  #add 1 to the nucleus counter

                            connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus, round(nucleus[1]),round(nucleus[0]))
                            edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method='edge')
                            pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method= 'frog_jump')

                            nucleus_shape = shapely.geometry.Polygon(pixel_sorted_nucleus)

                            if nucleus_shape.within(cell_shape) == True:  #checks if the area of the nucleus is within the cell area
                                nucleus_list.append(count_nucleus)

                        if len(nucleus_list) == 1 and Nucleus_count == 1:

                            connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                            edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method = 'edge')
                            pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method= 'frog_jump')

                            new_row = {}  #creates a new row for all the localizations

                            #adds all the relevant information
                            ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                            Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                            Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                            AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                            AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]
                            new_row.update({'FileDirectory': FileDirectory,
                                            'ImageName': ImageName,
                                            'ImageNumber': ImageNumber,
                                            'Location_Center_X': Location_Center_X,
                                            'Location_Center_Y': Location_Center_Y,
                                            'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                            'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter,
                                            'Cell_outline': pixel_sorted_cell,
                                            'Cell_pixel_area': connected_black_pixels_cell,
                                            'Nuclear_outline': pixel_sorted_nucleus,
                                            'Nuclear_pixel_area': connected_black_pixels_nucleus},
                                           )

                            #adds the mean intensities to the new row
                            for marker in get_markers_from_segmentation(Cells):
                                PathName = DataFrameDict_Full_cell[key]['PathName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                FileName = DataFrameDict_Full_cell[key]['FileName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                                MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                                MeanIntensity_Nucleus = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]
                                new_row.update({
                                    'PathName_' + marker.replace('MeanIntensity_', ''): PathName,
                                    'FileName_' + marker.replace('MeanIntensity_', ''): FileName,
                                    'Intensity_' + marker + '_Cell': MeanIntensity_Cell,
                                    'Intensity_' + marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                    'Intensity_' + marker + '_Nucleus': MeanIntensity_Nucleus})

                            cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row  #adds the new row to the data frame

                        break

                        if len(nucleus_list) == 2 and Nucleus_count == 2:  # checks the lenth of the nucleus list

                            connected_black_pixels_nucleus_1 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                            edge_pixels_nucleus_1 = find_edge_pixels(connected_black_pixels_nucleus_1, method = 'edge')
                            pixel_sorted_nucleus_1 = pixel_sort(edge_pixels_nucleus_1, method= 'frog_jump')

                            connected_black_pixels_nucleus_2 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[1]][1]),round(position_nucleus[nucleus_list[1]][0]))
                            edge_pixels_nucleus_2 = find_edge_pixels(connected_black_pixels_nucleus_2, method = 'edge')
                            pixel_sorted_nucleus_2 = pixel_sort(edge_pixels_nucleus_2, method= 'frog_jump')

                            new_row = {}  #creates a new row for all the localizations

                            #adds all the relevant information
                            ImageNumber = DataFrameDict_Full_cell[key]['ImageNumber'].to_list()[count_cell]
                            Location_Center_X = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()[count_cell]
                            Location_Center_Y = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()[count_cell]
                            AreaShape_MinFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MinFeretDiameter'].to_list()[count_cell]
                            AreaShape_MaxFeretDiameter = DataFrameDict_Full_cell[key]['AreaShape_MaxFeretDiameter'].to_list()[count_cell]
                            new_row.update({'FileDirectory': FileDirectory,
                                            'ImageName': ImageName,
                                            'ImageNumber': ImageNumber,
                                            'Location_Center_X': Location_Center_X,
                                            'Location_Center_Y': Location_Center_Y,
                                            'AreaShape_MinFeretDiameter': AreaShape_MinFeretDiameter,
                                            'AreaShape_MaxFeretDiameter': AreaShape_MaxFeretDiameter,
                                            'Cell_outline': pixel_sorted_cell,
                                            'Cell_pixel_area': connected_black_pixels_cell,
                                            'Nuclear_outline_1': pixel_sorted_nucleus_1,
                                            'Nuclear_pixel_area_1': connected_black_pixels_nucleus_1,
                                            'Nuclear_outline_2': pixel_sorted_nucleus_2,
                                            'Nuclear_pixel_area_2': connected_black_pixels_nucleus_2},
                                           )

                            #adds the mean intensities to the new row
                            for marker in get_markers_from_segmentation(Cells):
                                PathName = DataFrameDict_Full_cell[key]['PathName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                FileName = DataFrameDict_Full_cell[key]['FileName_' + marker.replace('MeanIntensity_', '')].to_list()[count_cell]
                                MeanIntensity_Cell = DataFrameDict_Full_cell[key]['Intensity_' + marker].to_list()[count_cell]
                                MeanIntensity_Cytoplasm = DataFrameDict_Cytoplasm[key]['Intensity_' + marker].to_list()[count_cytoplasm]
                                MeanIntensity_Nucleus_1 = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[0]]
                                MeanIntensity_Nucleus_2 = DataFrameDict_Nucleus[key]['Intensity_' + marker].to_list()[nucleus_list[1]]
                                new_row.update({
                                    marker.replace('MeanIntensity_', '') + '_PathName': PathName,
                                    marker.replace('MeanIntensity_', '') + '_FileName': FileName,
                                    marker + '_Cell': MeanIntensity_Cell,
                                    marker + '_Cytoplasm': MeanIntensity_Cytoplasm,
                                    marker + '_Nucleus_1': MeanIntensity_Nucleus_1,
                                    marker + '_Nucleus_2': MeanIntensity_Nucleus_2})

                            cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row  #adds the new row to the data frame
                        break

    return (cells_with_sub_cell_loc)