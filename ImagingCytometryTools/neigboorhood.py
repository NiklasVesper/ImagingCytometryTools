import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import scandir as sd
import os
import ast
from tqdm import tqdm

from ImagingCytometryTools.get_stuff import get_markers_from_segmentation
from ImagingCytometryTools.get_stuff import get_pixel_values

from ImagingCytometryTools.object_pixel_extraction import get_connected_black_pixels
from ImagingCytometryTools.object_pixel_extraction import find_edge_pixels
from ImagingCytometryTools.object_pixel_extraction import pixel_sort

warnings.simplefilter(action='ignore', category=FutureWarning)

'''
Functions that analyze the surrounding cellular neighborhood of cells identified with the segmentation.
'''

#calculates the neighboring cells based on distance
def neigboorhood_add_outline(Cells, Radius, FileDirectory, ImageName, outline_cell_image_dir):

    #creates a list for all the information that is kept
    columns = ['ImageNumber',
               'Location_Center_X',
               'Location_Center_Y',
               'AreaShape_MinFeretDiameter',
               'AreaShape_MaxFeretDiameter']
    for markers in get_markers_from_segmentation(Cells):  #gets protein markers from the file
        columns.append('PathName_' + markers.replace('MeanIntensity_', '')) #path name of the image
        columns.append('FileName_' + markers.replace('MeanIntensity_', '')) #file name of the image
        columns.append('Intensity_' + markers)

    Cells_clean = Cells.loc[:, columns] #removes unwanted columbs
    Cells_clean.insert(0,'ImageName',ImageName) #adds the image name
    Cells_clean.insert(0,'FileDirectory', FileDirectory)  # adds the image name

    UniqueNames_Full_cell = Cells_clean.ImageNumber.unique()  #creates a list of all unique images
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  #creates a dictionary of all unique images

    #resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells_clean[:][Cells_clean.ImageNumber == key].reset_index()

    file_dir_list_cell = [] #list for the image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  #goes throw all files and folders in given directory

        for file in os.listdir(paths):  #goes throw all files in a folder
            filedir = os.path.join(paths, file)  #returns full file directory

            if filedir.endswith(".tiff"):  #returns all files that end with .tiff
                file_dir_list_cell.append(filedir) #appends the found .tiff files into the list

    Cell_pixel_area = [] #empty list for cell pixel area
    Cell_outline = [] #empty list for cell outlines

    for key, value in DataFrameDict_Full_cell.items(): #takes connected images and they identified objects form from the dictionary

        #extracts x and y position for each individual cell
        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1]) #gets the outline image associated with the image key

        #goes throw each cell position
        for cell in position_cells:

            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0])) #gets the connected black pixels within the cell outline

            #if the start pixel is white it could happen that we jump outside the cell boundaries so we need to exclude these cells
            if connected_black_pixels_cell == None:
                Cell_pixel_area.append('to_remove')
                Cell_outline.append('to_remove')

            else:
                Cell_pixel_area.append(connected_black_pixels_cell) #appends all found black pixels as cell area
                edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method = 'edge') #finds all pixels that surround the cell area that mark the outline
                pixel_sorted_cell = pixel_sort(edge_pixels_cell, method = 'frog_jump') #sorts the pixels clockwise so they can be uses as a polygon in shapely
                Cell_outline.append(pixel_sorted_cell) #appends outline pixels

    Cells_clean['Cell_outline'] = Cell_outline #adds outline pixels to the data frame
    Cells_clean['Cell_pixel_area'] = Cell_pixel_area #adds area pixels to the data frame
    Cells_clean = Cells_clean.drop(Cells_clean[Cells_clean['Cell_pixel_area'] == 'to_remove'].index) #adds outline pixels to the data frame

    #reinitializes the data frame dictionary
    UniqueNames_Full_cell = Cells_clean.ImageNumber.unique()  #creates a list of all unique images
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  #creates a dictionary of all unique images

    #resets indexes for conected images and they identefied object
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells_clean[:][Cells_clean.ImageNumber == key].reset_index()

    Cell_number = [] #list for individual cell number
    Neigboorhood = [] #list for all the neighboring cells
    Neigboorhood_coordinates = [] #list for the neigboorhood coordinates

    for key, value in tqdm(DataFrameDict_Full_cell.items()): #takes connected images and they identified objects form from the dictionary

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list() #makes a list out of the cell outlines

        count_cell = -1 #cell counter
        #goes throw the cellular outlines
        for cell in outline_cell:
            count_cell = count_cell + 1 #add 1 to the cell counter
            Cell_number.append(count_cell) #assigns the cell a number

            #excludes the cell if the area is below 4 because you can not construct a polygon with it
            if len(cell) < 4:
                Neigboorhood.append('Cell size below 4 pixels')
                Neigboorhood_coordinates.append('Cell size below 4 pixels')

            elif len(cell) >=4:
                cell_area = shapely.geometry.Polygon(np.array(cell)) #turns the sorted edge pixels into a polygon

                #adds the neighborhood radius to the cellular area
                neighborhood_radius = Radius
                cell_neighborhood = cell_area.buffer(neighborhood_radius)

                #gets the buffer coordinates and stores them in a separate list
                cell_neighborhood_coordinates = [tuple(round(value) for value in tup) for tup in tuple(cell_neighborhood.exterior.coords)]
                Neigboorhood_coordinates.append(cell_neighborhood_coordinates)

                Whats_in_the_hood = [] #list for the individual neighborhood of each cell
                count_other_cell = -1 #counter for the other cells
                for other_cell in outline_cell:

                    count_other_cell = count_other_cell + 1 #add 1 to the other cell counter

                    #excludes the cell if the area is below 4 because you can not construct a polygon with it
                    if len(other_cell) < 4:
                        continue

                    elif len(other_cell) >=4:

                        #turns the sorted edge pixels into a polygon
                        outline_cell_np = np.array(other_cell)
                        other_cell_outline = shapely.geometry.Polygon(outline_cell_np)

                        #excludes the cell itself within the neighborhood
                        if other_cell_outline == cell_area:
                            continue

                        #adds the cell to the neighborhood if it intersects with the cell neighborhood
                        elif cell_neighborhood.intersects(other_cell_outline) == True:
                            Whats_in_the_hood.append(count_other_cell)

                Neigboorhood.append(Whats_in_the_hood)

    #adds the neighborhood to the data frame
    Cells_clean['Neighborhood'] = Neigboorhood
    Cells_clean['Cell_number'] = Cell_number
    Cells_clean['Neighborhood_coordinates'] = Neigboorhood_coordinates

    return(Cells_clean)

#calculates the neighboring cells based on distance but does not add any other information
def neigboorhood(Cells,Radius):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()  #creates a list of all unique images

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  #creates a dictionary of all unique images

    #resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    Cell_number = []  #list for individual cell number
    Neigboorhood = []  #list for all the neighboring cells
    Neigboorhood_coordinates = []  # list for the neigboorhood coordinates

    for key, value in tqdm(DataFrameDict_Full_cell.items()):

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list()  # makes a list out of the cell outlines

        count_cell = -1 #cell counter
        #goes throw the cellular outlines
        for cell in outline_cell:
            cell = ast.literal_eval(cell)
            count_cell = count_cell + 1 #add 1 to the cell counter
            Cell_number.append(count_cell) #assigns the cell a number

            #excludes the cell if the area is below 4 because you can not construct a polygon with it
            if len(cell) < 4:
                Neigboorhood.append('Cell size below 4 pixels')
                Neigboorhood_coordinates.append('Cell size below 4 pixels')

            elif len(cell) >= 4:

                cell_area = shapely.geometry.Polygon(np.array(cell)) #turns the sorted edge pixels into a polygon

                #adds the neighborhood radius to the cellular area
                neighborhood_radius = Radius
                cell_neighborhood = cell_area.buffer(neighborhood_radius)

                cell_neighborhood_coordinates = [tuple(round(value) for value in tup) for tup in tuple(cell_neighborhood.exterior.coords)]
                Neigboorhood_coordinates.append(cell_neighborhood_coordinates)

                Whats_in_the_hood = []  #list for the individual neighborhood of each cell
                count_other_cell = -1  #counter for the other cells
                for other_cell in outline_cell:
                    other_cell = ast.literal_eval(other_cell)
                    count_other_cell = count_other_cell + 1 #add 1 to the other cell counter

                    #excludes the cell if the area is below 4 because you can not construct a polygon with it
                    if len(other_cell) < 4:
                        continue

                    elif len(other_cell) >= 4:

                        #turns the sorted edge pixels into a polygon
                        outline_cell_np = np.array(other_cell)
                        other_cell_outline = shapely.geometry.Polygon(outline_cell_np)

                        #excludes the cell itself within the neighborhood
                        if other_cell_outline == cell_area:
                            continue

                        #adds the cell to the neighborhood if it intersects with the cell neighborhood
                        elif cell_neighborhood.intersects(other_cell_outline) == True:
                            Whats_in_the_hood.append(count_other_cell)

                Neigboorhood.append(Whats_in_the_hood)

    #adds the neighborhood to the data frame
    Cells['Neighborhood'] = Neigboorhood
    Cells['Cell_number'] = Cell_number
    Cells['Neighborhood_coordinates'] = Neigboorhood_coordinates

    return(Cells)