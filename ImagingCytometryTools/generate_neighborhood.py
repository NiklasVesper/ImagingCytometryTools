import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import scandir as sd
import os
import ast
from tqdm import tqdm

from ImagingCytometryTools.get_data_from_files import get_markers_from_segmentation
from ImagingCytometryTools.get_data_from_files import get_pixel_values

from ImagingCytometryTools.object_pixel_extraction import get_connected_black_pixels
from ImagingCytometryTools.object_pixel_extraction import find_edge_pixels
from ImagingCytometryTools.object_pixel_extraction import pixel_sort

warnings.simplefilter(action='ignore', category=FutureWarning)

'''
Functions that analyzes the surrounding cellular neighborhood of cells identified with the segmentation.
'''

#calculates the neighboring cells based on distance
def generate_neighborhood_add_outline(Cells, Radius, FileDirectory, ImageName, outline_cell_image_dir):

    columns = ['ImageNumber',
               'Location_Center_X',
               'Location_Center_Y',
               'AreaShape_MinFeretDiameter',
               'AreaShape_MaxFeretDiameter']
    for markers in get_markers_from_segmentation(Cells):
        columns.append('PathName_' + markers.replace('MeanIntensity_', ''))
        columns.append('FileName_' + markers.replace('MeanIntensity_', ''))
        columns.append('Intensity_' + markers)

    Cells_clean = Cells.loc[:, columns]
    Cells_clean.insert(0,'ImageName',ImageName)
    Cells_clean.insert(0,'FileDirectory', FileDirectory)

    UniqueNames_Full_cell = Cells_clean.ImageNumber.unique()
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells_clean[:][Cells_clean.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    Cell_pixel_area = []
    Cell_outline = []

    for key, value in DataFrameDict_Full_cell.items():

        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])

        for cell in position_cells:

            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0]))

            if connected_black_pixels_cell == None:
                Cell_pixel_area.append('to_remove')
                Cell_outline.append('to_remove')

            else:
                edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method = 'edge')
                pixel_sorted_cell = pixel_sort(edge_pixels_cell, method = 'frog_jump')

                if pixel_sorted_cell == None:
                    Cell_pixel_area.append('to_remove')
                    Cell_outline.append('to_remove')

                else:
                    Cell_pixel_area.append(connected_black_pixels_cell)
                    Cell_outline.append(pixel_sorted_cell)

    Cells_clean['Cell_outline'] = Cell_outline
    Cells_clean['Cell_pixel_area'] = Cell_pixel_area
    Cells_clean = Cells_clean.drop(Cells_clean[Cells_clean['Cell_pixel_area'] == 'to_remove'].index)

    UniqueNames_Full_cell = Cells_clean.ImageNumber.unique()
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells_clean[:][Cells_clean.ImageNumber == key].reset_index()

    Cell_number = []
    Neigboorhood = []
    Neigboorhood_coordinates = []

    for key, value in tqdm(DataFrameDict_Full_cell.items()):

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list()

        count_cell = -1
        for cell in outline_cell:
            count_cell = count_cell + 1
            Cell_number.append(count_cell)

            if len(cell) < 4:
                Neigboorhood.append('Cell size below 4 pixels')
                Neigboorhood_coordinates.append('Cell size below 4 pixels')

            elif len(cell) >=4:
                cell_area = shapely.geometry.Polygon(np.array(cell))

                neighborhood_radius = Radius
                cell_neighborhood = cell_area.buffer(neighborhood_radius)

                cell_neighborhood_coordinates = [tuple(round(value) for value in tup) for tup in tuple(cell_neighborhood.exterior.coords)]
                Neigboorhood_coordinates.append(cell_neighborhood_coordinates)

                Whats_in_the_hood = []
                count_other_cell = -1
                for other_cell in outline_cell:

                    count_other_cell = count_other_cell + 1

                    if len(other_cell) < 4:
                        continue

                    elif len(other_cell) >=4:

                        outline_cell_np = np.array(other_cell)
                        other_cell_outline = shapely.geometry.Polygon(outline_cell_np)

                        if other_cell_outline == cell_area:
                            continue

                        elif cell_neighborhood.intersects(other_cell_outline) == True:
                            Whats_in_the_hood.append(count_other_cell)

                Neigboorhood.append(Whats_in_the_hood)

    Cells_clean['Neighborhood'] = Neigboorhood
    Cells_clean['Cell_number'] = Cell_number
    Cells_clean['Neighborhood_coordinates'] = Neigboorhood_coordinates

    return(Cells_clean)

#calculates the neighboring cells based on distance but does not add any other information
def generate_neighborhood(Cells,Radius):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    Cell_number = []
    Neigboorhood = []
    Neigboorhood_coordinates = []

    for key, value in tqdm(DataFrameDict_Full_cell.items()):

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list()

        count_cell = -1
        for cell in outline_cell:
            cell = ast.literal_eval(cell)
            count_cell = count_cell + 1
            Cell_number.append(count_cell)

            if len(cell) < 4:
                Neigboorhood.append('Cell size below 4 pixels')
                Neigboorhood_coordinates.append('Cell size below 4 pixels')

            elif len(cell) >= 4:

                cell_area = shapely.geometry.Polygon(np.array(cell))

                neighborhood_radius = Radius
                cell_neighborhood = cell_area.buffer(neighborhood_radius)

                cell_neighborhood_coordinates = [tuple(round(value) for value in tup) for tup in tuple(cell_neighborhood.exterior.coords)]
                Neigboorhood_coordinates.append(cell_neighborhood_coordinates)

                Whats_in_the_hood = []
                count_other_cell = -1
                for other_cell in outline_cell:
                    other_cell = ast.literal_eval(other_cell)
                    count_other_cell = count_other_cell + 1

                    if len(other_cell) < 4:
                        continue

                    elif len(other_cell) >= 4:

                        outline_cell_np = np.array(other_cell)
                        other_cell_outline = shapely.geometry.Polygon(outline_cell_np)

                        if other_cell_outline == cell_area:
                            continue

                        elif cell_neighborhood.intersects(other_cell_outline) == True:
                            Whats_in_the_hood.append(count_other_cell)

                Neigboorhood.append(Whats_in_the_hood)

    Cells['Neighborhood'] = Neigboorhood
    Cells['Cell_number'] = Cell_number
    Cells['Neighborhood_coordinates'] = Neigboorhood_coordinates

    return(Cells)