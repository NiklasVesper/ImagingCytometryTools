import pandas as pd
import numpy as np
import shapely.geometry
import warnings
import os
import scandir as sd
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=FutureWarning)

from ImagingCytometryTools.get_data_from_files import get_markers_from_segmentation
from ImagingCytometryTools.get_data_from_files import get_pixel_values

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
def generate_subcellular_localization_basic(Cells, Cytoplasm, Nucleus, FileDirectory, ImageName, outline_cell_image_dir, outline_nucleus_image_dir, nucleus_count = 1):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique()
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique()

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm}
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus}

    if nucleus_count == 1:
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
        for markers in get_markers_from_segmentation(Cells):
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))
            columns.append('Intensity_' + markers + '_Cell')
            columns.append('Intensity_' + markers + '_Cytoplasm')
            columns.append('Intensity_' + markers + '_Nucleus')

    elif nucleus_count == 2:
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
        for markers in get_markers_from_segmentation(Cells):
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))
            columns.append('Intensity_' + markers + '_Cell')
            columns.append('Intensity_' + markers + '_Cytoplasm')
            columns.append('Intensity_' + markers + '_Nucleus_1')
            columns.append('Intensity_' + markers + '_Nucleus_2')

    cells_with_sub_cell_loc = pd.DataFrame(columns=columns)

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()
    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()
    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    file_dir_list_nucleus = []
    for paths, dirs, files in sd.walk(outline_nucleus_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_nucleus.append(filedir)

    for key, value in tqdm(DataFrameDict_Full_cell.items()):

        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list()
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list()
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t)))

        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list()
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list()
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])
        outline_image_nucleus = get_pixel_values(file_dir_list_nucleus[key - 1])

        count_cell = -1
        for cell in position_cells:
            count_cell = count_cell + 1
            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0]))

            if connected_black_pixels_cell == None:
                continue

            edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method='edge')

            if len(edge_pixels_cell) < 4:
                continue

            elif len(edge_pixels_cell) >=4:

                pixel_sorted_cell = pixel_sort(edge_pixels_cell, method='frog_jump')

                if pixel_sorted_cell == None:
                    continue
                
                else:
                    cell_shape = shapely.geometry.Polygon(pixel_sorted_cell)
    
                    count_cytoplasm = -1
                    for cytoplasm in position_cytoplasm:
                        count_cytoplasm = count_cytoplasm + 1
                        cytoplasm_position = shapely.geometry.Point(cytoplasm)
    
                        if cytoplasm_position.within(cell_shape) == True:
    
                            count_nucleus = -1
                            nucleus_list = []
                            for nucleus in position_nucleus:
                                count_nucleus = count_nucleus + 1
                                nucleus_position = shapely.geometry.Point(nucleus)

                                if nucleus_position.within(cell_shape) == True:
                                    nucleus_list.append(count_nucleus)
    
                            if len(nucleus_list) == 1 and nucleus_count == 1:
    
                                connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))

                                if connected_black_pixels_nucleus == None:
                                    continue

                                else:
                                    edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method='edge')
                                    pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method='frog_jump')

                                    if pixel_sorted_nucleus == None:
                                        continue

                                    else:

                                        new_row = {}

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

                                        cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row
    
                            break
    
                            if len(nucleus_list) == 2 and nucleus_count == 2:
    
                                connected_black_pixels_nucleus_1 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                                connected_black_pixels_nucleus_2 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[1]][1]),round(position_nucleus[nucleus_list[1]][0]))

                                if connected_black_pixels_nucleus_1 == None or connected_black_pixels_nucleus_2 == None:
                                    continue

                                else:

                                    edge_pixels_nucleus_1 = find_edge_pixels(connected_black_pixels_nucleus_1, method='edge')
                                    pixel_sorted_nucleus_1 = pixel_sort(edge_pixels_nucleus_1, method='frog_jump')

                                    edge_pixels_nucleus_2 = find_edge_pixels(connected_black_pixels_nucleus_2, method='edge')
                                    pixel_sorted_nucleus_2 = pixel_sort(edge_pixels_nucleus_2, method='frog_jump')

                                    if pixel_sorted_nucleus_1 == None or pixel_sorted_nucleus_2 == None:
                                        continue

                                    else:
                                
                                        new_row = {}

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

                                        cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row
    
                            break

    return(cells_with_sub_cell_loc)

def generate_subcellular_localization_advanced(Cells, Cytoplasm, Nucleus, FileDirectory, ImageName, outline_cell_image_dir, outline_nucleus_image_dir,nucleus_count):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    UniqueNames_Cytoplasm = Cytoplasm.ImageNumber.unique()
    UniqueNames_Nucleus = Nucleus.ImageNumber.unique()

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}
    DataFrameDict_Cytoplasm = {elem: pd.DataFrame() for elem in UniqueNames_Cytoplasm}
    DataFrameDict_Nucleus = {elem: pd.DataFrame() for elem in UniqueNames_Nucleus}

    if nucleus_count == 1:
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
        for markers in get_markers_from_segmentation(Cells):
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))
            columns.append('Intensity_' + markers + '_Cell')
            columns.append('Intensity_' + markers + '_Cytoplasm')
            columns.append('Intensity_' + markers + '_Nucleus')

    elif nucleus_count == 2:
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
        for markers in get_markers_from_segmentation(Cells):
            columns.append('PathName_' + markers.replace('MeanIntensity_', ''))
            columns.append('FileName_' + markers.replace('MeanIntensity_', ''))
            columns.append('Intensity_' + markers + '_Cell')
            columns.append('Intensity_' + markers + '_Cytoplasm')
            columns.append('Intensity_' + markers + '_Nucleus_1')
            columns.append('Intensity_' + markers + '_Nucleus_2')

    cells_with_sub_cell_loc = pd.DataFrame(columns=columns)

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()
    for key in DataFrameDict_Cytoplasm.keys():
        DataFrameDict_Cytoplasm[key] = Cytoplasm[:][Cytoplasm.ImageNumber == key].reset_index()
    for key in DataFrameDict_Nucleus.keys():
        DataFrameDict_Nucleus[key] = Nucleus[:][Nucleus.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    file_dir_list_nucleus = []
    for paths, dirs, files in sd.walk(outline_nucleus_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_nucleus.append(filedir)

    for key, value in tqdm(DataFrameDict_Full_cell.items(), colour='white'):

        x_cell_t = DataFrameDict_Full_cell[key]['Location_Center_X'].to_list()
        y_cell_t = DataFrameDict_Full_cell[key]['Location_Center_Y'].to_list()
        position_cells = np.array(list(zip(x_cell_t, y_cell_t)))

        x_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_X'].to_list()
        y_cytoplasm_t = DataFrameDict_Cytoplasm[key]['Location_Center_Y'].to_list()
        position_cytoplasm = np.array(list(zip(x_cytoplasm_t, y_cytoplasm_t)))

        x_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_X'].to_list()
        y_nucleus_t = DataFrameDict_Nucleus[key]['Location_Center_Y'].to_list()
        position_nucleus = np.array(list(zip(x_nucleus_t, y_nucleus_t)))

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])
        outline_image_nucleus = get_pixel_values(file_dir_list_nucleus[key - 1])

        count_cell = -1
        for cell in position_cells:
            count_cell = count_cell + 1

            connected_black_pixels_cell = get_connected_black_pixels(outline_image_cell, round(cell[1]), round(cell[0]))
            
            if connected_black_pixels_cell == None:
                continue
                
            else:
                edge_pixels_cell = find_edge_pixels(connected_black_pixels_cell, method='edge')
    
                if len(edge_pixels_cell) < 4:
                    continue
    
                elif len(edge_pixels_cell) >=4:
    
                    pixel_sorted_cell = pixel_sort(edge_pixels_cell, method='frog_jump')
                    
                    if pixel_sorted_cell == None:
                        continue
                        
                    else: 
                        cell_shape = shapely.geometry.Polygon(pixel_sorted_cell)
        
                        count_cytoplasm = -1
                        for cytoplasm in position_cytoplasm:
                            count_cytoplasm = count_cytoplasm + 1
                            cytoplasm_position = shapely.geometry.Point(cytoplasm)
        
                            if cytoplasm_position.within(cell_shape) == True:
        
                                count_nucleus = -1
                                nucleus_list = []
                                for nucleus in position_nucleus:
                                    count_nucleus = count_nucleus + 1
        
                                    connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus, round(nucleus[1]),round(nucleus[0]))
                                    
                                    if connected_black_pixels_nucleus == None:
                                        continue
                                        
                                    else:
                                        edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method='edge')
                                        pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method='frog_jump')
                                        
                                        if pixel_sorted_nucleus == None:
                                            continue
                                        
                                        else:
                                            nucleus_shape = shapely.geometry.Polygon(pixel_sorted_nucleus)
                
                                            if nucleus_shape.within(cell_shape) == True:
                                                nucleus_list.append(count_nucleus)
        
                                if len(nucleus_list) == 1 and nucleus_count == 1:
        
                                    connected_black_pixels_nucleus = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                                    edge_pixels_nucleus = find_edge_pixels(connected_black_pixels_nucleus, method='edge')
                                    pixel_sorted_nucleus = pixel_sort(edge_pixels_nucleus, method='frog_jump')
        
                                    new_row = {}
        
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
        
                                    cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row
        
                                break
        
                                if len(nucleus_list) == 2 and nucleus_count == 2:
        
                                    connected_black_pixels_nucleus_1 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[0]][1]),round(position_nucleus[nucleus_list[0]][0]))
                                    edge_pixels_nucleus_1 = find_edge_pixels(connected_black_pixels_nucleus_1, method='edge')
                                    pixel_sorted_nucleus_1 = pixel_sort(edge_pixels_nucleus_1, method='frog_jump')
        
                                    connected_black_pixels_nucleus_2 = get_connected_black_pixels(outline_image_nucleus,round(position_nucleus[nucleus_list[1]][1]),round(position_nucleus[nucleus_list[1]][0]))
                                    edge_pixels_nucleus_2 = find_edge_pixels(connected_black_pixels_nucleus_2, method='edge')
                                    pixel_sorted_nucleus_2 = pixel_sort(edge_pixels_nucleus_2, method='frog_jump')
        
                                    new_row = {}
        
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
        
                                    cells_with_sub_cell_loc.loc[len(cells_with_sub_cell_loc)] = new_row
                                break

    return(cells_with_sub_cell_loc)