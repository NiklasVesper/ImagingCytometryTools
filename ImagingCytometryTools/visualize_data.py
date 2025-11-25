import pandas as pd
import scandir as sd
from PIL import Image, ImageOps
import os
import ast
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import random

from ImagingCytometryTools.get_data_from_files import get_pixel_values
from ImagingCytometryTools.helper_functions import crop_image

'''
A collection of functions to visualise data as images.
'''

#check if the outlines and or subcellular matching is correct
def test_outline_identification_and_matching(Cells, new_folder_name, outline_cell_image_dir, show_individual_cells = False, show_nucleus = False):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    for key, value in DataFrameDict_Full_cell.items():

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])
        outline_image_cell[outline_image_cell > 255] = 255

        filename = os.path.basename(file_dir_list_cell[key - 1])
        filename_string = str(filename)

        dir_name = os.path.dirname(file_dir_list_cell[key - 1])
        dir_name_string = str(dir_name)

        folder_name = os.path.basename(dir_name)
        folder_name_string = str(folder_name)

        filedir_test_image = dir_name_string[:-len(folder_name_string)-17] + 'images/' + new_folder_name

        if os.path.isdir(filedir_test_image) == True:
            pass
        else:
            os.makedirs(filedir_test_image)

        export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list()

        if show_individual_cells == False:

            cell_pixel_area = DataFrameDict_Full_cell[key]['Cell_pixel_area'].to_list()

            if show_nucleus == False:

                for cell_o in outline_cell:
                    cell_o = ast.literal_eval(cell_o)

                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = 200

                for cell_a in cell_pixel_area:
                    cell_a = ast.literal_eval(cell_a)

                    for pixel in cell_a:
                        outline_image_cell[pixel[1], pixel[0]] = 100

            if show_nucleus == True:

                nucleus_pixel_area = DataFrameDict_Full_cell[key]['Nuclear_pixel_area'].to_list()

                for cell_o in outline_cell:
                    cell_o = ast.literal_eval(cell_o)

                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = 200

                for nucleus_a in nucleus_pixel_area:
                    nucleus_a = ast.literal_eval(nucleus_a)

                    for pixel in nucleus_a:
                        outline_image_cell[pixel[1], pixel[0]] = 100

            data = Image.fromarray(outline_image_cell.astype(np.uint8))
            data.save(export_file_path)

        if show_individual_cells == True:

            if show_nucleus == False:

                cell_pixel_area = DataFrameDict_Full_cell[key]['Cell_pixel_area'].to_list()

                cell_count = -1
                for cell_o in outline_cell:
                    color = round(random.uniform(50, 200))

                    cell_count = cell_count + 1

                    cell_o = ast.literal_eval(cell_o)
                    cell_a = ast.literal_eval(cell_pixel_area[cell_count])
                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = color

                    for pixel in cell_a:
                        outline_image_cell[pixel[1], pixel[0]] = color
            
            if show_nucleus == True:

                nucleus_pixel_area = DataFrameDict_Full_cell[key]['Nuclear_pixel_area'].to_list()

                cell_count = -1
                for cell_o in outline_cell:

                    color = round(random.uniform(50, 200))
                    cell_count = cell_count + 1
                    cell_o = ast.literal_eval(cell_o)
                    nucleus_a = ast.literal_eval(nucleus_pixel_area[cell_count])

                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = color

                    for pixel in nucleus_a:
                        outline_image_cell[pixel[1], pixel[0]] = color

            image = Image.fromarray(outline_image_cell.astype(np.uint8))
            image.save(export_file_path)

#generates image galaries
def generate_image_galaries(Cells, df_column, output_folder_phenotypes, outline_cell_image_dir, proteins, max_pixel, contrast_multiplier = [1,1,1], select_cell_type_and_state = [False], add_mixed_cells = True, generate_crops = True, crop_size = 40, show_neighborhood_radius = True, select_neighboring_cell_type_and_state=[False], show_neighboring_cells = False):

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    filename = os.path.basename(Cells)
    filename_string = str(filename)

    dir_name = os.path.dirname(Cells)
    dir_name_string = str(dir_name)

    folder_name = os.path.basename(dir_name)
    folder_name_string = str(folder_name)

    Cells = pd.read_csv(Cells)

    if generate_crops == False:

        if len(proteins) == 1:

            UniqueImageNumber = Cells.ImageNumber.unique()
            
            DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

            for key in DataFrameDict_Full_cell.keys():
                DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

            path_name_1 = 'PathName_' + proteins[0]
            file_name_1 = 'FileName_' + proteins[0]

            UniquePathName_1 = Cells[path_name_1].unique()
            UniqueFileName_1 = Cells[file_name_1].unique()

            for number in UniqueImageNumber:
                
                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_red[outline_image_cell_red > 255] = 255
                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_green[outline_image_cell_green > 255] = 255
                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))

                image_normalized_1 = np.divide(image_1, max_pixel[0])

                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])

                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)

                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255

                Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                if select_cell_type_and_state[0] == True:

                    for index, cell in Cells_on_image.iterrows():

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        outline_image_cell_green[pixel[1], pixel[0]] = 0

                                if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):
                                            
                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                    
                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_green[pixel[1], pixel[0]] = 255
                                                        
                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_green[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255
                                        
                        if add_mixed_cells == False:
                            
                            if len(cell_types_and_states) == 2:

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        outline_image_cell_green[pixel[1], pixel[0]] = 0

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_green[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_green[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255

                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = outline_image_cell_green.astype(np.uint8)
                image_normalized_b = outline_image_cell_blue.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(number) + '.png'

                image.save(export_file_path_full)
                
        if len(proteins) == 2:

            UniqueImageNumber = Cells.ImageNumber.unique()

            DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

            for key in DataFrameDict_Full_cell.keys():
                DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

            path_name_1 = 'PathName_' + proteins[0]
            file_name_1 = 'FileName_' + proteins[0]
            path_name_2 = 'PathName_' + proteins[1]
            file_name_2 = 'FileName_' + proteins[1]

            UniquePathName_1 = Cells[path_name_1].unique()
            UniqueFileName_1 = Cells[file_name_1].unique()
            UniquePathName_2 = Cells[path_name_2].unique()
            UniqueFileName_2 = Cells[file_name_2].unique()

            for number in UniqueImageNumber:

                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_red[outline_image_cell_red > 255] = 255
                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_green[outline_image_cell_green > 255] = 255
                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                image_normalized_1 = np.divide(image_1, max_pixel[0])
                image_normalized_2 = np.divide(image_2, max_pixel[1])

                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)

                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                if select_cell_type_and_state[0] == True:

                    for index, cell in Cells_on_image.iterrows():

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 0

                                if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255

                        if add_mixed_cells == False:

                            if len(cell_types_and_states) == 2:

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 0

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255

                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                image_normalized_b = outline_image_cell_blue.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(proteins[1]) + '_' + str(number) + '.png'

                image.save(export_file_path_full)

        if len(proteins) == 3:

            UniqueImageNumber = Cells.ImageNumber.unique()

            DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

            for key in DataFrameDict_Full_cell.keys():
                DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

            path_name_1 = 'PathName_' + proteins[0]
            file_name_1 = 'FileName_' + proteins[0]
            path_name_2 = 'PathName_' + proteins[1]
            file_name_2 = 'FileName_' + proteins[1]
            path_name_3 = 'PathName_' + proteins[2]
            file_name_3 = 'FileName_' + proteins[2]

            UniquePathName_1 = Cells[path_name_1].unique()
            UniqueFileName_1 = Cells[file_name_1].unique()
            UniquePathName_2 = Cells[path_name_2].unique()
            UniqueFileName_2 = Cells[file_name_2].unique()
            UniquePathName_3 = Cells[path_name_3].unique()
            UniqueFileName_3 = Cells[file_name_3].unique()

            for number in UniqueImageNumber:

                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_red[outline_image_cell_red > 255] = 255
                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_green[outline_image_cell_green > 255] = 255
                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                image_normalized_1 = np.divide(image_1, max_pixel[0])
                image_normalized_2 = np.divide(image_2, max_pixel[1])
                image_normalized_3 = np.divide(image_3, max_pixel[2])

                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)
                image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell_blue)

                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                if select_cell_type_and_state[0] == True:

                    for index, cell in Cells_on_image.iterrows():

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 0

                                if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255

                        if add_mixed_cells == False:

                            if len(cell_types_and_states) == 2:

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0]:

                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_1[pixel[1], pixel[0]] = 0

                                    for pixel in outline_individuall_cell:
                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 0

                                if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                    outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])
                                                    outline_individuall_neighboring_cell = ast.literal_eval(Cells_on_image.iloc[neighboring_cell]['Cell_outline'])

                                                    for pixel in outline_individuall_cell:
                                                        image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 255

                                                    for pixel in outline_individuall_neighboring_cell:
                                                        outline_image_cell_blue[pixel[1], pixel[0]] = 0

                                    if select_neighboring_cell_type_and_state[0] == False:

                                        outline_individuall_cell = ast.literal_eval(cell['Cell_outline'])

                                        for pixel in outline_individuall_cell:
                                            image_normalized_contrasted_add_outline_2[pixel[1], pixel[0]] = 255

                                        for pixel in outline_individuall_cell:
                                            outline_image_cell_blue[pixel[1], pixel[0]] = 255

                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(proteins[1]) + '_' + str(proteins[2]) + '_' + str(number) + '.png'

                image.save(export_file_path_full)

    if generate_crops == True:

        output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes

        if os.path.isdir(output_folder) == True:
            pass
        else:
            os.makedirs(output_folder)

        try:
            for file in os.listdir(output_folder):
                if file.endswith('.png'):
                    os.remove(str(output_folder) + '/' + str(file))

        finally:

            image_counter = 0

            if len(proteins) == 1:

                UniqueImageNumber = Cells.ImageNumber.unique()

                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                path_name_1 = 'PathName_' + proteins[0]
                file_name_1 = 'FileName_' + proteins[0]

                UniquePathName_1 = Cells[path_name_1].unique()
                UniqueFileName_1 = Cells[file_name_1].unique()

                for number in UniqueImageNumber:

                    image_counter = image_counter + 1

                    Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                    if select_cell_type_and_state[0] == True:

                        cell_counter = 0

                        for index, cell in Cells_on_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))

                                                image_normalized_1 = np.divide(image_1, max_pixel[0])

                                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])

                                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)

                                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = outline_image_cell_green.astype(np.uint8)
                                                    image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(outline_image_cell_green, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)  # Save to a file

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                            else:
                                                continue

                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                            image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))

                                            image_normalized_1 = np.divide(image_1, max_pixel[0])

                                            image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])

                                            image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)

                                            image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255

                                            if show_neighborhood_radius == False:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                image_normalized_g = outline_image_cell_green.astype(np.uint8)
                                                image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                            if show_neighborhood_radius == True:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                image_g = np.add(outline_image_cell_green, outline_image_neighborhood)
                                                image_g[image_g > 255] = 255
                                                image_normalized_g = image_g.astype(np.uint8)

                                                image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                image_b[image_b > 255] = 255
                                                image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y,crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                            croped_final_image.save(export_file_path_full)

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))

                                                image_normalized_1 = np.divide(image_1, max_pixel[0])

                                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])

                                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)

                                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = outline_image_cell_green.astype(np.uint8)
                                                    image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(outline_image_cell_green, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)  # Save to a file

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)

                                            else:
                                                continue

                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                            image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))

                                            image_normalized_1 = np.divide(image_1, max_pixel[0])

                                            image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])

                                            image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)

                                            image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255

                                            if show_neighborhood_radius == False:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                image_normalized_g = outline_image_cell_green.astype(np.uint8)
                                                image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                            if show_neighborhood_radius == True:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                image_g = np.add(outline_image_cell_green, outline_image_neighborhood)
                                                image_g[image_g > 255] = 255
                                                image_normalized_g = image_g.astype(np.uint8)

                                                image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                image_b[image_b > 255] = 255
                                                image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                            croped_final_image.save(export_file_path_full)

            if len(proteins) == 2:

                UniqueImageNumber = Cells.ImageNumber.unique()

                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                path_name_1 = 'PathName_' + proteins[0]
                file_name_1 = 'FileName_' + proteins[0]
                path_name_2 = 'PathName_' + proteins[1]
                file_name_2 = 'FileName_' + proteins[1]

                UniquePathName_1 = Cells[path_name_1].unique()
                UniqueFileName_1 = Cells[file_name_1].unique()
                UniquePathName_2 = Cells[path_name_2].unique()
                UniqueFileName_2 = Cells[file_name_2].unique()

                for number in UniqueImageNumber:

                    image_counter = image_counter + 1

                    Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                    if select_cell_type_and_state[0] == True:

                        cell_counter = 0

                        for index, cell in Cells_on_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                                                image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                image_normalized_2 = np.divide(image_2, max_pixel[1])

                                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)

                                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                    image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)

                                            else:
                                                continue

                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                    image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                                                    if show_neighborhood_radius == False:
                                                        
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                        image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                    if show_neighborhood_radius == True:
                                                        
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y,crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'  # creates export file path
                                            croped_final_image.save(export_file_path_full)


                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                        image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))  # imports the proper image
                                                        image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                                                        image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                        image_normalized_2 = np.divide(image_2, max_pixel[1])

                                                        image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                        image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                                                        image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                        image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)

                                                        image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                        image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                                                        if show_neighborhood_radius == False:
                                                            image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                            image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                            image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                        if show_neighborhood_radius == True:
                                                            image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                            image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                            image_g[image_g > 255] = 255
                                                            image_normalized_g = image_g.astype(np.uint8)

                                                            image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                            image_b[image_b > 255] = 255
                                                            image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)  # Save to a file

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)

                                            else:
                                                continue


                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            # imports the cellular outlines
                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))  # imports the proper image
                                                    image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                                                    if show_neighborhood_radius == False:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                        image_normalized_b = outline_image_cell_blue.astype(np.uint8)

                                                    if show_neighborhood_radius == True:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(outline_image_cell_blue, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'  # creates export file path
                                            croped_final_image.save(export_file_path_full)  # Save to a file

            if len(proteins) == 3:

                UniqueImageNumber = Cells.ImageNumber.unique()

                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                path_name_1 = 'PathName_' + proteins[0]
                file_name_1 = 'FileName_' + proteins[0]
                path_name_2 = 'PathName_' + proteins[1]
                file_name_2 = 'FileName_' + proteins[1]
                path_name_3 = 'PathName_' + proteins[2]
                file_name_3 = 'FileName_' + proteins[2]

                UniquePathName_1 = Cells[path_name_1].unique()
                UniqueFileName_1 = Cells[file_name_1].unique()
                UniquePathName_2 = Cells[path_name_2].unique()
                UniqueFileName_2 = Cells[file_name_2].unique()
                UniquePathName_3 = Cells[path_name_3].unique()
                UniqueFileName_3 = Cells[file_name_3].unique()

                for number in UniqueImageNumber:

                    image_counter = image_counter + 1

                    Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                    if select_cell_type_and_state[0] == True:

                        cell_counter = 0

                        for index, cell in Cells_on_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                        image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                        image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                                                        image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                                                        image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                        image_normalized_2 = np.divide(image_2, max_pixel[1])
                                                        image_normalized_3 = np.divide(image_3, max_pixel[2])

                                                        image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                        image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                                                        image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                                                        image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                        image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)
                                                        image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell_blue)

                                                        image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                        image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                                                        image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                    image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)  # Save to a file

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)

                                            else:
                                                continue

                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                    image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                                                    image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])
                                                    image_normalized_3 = np.divide(image_3, max_pixel[2])

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                                                    image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)
                                                    image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell_blue)

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                                                    image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                    image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                            croped_final_image.save(export_file_path_full)

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            is_there_a_propper_neighbor = 0
                                            neighbor_numbers = []
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        is_there_a_propper_neighbor = is_there_a_propper_neighbor + 1
                                                        neighbor_numbers.append(neighboring_cell)

                                            if is_there_a_propper_neighbor > 0:

                                                cell_counter = cell_counter + 1

                                                postition_x = round(cell['Location_Center_X'])
                                                postition_y = round(cell['Location_Center_Y'])

                                                outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_red[outline_image_cell_red > 255] = 255
                                                outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_green[outline_image_cell_green > 255] = 255
                                                outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                                if show_neighborhood_radius == True:

                                                    outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])

                                                    for pixel in outline_image_neighborhood_cell:

                                                        pixel = list(pixel)

                                                        if pixel[0] < 0:
                                                            pixel[0] = 0
                                                        if pixel[1] < 0:
                                                            pixel[1] = 0
                                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[0] = (sub_list[1] - 2)
                                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                            sub_list = list(outline_image_neighborhood.shape)
                                                            pixel[1] = (sub_list[0] - 2)

                                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                        image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                        image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                                                        image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                                                        image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                        image_normalized_2 = np.divide(image_2, max_pixel[1])
                                                        image_normalized_3 = np.divide(image_3, max_pixel[2])

                                                        image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                        image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                                                        image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                                                        image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                        image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)
                                                        image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell_blue)

                                                        image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                        image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                                                        image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                                                    if show_neighborhood_radius == False:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                        image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)

                                                    if show_neighborhood_radius == True:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                image = Image.fromarray(rgb_array, mode='RGB')
                                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                                croped_final_image.save(export_file_path_full)  # Save to a file

                                                if show_neighboring_cells == True:

                                                    outline_image_base = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[number - 1])

                                                    cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                                                    for pixel in cell_area:
                                                        outline_image_cell[pixel[1], pixel[0]] = 255

                                                    for neighbor_number in neighbor_numbers:

                                                        neighboring_cell_area = ast.literal_eval(Cells_on_image.iloc[neighbor_number]['Cell_pixel_area'])

                                                        for pixel in neighboring_cell_area:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                    rgb_array = np.stack((outline_image_cell.astype(np.uint8), outline_image_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'

                                                    croped_final_image_neigboorhood.save(export_file_path_full)

                                            else:
                                                continue

                                        if select_neighboring_cell_type_and_state[0] == False:

                                            cell_counter = cell_counter + 1

                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            outline_image_cell_red = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_red[outline_image_cell_red > 255] = 255
                                            outline_image_cell_green = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_green[outline_image_cell_green > 255] = 255
                                            outline_image_cell_blue = get_pixel_values(file_dir_list_cell[number - 1])
                                            outline_image_cell_blue[outline_image_cell_blue > 255] = 255

                                            if show_neighborhood_radius == True:

                                                outline_image_neighborhood = get_pixel_values(file_dir_list_cell[number - 1])
                                                outline_image_neighborhood[outline_image_neighborhood > 0] = 0

                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_neighborhood.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                                                    image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                                                    image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])
                                                    image_normalized_3 = np.divide(image_3, max_pixel[2])

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                                                    image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell_red)
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell_green)
                                                    image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell_blue)

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                                                    image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                                                if show_neighborhood_radius == False:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                    image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                                                    image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)

                                                if show_neighborhood_radius == True:
                                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                    image_g = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                    image_g[image_g > 255] = 255
                                                    image_normalized_g = image_g.astype(np.uint8)

                                                    image_b = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                    image_b[image_b > 255] = 255
                                                    image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(image_counter) + '_' + str(cell_counter) + '.png'
                                            croped_final_image.save(export_file_path_full)

#tests the ISO
def test_image_analysis_of_spatial_overlap(Cells, outline_cell_image_dir, output_folder, dir_name_string, folder_name, crop_size = 30):

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder

    if os.path.isdir(output_folder) == True:
        pass
    else:
        os.makedirs(output_folder)

    try:
        for file in os.listdir(output_folder):
            if file.endswith('.png'):
                os.remove(str(output_folder) + '/' + str(file))

    finally:

        image_counter = -1
        for index, cell in Cells.iterrows():

            if isinstance(cell['ISO_visualisation_and_information'],str):

                image_counter = image_counter + 1

                cell_type_ISO = ast.literal_eval(cell['cell_types_and_states_ISO'])

                postition_x = round(cell['Location_Center_X'])
                postition_y = round(cell['Location_Center_Y'])

                ISO_visualisation = ast.literal_eval(cell['ISO_visualisation_and_information'])

                image_area = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                cell_area = ast.literal_eval(cell['Cell_pixel_area'])

                for pixel in cell_area:
                    image_area[pixel[1], pixel[0]] = 170
                image_area = Image.fromarray(image_area.astype(np.uint8))
                croped_image_area = crop_image(image_area, postition_x, postition_y, crop_size)

                protein = str(cell_type_ISO[0][0])[:-1]
                path_name = 'PathName_' + protein
                file_name = 'FileName_' + protein

                ppas = ast.literal_eval(cell['pixel_positive_area_'+protein])
                image_ppa = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                image_ppa[image_ppa > 0] = 0
                for ppa in ppas:
                    for pixel in ppa:
                        image_ppa[pixel[0], pixel[1]] = 200
                image_ppa = Image.fromarray(image_ppa.astype(np.uint8))

                image_cell_type_ISO = get_pixel_values(cell[path_name] + '/' + cell[file_name])
                image_cell_type_ISO = Image.fromarray(image_cell_type_ISO.astype(np.uint8))
                croped_image_cell_type_ISO = crop_image(image_cell_type_ISO, postition_x, postition_y, crop_size)
                croped_image_ppa_og = crop_image(image_ppa, postition_x, postition_y, crop_size)

                fig, axs = plt.subplots(nrows=len(ISO_visualisation)+1, ncols=3, figsize=(17.5, 10))

                axs[0, 0].imshow(croped_image_area, cmap='gray')
                axs[0, 0].set_title(cell_type_ISO[0][0])

                axs[0, 1].imshow(croped_image_cell_type_ISO, cmap='gray')
                axs[0, 1].set_title(cell_type_ISO[0][0])

                axs[0, 2].imshow(croped_image_ppa_og, cmap='magma')
                axs[0, 2].set_title(cell_type_ISO[0][0])

                light_gray_patch = mpatches.Patch(color='#CCCCCC', label='Overlap corrected cell')
                dark_gray_patch = mpatches.Patch(color='#666666', label='Overlaping cell')
                fig.legend(handles=[light_gray_patch, dark_gray_patch], loc="upper right")

                for x in range(0, len(ISO_visualisation)):

                    protein = ISO_visualisation[x][0]
                    path_name = 'PathName_' + protein
                    file_name = 'FileName_' + protein
                    image_protein = get_pixel_values(cell[path_name] + '/' + cell[file_name])
                    outline_image_protein = Image.fromarray(image_protein.astype(np.uint8))
                    croped_image_protein = crop_image(outline_image_protein, postition_x, postition_y, crop_size)


                    index_cell = ISO_visualisation[x][1]
                    other_cell_area = ast.literal_eval(Cells.loc[index_cell, 'Cell_pixel_area'])
                    image_area = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                    for pixel in other_cell_area:
                        image_area[pixel[1], pixel[0]] = 100
                    for pixel in cell_area:
                        image_area[pixel[1], pixel[0]] = 170

                    image_area = Image.fromarray(image_area.astype(np.uint8))
                    croped_image_area= crop_image(image_area, postition_x, postition_y, crop_size)

                    ppa = ISO_visualisation[x][2]
                    image_ppa = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                    image_ppa[image_ppa > 0] = 0
                    for pixel in ppa:
                        image_ppa[pixel[1], pixel[0]] = 200
                    image_ppa = Image.fromarray(image_ppa.astype(np.uint8))
                    croped_image_ppa = crop_image(image_ppa, postition_x, postition_y, crop_size)

                    axs[x + 1, 0].imshow(croped_image_area, cmap = 'gray')
                    axs[x + 1, 0].set_title(protein)

                    axs[x+1, 1].imshow(croped_image_protein, cmap = 'gray')
                    axs[x+1, 1].set_title(protein)

                    axs[x+1, 2].imshow(croped_image_ppa, cmap = 'magma')
                    axs[x+1, 2].set_title(protein)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder + '/' + 'ISO_test_' + str(image_counter) + '.png' #generates the output file
                fig.savefig(export_file_path_full, bbox_inches='tight', dpi=300)
                plt.close()

#shows the ISO
def show_image_analysis_of_spatial_overlap(Cells, new_folder_name, outline_cell_image_dir):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    for key, value in DataFrameDict_Full_cell.items():

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])
        outline_image_cell[outline_image_cell > 255] = 255

        filename = os.path.basename(file_dir_list_cell[key - 1])
        filename_string = str(filename)

        dir_name = os.path.dirname(file_dir_list_cell[key - 1])
        dir_name_string = str(dir_name)

        folder_name = os.path.basename(dir_name)
        folder_name_string = str(folder_name)

        filedir_test_image = dir_name_string[:-len(folder_name_string) - 17] + 'images/' + new_folder_name

        if os.path.isdir(filedir_test_image) == True:
            pass
        else:
            os.makedirs(filedir_test_image)

        export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'

        for index, cell in Cells.iterrows():

            if int(cell['ImageNumber']) == key:

                cell_type = ast.literal_eval(cell['cell_types_and_states'])
                cell_type_ISO = ast.literal_eval(cell['cell_types_and_states_ISO'])

                cell_pixel_area = ast.literal_eval(cell['Cell_pixel_area'])

                for pixel in cell_pixel_area:
                    outline_image_cell[pixel[1], pixel[0]] = 200

                if len(cell_type) > 2:

                    for pixel in cell_pixel_area:
                        outline_image_cell[pixel[1], pixel[0]] = 120

                if len(cell_type) > 2 and len(cell_type_ISO) == 2:

                    pixel_counter = 0
                    for pixel in cell_pixel_area:

                        pixel_counter = pixel_counter + 1
                        if pixel_counter % 2 == 0:
                            outline_image_cell[pixel[1], pixel[0]] = 80
                        else:
                            outline_image_cell[pixel[1], pixel[0]] = 120

        data = Image.fromarray(outline_image_cell.astype(np.uint8))
        data.save(export_file_path)

#test the area manipulation
def area_manipulation_test(Cells, new_folder_name, outline_cell_image_dir):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}

    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    file_dir_list_cell = []
    for paths, dirs, files in sd.walk(outline_cell_image_dir):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):
                file_dir_list_cell.append(filedir)

    for key, value in DataFrameDict_Full_cell.items():
        
        Cells_on_image = pd.DataFrame(value)

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])
        outline_image_cell[outline_image_cell > 255] = 255

        filename = os.path.basename(file_dir_list_cell[key - 1])
        filename_string = str(filename)

        dir_name = os.path.dirname(file_dir_list_cell[key - 1])
        dir_name_string = str(dir_name)

        folder_name = os.path.basename(dir_name)
        folder_name_string = str(folder_name)

        filedir_test_image = dir_name_string[:-len(folder_name_string) - 17] + 'images/' + new_folder_name

        if os.path.isdir(filedir_test_image) == True:
            pass
        else:
            os.makedirs(filedir_test_image)

        export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'

        for index, cell in Cells_on_image.iterrows():
            
            if cell['Selected Areas'] == 'in the area':
        
                cell_pixel_area = ast.literal_eval(cell['Cell_pixel_area'])
        
                for pixel in cell_pixel_area:
                    outline_image_cell[pixel[1], pixel[0]] = 100

        data = Image.fromarray(outline_image_cell.astype(np.uint8))
        data.save(export_file_path)