import pandas as pd
import scandir as sd
from PIL import Image, ImageOps
import os
import ast
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import random

from ImagingCytometryTools.get_stuff import get_pixel_values

'''
A collection of functions to visualise data as images.
'''

#for easier image copping
def crop_image(image, pos_x, pos_y, crop_size):

    left = pos_x - crop_size
    right = pos_x + crop_size

    upper = pos_y + crop_size
    lower = pos_y - crop_size

    image_croped = image.crop(((left, lower, right, upper)))

    return(image_croped)

#check if the outlines and or subcellular matching is correct
def test_outline_identification_and_matching(Cells, new_folder_name, outline_cell_image_dir, show_individual_cells = False, show_nucleus = False):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()  # creates a list of all unique images

    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  # creates a dictionary of all unique images

    # resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    file_dir_list_cell = []  # list for the cell image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".tiff"):  # returns all files that end with .tiff
                file_dir_list_cell.append(filedir)  # appends the found .tiff files into the list

    #goes throw each image
    for key, value in DataFrameDict_Full_cell.items():

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])  # gets the cell outline image associated with the image key
        outline_image_cell[outline_image_cell > 255] = 255

        # gives you the file name
        filename = os.path.basename(file_dir_list_cell[key - 1])
        filename_string = str(filename)

        # gives you the directory name
        dir_name = os.path.dirname(file_dir_list_cell[key - 1])
        dir_name_string = str(dir_name)

        # gives you the folder name
        folder_name = os.path.basename(dir_name)
        folder_name_string = str(folder_name)

        filedir_test_image = dir_name_string[:-len(folder_name_string)-17] + 'images/' + new_folder_name  # creates a new directory for the neighborhood analysis

        # checks if the new directory already exists
        if os.path.isdir(filedir_test_image) == True:
            pass
        else:
            os.makedirs(filedir_test_image)  # creates a new directory if necessary

        export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'  # creates a export file path

        outline_cell = DataFrameDict_Full_cell[key]['Cell_outline'].to_list() #makes a list out of the cell outlines

        if show_individual_cells == False:

            cell_pixel_area = DataFrameDict_Full_cell[key]['Cell_pixel_area'].to_list()  # makes a list out of the cell outlines

            if show_nucleus == False:

                #assigns the cell outline a gray scale value
                for cell_o in outline_cell:
                    cell_o = ast.literal_eval(cell_o)

                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = 200

                # assigns the cell area a gray scale value
                for cell_a in cell_pixel_area:
                    cell_a = ast.literal_eval(cell_a)

                    for pixel in cell_a:
                        outline_image_cell[pixel[1], pixel[0]] = 100

            if show_nucleus == True:

                nucleus_pixel_area = DataFrameDict_Full_cell[key]['Nuclear_pixel_area'].to_list()  # makes a list out of the cell outlines

                # assigns the cell outline a gray scale value
                for cell_o in outline_cell:
                    cell_o = ast.literal_eval(cell_o)

                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = 200

                # assigns the nucleus area a gray scale value
                for nucleus_a in nucleus_pixel_area:
                    nucleus_a = ast.literal_eval(nucleus_a)

                    for pixel in nucleus_a:
                        outline_image_cell[pixel[1], pixel[0]] = 100

            data = Image.fromarray(outline_image_cell.astype(np.uint8))
            data.save(export_file_path)

        if show_individual_cells == True:

            if show_nucleus == False:

                cell_pixel_area = DataFrameDict_Full_cell[key]['Cell_pixel_area'].to_list()  # makes a list out of the cell outlines

                cell_count = -1
                #assigns each individual cell outline and area a gray scale value
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

                nucleus_pixel_area = DataFrameDict_Full_cell[key]['Nuclear_pixel_area'].to_list()  # makes a list out of the cell outlines

                cell_count = -1
                # assigns each individual cell outline and nuclear area a gray scale value
                for cell_o in outline_cell:
                    color = round(random.uniform(50, 200))

                    cell_count = cell_count + 1

                    cell_o = ast.literal_eval(cell_o)
                    nucleus_a = ast.literal_eval(nucleus_pixel_area[cell_count])
                    for pixel in cell_o:
                        outline_image_cell[pixel[1], pixel[0]] = color

                    for pixel in nucleus_a:
                        outline_image_cell[pixel[1], pixel[0]] = color

            #saves the image
            image = Image.fromarray(outline_image_cell.astype(np.uint8))
            image.save(export_file_path)

#generates image galaries
def generate_image_galaries(Cells, output_folder_phenotypes, outline_cell_image_dir, proteins, cell_types_and_state, max_pixel, contrast_multiplier = [1,1,1], crop_size = 40, add_mixed_cells = True, show_neighborhood = True, select_neighbors=[False], show_neighboring_cells = False, generate_crops = True):

    file_dir_list_cell = []  # list for the image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".tiff"):  # returns all files that end with .tiff
                file_dir_list_cell.append(filedir)  # appends the found .tiff files into the list

    # gives you the file name
    filename = os.path.basename(Cells)
    filename_string = str(filename)

    # gives you the directory name
    dir_name = os.path.dirname(Cells)
    dir_name_string = str(dir_name)

    # gives you the folder name
    folder_name = os.path.basename(dir_name)
    folder_name_string = str(folder_name)

    Cells = pd.read_csv(Cells) #imports the data frame

    #generates the full image
    if generate_crops == False:

        if len(proteins) == 1:

            UniqueImageNumber = Cells.ImageNumber.unique() #gets the image numbers

            path_name_1 = 'PathName_' + proteins[0] #creates the path name for the protein
            file_name_1 = 'FileName_' + proteins[0] #creates the file name for the protein

            UniquePathName_1 = Cells[path_name_1].unique() #gets the unique path names
            UniqueFileName_1 = Cells[file_name_1].unique() #gets the unique file names

            #goes throw each image
            for number in UniqueImageNumber:

                #imports the cellular outlines
                outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell[outline_image_cell > 255] = 255

                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1])) #imports the proper image

                image_normalized_1 = np.divide(image_1, max_pixel[0]) #normalizes the image to the choosen max pixel value

                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0]) #creates the contrast

                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell) #adds the outlines

                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255 #sets everything above 255 to 255

                #generates the rgb image
                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = outline_image_cell.astype(np.uint8)
                image_normalized_b = outline_image_cell.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes #generates the output folder

                #checks if the folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(number) + '.png' #generates the output file

                image.save(export_file_path_full) #exports the image

        if len(proteins) == 2:

            UniqueImageNumber = Cells.ImageNumber.unique() #gets the image numbers

            path_name_1 = 'PathName_' + proteins[0] #creates the path name for the protein
            file_name_1 = 'FileName_' + proteins[0] #creates the file name for the protein
            path_name_2 = 'PathName_' + proteins[1] #creates the path name for the protein
            file_name_2 = 'FileName_' + proteins[1] #creates the file name for the protein

            UniquePathName_1 = Cells[path_name_1].unique() #gets the unique path names
            UniqueFileName_1 = Cells[file_name_1].unique() #gets the unique file names
            UniquePathName_2 = Cells[path_name_2].unique() #gets the unique path names
            UniqueFileName_2 = Cells[file_name_2].unique() #gets the unique file names

            #goes throw each image
            for number in UniqueImageNumber:

                # imports the cellular outlines
                outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell[outline_image_cell > 255] = 255

                # imports the proper image
                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))

                # normalizes the image to the choosen max pixel value
                image_normalized_1 = np.divide(image_1, max_pixel[0])
                image_normalized_2 = np.divide(image_2, max_pixel[1])

                # creates the contrast
                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])

                # adds the outlines
                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)
                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)

                # sets everything above 255 to 255
                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255

                # generates the rgb image
                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = outline_image_cell.astype(np.uint8)
                image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes #generates the output folder

                # checks if the folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(proteins[1]) + '_' + str(number) + '.png' #generates the output file

                image.save(export_file_path_full) #exports the image

        if len(proteins) == 3:

            UniqueImageNumber = Cells.ImageNumber.unique() #gets the image numbers

            path_name_1 = 'PathName_' + proteins[0] #creates the path name for the protein
            file_name_1 = 'FileName_' + proteins[0] #creates the file name for the protein
            path_name_2 = 'PathName_' + proteins[1] #creates the path name for the protein
            file_name_2 = 'FileName_' + proteins[1] #creates the file name for the protein
            path_name_3 = 'PathName_' + proteins[2] #creates the path name for the protein
            file_name_3 = 'FileName_' + proteins[2] #creates the file name for the protein

            UniquePathName_1 = Cells[path_name_1].unique() #gets the unique path names
            UniqueFileName_1 = Cells[file_name_1].unique() #gets the unique file names
            UniquePathName_2 = Cells[path_name_2].unique() #gets the unique path names
            UniqueFileName_2 = Cells[file_name_2].unique() #gets the unique file names
            UniquePathName_3 = Cells[path_name_3].unique() #gets the unique path names
            UniqueFileName_3 = Cells[file_name_3].unique() #gets the unique file names

            # goes throw each image
            for number in UniqueImageNumber:

                # imports the cellular outlines
                outline_image_cell = get_pixel_values(file_dir_list_cell[number - 1])
                outline_image_cell[outline_image_cell > 255] = 255

                # imports the proper image
                image_1 = get_pixel_values(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
                image_2 = get_pixel_values(str(UniquePathName_2[number - 1]) + '/' + str(UniqueFileName_2[number - 1]))
                image_3 = get_pixel_values(str(UniquePathName_3[number - 1]) + '/' + str(UniqueFileName_3[number - 1]))

                # normalizes the image to the choosen max pixel value
                image_normalized_1 = np.divide(image_1, max_pixel[0])
                image_normalized_2 = np.divide(image_2, max_pixel[1])
                image_normalized_3 = np.divide(image_3, max_pixel[2])

                # creates the contrast
                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])
                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])
                image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])

                # adds the outlines
                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)
                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)
                image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell)

                # sets everything above 255 to 255
                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255
                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255
                image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255

                # generates the rgb image
                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                image_normalized_g = image_normalized_contrasted_add_outline_2.astype(np.uint8)
                image_normalized_b = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                image = Image.fromarray(rgb_array, mode='RGB')

                output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes #generates the output folder

                # checks if the folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(proteins[1]) + '_' + str(proteins[2]) + '_' + str(number) + '.png' #generates the output file

                image.save(export_file_path_full) #exports the image

    # generates croped images of individual cells
    if generate_crops == True:

        output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes  # generates output folder

        # checks if the folder exists
        if os.path.isdir(output_folder) == True:
            pass
        else:
            os.makedirs(output_folder)

        try:
            for file in os.listdir(output_folder):
                if file.endswith('.png'):
                    os.remove(str(output_folder) + '/' + str(file))

        finally:
            image_counter = 0 #counter to name the images

            # loops through all cells
            for index, cell in Cells.iterrows():

                cell_type_and_states = ast.literal_eval(cell['cell_types_and_states']) #imports the assigned cell type

                #includes all cells that have that type
                if add_mixed_cells == True:

                    #goes throw all assigned cell types
                    for assigned_cell_type_and_state in cell_type_and_states:

                        if select_neighbors[0] == False:

                            #checks if one of the assigned cell types matches the wnated type
                            if assigned_cell_type_and_state[0] == cell_types_and_state[0] and set(cell_types_and_state[1:]).issubset(assigned_cell_type_and_state[1:]):

                                image_counter = image_counter + 1 #ticks up the counter for the image name

                                if len(proteins) == 1:

                                    #gets the position of the cell to crop the image
                                    postition_x = round(cell['Location_Center_X'])
                                    postition_y = round(cell['Location_Center_Y'])

                                    #gets the outline
                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell[outline_image_cell > 255] = 255

                                    #adds the neighborhood radius
                                    if show_neighborhood == True:

                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0) #creates a blank image
                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                        #adds the neighborhood coordinates to the blanck image
                                        for pixel in outline_image_neighborhood_cell:

                                            pixel = list(pixel)

                                            if pixel[0] < 0:
                                                pixel[0] = 0
                                            if pixel[1] < 0:
                                                pixel[1] = 0
                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[0] = (sub_list[1] - 2)
                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[1] = (sub_list[0] - 2)

                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                    path_name_1 = 'PathName_' + proteins[0]
                                    file_name_1 = 'FileName_' + proteins[0]

                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])

                                    image_1 = get_pixel_values(path_1)

                                    image_normalized_1 = np.divide(image_1, max_pixel[0]) #normalizes the image to the choosen max pixel value

                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0]) #creates the contrast

                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell) #adds the outlines

                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255 #sets everything above 255 to 255

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == False:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                        image_normalized_g = outline_image_cell.astype(np.uint8)
                                        image_normalized_b = outline_image_cell.astype(np.uint8)

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == True:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                        image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                        image_g[image_g > 255] = 255
                                        image_normalized_g = image_g.astype(np.uint8)

                                        image_b = np.add(outline_image_cell, outline_image_neighborhood)
                                        image_b[image_b > 255] = 255
                                        image_normalized_b = image_b.astype(np.uint8)

                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                    image = Image.fromarray(rgb_array, mode='RGB')
                                    croped_final_image = crop_image(image, postition_x, postition_y,crop_size)
                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                    if show_neighboring_cells == True:

                                        outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                        for pixel in cell_a_og:
                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        for index, cell_n in Cells_n.iterrows():

                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                for pixel in cell_a:
                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                        outline_image_cell_base[outline_image_cell_base > 255] = 255

                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                if len(proteins) == 2:

                                    # gets the position of the cell to crop the image
                                    postition_x = round(cell['Location_Center_X'])
                                    postition_y = round(cell['Location_Center_Y'])

                                    # gets the outline
                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell[outline_image_cell > 255] = 255

                                    # adds the neighborhood radius
                                    if show_neighborhood == True:

                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                        # adds the neighborhood coordinates to the blanck image
                                        for pixel in outline_image_neighborhood_cell:

                                            pixel = list(pixel)

                                            if pixel[0] < 0:
                                                pixel[0] = 0
                                            if pixel[1] < 0:
                                                pixel[1] = 0
                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[0] = (sub_list[1] - 2)
                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[1] = (sub_list[0] - 2)

                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                    path_name_1 = 'PathName_' + proteins[0]
                                    file_name_1 = 'FileName_' + proteins[0]
                                    path_name_2 = 'PathName_' + proteins[1]
                                    file_name_2 = 'FileName_' + proteins[1]

                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                    path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])

                                    image_1 = get_pixel_values(path_1)
                                    image_2 = get_pixel_values(path_2)

                                    image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                    image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value

                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast

                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines

                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == False:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                        image_normalized_g = outline_image_cell.astype(np.uint8)
                                        image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == True:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                        image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                        image_g[image_g > 255] = 255
                                        image_normalized_g = image_g.astype(np.uint8)

                                        image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                        image_b[image_b > 255] = 255
                                        image_normalized_b = image_b.astype(np.uint8)

                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                    image = Image.fromarray(rgb_array, mode='RGB')
                                    croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                    if show_neighboring_cells == True:

                                        outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                        for pixel in cell_a_og:
                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        for index, cell_n in Cells_n.iterrows():

                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                for pixel in cell_a:
                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                        outline_image_cell_base[outline_image_cell_base > 255] = 255

                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                if len(proteins) == 3:

                                    # gets the position of the cell to crop the image
                                    postition_x = round(cell['Location_Center_X'])
                                    postition_y = round(cell['Location_Center_Y'])

                                    # gets the outline
                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell[outline_image_cell > 255] = 255

                                    # adds the neighborhood radius
                                    if show_neighborhood == True:

                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                        # adds the neighborhood coordinates to the blanck image
                                        for pixel in outline_image_neighborhood_cell:

                                            pixel = list(pixel)

                                            if pixel[0] < 0:
                                                pixel[0] = 0
                                            if pixel[1] < 0:
                                                pixel[1] = 0
                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[0] = (sub_list[1] - 2)
                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                sub_list = list(outline_image_cell.shape)
                                                pixel[1] = (sub_list[0] - 2)

                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                    path_name_1 = 'PathName_' + proteins[0]
                                    file_name_1 = 'FileName_' + proteins[0]
                                    path_name_2 = 'PathName_' + proteins[1]
                                    file_name_2 = 'FileName_' + proteins[1]
                                    path_name_3 = 'PathName_' + proteins[2]
                                    file_name_3 = 'FileName_' + proteins[2]

                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                    path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])
                                    path_3 = str(cell[path_name_3]) + '/' + str(cell[file_name_3])

                                    image_1 = get_pixel_values(path_1)
                                    image_2 = get_pixel_values(path_2)
                                    image_3 = get_pixel_values(path_3)

                                    image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                    image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value
                                    image_normalized_3 = np.divide(image_3, max_pixel[2])  # normalizes the image to the choosen max pixel value

                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast
                                    image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])  # creates the contrast

                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines
                                    image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell)  # adds the outlines

                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255
                                    image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255  # sets everything above 255 to 255

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == False:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                        image_normalized_g = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                                        image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                    # generates the rgb image without neighborhood
                                    if show_neighborhood == True:
                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                        image_g = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                        image_g[image_g > 255] = 255
                                        image_normalized_g = image_g.astype(np.uint8)

                                        image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                        image_b[image_b > 255] = 255
                                        image_normalized_b = image_b.astype(np.uint8)

                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                    image = Image.fromarray(rgb_array, mode='RGB')
                                    croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                    if show_neighboring_cells == True:

                                        outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                        for pixel in cell_a_og:
                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                        for index, cell_n in Cells_n.iterrows():

                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                for pixel in cell_a:
                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                        outline_image_cell_base[outline_image_cell_base > 255] = 255

                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                        if select_neighbors[0] == True:

                            if assigned_cell_type_and_state[0] == cell_types_and_state[0] and set(cell_types_and_state[1:]).issubset(assigned_cell_type_and_state[1:]):

                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]

                                for index, cell_n in Cells_n.iterrows():

                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        cell_n_types = ast.literal_eval(cell_n['cell_types_and_states'])

                                        for cell_n_type in cell_n_types:

                                            if select_neighbors[1] == cell_n_type[0]:

                                                image_counter = image_counter + 1  # ticks up the counter for the image name

                                                if len(proteins) == 1:

                                                    # gets the position of the cell to crop the image
                                                    postition_x = round(cell['Location_Center_X'])
                                                    postition_y = round(cell['Location_Center_Y'])

                                                    # gets the outline
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                    outline_image_cell[outline_image_cell > 255] = 255

                                                    # adds the neighborhood radius
                                                    if show_neighborhood == True:

                                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                        # adds the neighborhood coordinates to the blanck image
                                                        for pixel in outline_image_neighborhood_cell:

                                                            pixel = list(pixel)

                                                            if pixel[0] < 0:
                                                                pixel[0] = 0
                                                            if pixel[1] < 0:
                                                                pixel[1] = 0
                                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[0] = (sub_list[1] - 2)
                                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[1] = (sub_list[0] - 2)

                                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    path_name_1 = 'PathName_' + proteins[0]
                                                    file_name_1 = 'FileName_' + proteins[0]

                                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])

                                                    image_1 = get_pixel_values(path_1)

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == False:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = outline_image_cell.astype(np.uint8)
                                                        image_normalized_b = outline_image_cell.astype(np.uint8)

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == True:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(outline_image_cell, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                    image = Image.fromarray(rgb_array, mode='RGB')
                                                    croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                                    if show_neighboring_cells == True:

                                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                        for pixel in cell_a_og:
                                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                        outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        for index, cell_n in Cells_n.iterrows():

                                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                                cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                                for pixel in cell_a:
                                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                                for c_t_a_s in cell_p:

                                                                    if select_neighbors[1] == c_t_a_s[0]:

                                                                        cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                                        for pixel in cell_p_a:
                                                                            outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                        outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                                if len(proteins) == 2:

                                                    # gets the position of the cell to crop the image
                                                    postition_x = round(cell['Location_Center_X'])
                                                    postition_y = round(cell['Location_Center_Y'])

                                                    # gets the outline
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                    outline_image_cell[outline_image_cell > 255] = 255

                                                    # adds the neighborhood radius
                                                    if show_neighborhood == True:

                                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                        # adds the neighborhood coordinates to the blanck image
                                                        for pixel in outline_image_neighborhood_cell:

                                                            pixel = list(pixel)

                                                            if pixel[0] < 0:
                                                                pixel[0] = 0
                                                            if pixel[1] < 0:
                                                                pixel[1] = 0
                                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[0] = (sub_list[1] - 2)
                                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[1] = (sub_list[0] - 2)

                                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    path_name_1 = 'PathName_' + proteins[0]
                                                    file_name_1 = 'FileName_' + proteins[0]
                                                    path_name_2 = 'PathName_' + proteins[1]
                                                    file_name_2 = 'FileName_' + proteins[1]

                                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                                    path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])

                                                    image_1 = get_pixel_values(path_1)
                                                    image_2 = get_pixel_values(path_2)

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == False:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = outline_image_cell.astype(np.uint8)
                                                        image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == True:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                    image = Image.fromarray(rgb_array, mode='RGB')
                                                    croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                                    if show_neighboring_cells == True:

                                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                        for pixel in cell_a_og:
                                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                        outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        for index, cell_n in Cells_n.iterrows():

                                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                                cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                                for pixel in cell_a:
                                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                                for c_t_a_s in cell_p:

                                                                    if select_neighbors[1] == c_t_a_s[0]:

                                                                        cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                                        for pixel in cell_p_a:
                                                                            outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                        outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                                if len(proteins) == 3:

                                                    # gets the position of the cell to crop the image
                                                    postition_x = round(cell['Location_Center_X'])
                                                    postition_y = round(cell['Location_Center_Y'])

                                                    # gets the outline
                                                    outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                    outline_image_cell[outline_image_cell > 255] = 255

                                                    # adds the neighborhood radius
                                                    if show_neighborhood == True:

                                                        outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                        outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                        # adds the neighborhood coordinates to the blanck image
                                                        for pixel in outline_image_neighborhood_cell:

                                                            pixel = list(pixel)

                                                            if pixel[0] < 0:
                                                                pixel[0] = 0
                                                            if pixel[1] < 0:
                                                                pixel[1] = 0
                                                            if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[0] = (sub_list[1] - 2)
                                                            if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                                sub_list = list(outline_image_cell.shape)
                                                                pixel[1] = (sub_list[0] - 2)

                                                            outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                                    path_name_1 = 'PathName_' + proteins[0]
                                                    file_name_1 = 'FileName_' + proteins[0]
                                                    path_name_2 = 'PathName_' + proteins[1]
                                                    file_name_2 = 'FileName_' + proteins[1]
                                                    path_name_3 = 'PathName_' + proteins[2]
                                                    file_name_3 = 'FileName_' + proteins[2]

                                                    path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                                    path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])
                                                    path_3 = str(cell[path_name_3]) + '/' + str(cell[file_name_3])

                                                    image_1 = get_pixel_values(path_1)
                                                    image_2 = get_pixel_values(path_2)
                                                    image_3 = get_pixel_values(path_3)

                                                    image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                                    image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value
                                                    image_normalized_3 = np.divide(image_3, max_pixel[2])  # normalizes the image to the choosen max pixel value

                                                    image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                                    image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast
                                                    image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])  # creates the contrast

                                                    image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                                    image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines
                                                    image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell)  # adds the outlines

                                                    image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                                    image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255
                                                    image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255  # sets everything above 255 to 255

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == False:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                        image_normalized_g = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                                                        image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                                    # generates the rgb image without neighborhood
                                                    if show_neighborhood == True:
                                                        image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                        image_g = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                        image_g[image_g > 255] = 255
                                                        image_normalized_g = image_g.astype(np.uint8)

                                                        image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                        image_b[image_b > 255] = 255
                                                        image_normalized_b = image_b.astype(np.uint8)

                                                    rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                                    image = Image.fromarray(rgb_array, mode='RGB')
                                                    croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                                    croped_final_image.save(export_file_path_full)  # Save to a file

                                                    if show_neighboring_cells == True:

                                                        outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                        for pixel in cell_a_og:
                                                            outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                        cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                        Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                        outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                        outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                        for index, cell_n in Cells_n.iterrows():

                                                            if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                                cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                                cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                                for pixel in cell_a:
                                                                    outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                                for c_t_a_s in cell_p:

                                                                    if select_neighbors[1] == c_t_a_s[0]:

                                                                        cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                                        for pixel in cell_p_a:
                                                                            outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                        outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                        outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                        outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                        rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                        image_n = Image.fromarray(rgb_array, mode='RGB')

                                                        croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                        export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                        croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                                break
                                        else:
                                            continue
                                        break

                # includes only cells that have that type
                if add_mixed_cells == False:

                    if select_neighbors[0] == False:

                        # checks if one of the assigned cell types matches the wnated type
                        if len(cell_type_and_states) == 2 and cell_type_and_states[0][0] == cell_types_and_state[0] and set(cell_types_and_state[1:]).issubset(cell_type_and_states[0][1:]):
                            print(cell_type_and_states)
                            image_counter = image_counter + 1  # ticks up the counter for the image name

                            if len(proteins) == 1:

                                # gets the position of the cell to crop the image
                                postition_x = round(cell['Location_Center_X'])
                                postition_y = round(cell['Location_Center_Y'])

                                # gets the outline
                                outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                outline_image_cell[outline_image_cell > 255] = 255

                                # adds the neighborhood radius
                                if show_neighborhood == True:

                                    outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                    # adds the neighborhood coordinates to the blanck image
                                    for pixel in outline_image_neighborhood_cell:

                                        pixel = list(pixel)

                                        if pixel[0] < 0:
                                            pixel[0] = 0
                                        if pixel[1] < 0:
                                            pixel[1] = 0
                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[0] = (sub_list[1] - 2)
                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[1] = (sub_list[0] - 2)

                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                path_name_1 = 'PathName_' + proteins[0]
                                file_name_1 = 'FileName_' + proteins[0]

                                path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])

                                image_1 = get_pixel_values(path_1)

                                image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value

                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast

                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines

                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255

                                # generates the rgb image without neighborhood
                                if show_neighborhood == False:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                    image_normalized_g = outline_image_cell.astype(np.uint8)
                                    image_normalized_b = outline_image_cell.astype(np.uint8)

                                # generates the rgb image without neighborhood
                                if show_neighborhood == True:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                    image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                    image_g[image_g > 255] = 255
                                    image_normalized_g = image_g.astype(np.uint8)

                                    image_b = np.add(outline_image_cell, outline_image_neighborhood)
                                    image_b[image_b > 255] = 255
                                    image_normalized_b = image_b.astype(np.uint8)

                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                image = Image.fromarray(rgb_array, mode='RGB')
                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                croped_final_image.save(export_file_path_full)  # Save to a file

                                if show_neighboring_cells == True:

                                    outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                    for pixel in cell_a_og:
                                        outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                    cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                    Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    for index, cell_n in Cells_n.iterrows():

                                        if int(cell_n['Cell_number']) in cell_neigboorhood:

                                            cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                            for pixel in cell_a:
                                                outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                    outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                    outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                    outline_image_cell_base[outline_image_cell_base > 255] = 255

                                    rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                    croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                            if len(proteins) == 2:

                                # gets the position of the cell to crop the image
                                postition_x = round(cell['Location_Center_X'])
                                postition_y = round(cell['Location_Center_Y'])

                                # gets the outline
                                outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                outline_image_cell[outline_image_cell > 255] = 255

                                # adds the neighborhood radius
                                if show_neighborhood == True:

                                    outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                    # adds the neighborhood coordinates to the blanck image
                                    for pixel in outline_image_neighborhood_cell:

                                        pixel = list(pixel)

                                        if pixel[0] < 0:
                                            pixel[0] = 0
                                        if pixel[1] < 0:
                                            pixel[1] = 0
                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[0] = (sub_list[1] - 2)
                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[1] = (sub_list[0] - 2)

                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                path_name_1 = 'PathName_' + proteins[0]
                                file_name_1 = 'FileName_' + proteins[0]
                                path_name_2 = 'PathName_' + proteins[1]
                                file_name_2 = 'FileName_' + proteins[1]

                                path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])

                                image_1 = get_pixel_values(path_1)
                                image_2 = get_pixel_values(path_2)

                                image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value

                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast

                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines

                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255

                                # generates the rgb image without neighborhood
                                if show_neighborhood == False:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                    image_normalized_g = outline_image_cell.astype(np.uint8)
                                    image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                # generates the rgb image without neighborhood
                                if show_neighborhood == True:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                    image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                    image_g[image_g > 255] = 255
                                    image_normalized_g = image_g.astype(np.uint8)

                                    image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                    image_b[image_b > 255] = 255
                                    image_normalized_b = image_b.astype(np.uint8)

                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                image = Image.fromarray(rgb_array, mode='RGB')
                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                croped_final_image.save(export_file_path_full)  # Save to a file

                                if show_neighboring_cells == True:

                                    outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                    for pixel in cell_a_og:
                                        outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                    cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                    Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    for index, cell_n in Cells_n.iterrows():

                                        if int(cell_n['Cell_number']) in cell_neigboorhood:

                                            cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                            for pixel in cell_a:
                                                outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                    outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                    outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                    outline_image_cell_base[outline_image_cell_base > 255] = 255

                                    rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                    croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                            if len(proteins) == 3:

                                # gets the position of the cell to crop the image
                                postition_x = round(cell['Location_Center_X'])
                                postition_y = round(cell['Location_Center_Y'])

                                # gets the outline
                                outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                outline_image_cell[outline_image_cell > 255] = 255

                                # adds the neighborhood radius
                                if show_neighborhood == True:

                                    outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                    outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                    # adds the neighborhood coordinates to the blanck image
                                    for pixel in outline_image_neighborhood_cell:

                                        pixel = list(pixel)

                                        if pixel[0] < 0:
                                            pixel[0] = 0
                                        if pixel[1] < 0:
                                            pixel[1] = 0
                                        if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[0] = (sub_list[1] - 2)
                                        if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                            sub_list = list(outline_image_cell.shape)
                                            pixel[1] = (sub_list[0] - 2)

                                        outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                path_name_1 = 'PathName_' + proteins[0]
                                file_name_1 = 'FileName_' + proteins[0]
                                path_name_2 = 'PathName_' + proteins[1]
                                file_name_2 = 'FileName_' + proteins[1]
                                path_name_3 = 'PathName_' + proteins[2]
                                file_name_3 = 'FileName_' + proteins[2]

                                path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])
                                path_3 = str(cell[path_name_3]) + '/' + str(cell[file_name_3])

                                image_1 = get_pixel_values(path_1)
                                image_2 = get_pixel_values(path_2)
                                image_3 = get_pixel_values(path_3)

                                image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value
                                image_normalized_3 = np.divide(image_3, max_pixel[2])  # normalizes the image to the choosen max pixel value

                                image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast
                                image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])  # creates the contrast

                                image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines
                                image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell)  # adds the outlines

                                image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255
                                image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255  # sets everything above 255 to 255

                                # generates the rgb image without neighborhood
                                if show_neighborhood == False:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                    image_normalized_g = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                                    image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                # generates the rgb image without neighborhood
                                if show_neighborhood == True:
                                    image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                    image_g = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                    image_g[image_g > 255] = 255
                                    image_normalized_g = image_g.astype(np.uint8)

                                    image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                    image_b[image_b > 255] = 255
                                    image_normalized_b = image_b.astype(np.uint8)

                                rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                image = Image.fromarray(rgb_array, mode='RGB')
                                croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                croped_final_image.save(export_file_path_full)  # Save to a file

                                if show_neighboring_cells == True:

                                    outline_image_cell_base = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                    outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                    for pixel in cell_a_og:
                                        outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                    cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                    Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                    outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                    for index, cell_n in Cells_n.iterrows():

                                        if int(cell_n['Cell_number']) in cell_neigboorhood:

                                            cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                            for pixel in cell_a:
                                                outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                    outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                    outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                    outline_image_cell_base[outline_image_cell_base > 255] = 255

                                    rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_base.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                    image_n = Image.fromarray(rgb_array, mode='RGB')

                                    croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                    export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                    croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                    if select_neighbors[0] == True:

                        if len(cell_type_and_states) == 2 and cell_type_and_states[0][0] == cell_types_and_state[0] and set(cell_types_and_state[1:]).issubset(cell_type_and_states[0][1:]):

                            cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                            Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]

                            for index, cell_n in Cells_n.iterrows():

                                if int(cell_n['Cell_number']) in cell_neigboorhood:

                                    cell_n_types = ast.literal_eval(cell_n['cell_types_and_states'])

                                    if len(cell_n_types) == 2 and select_neighbors[1] == cell_n_types[0][0]:

                                        image_counter = image_counter + 1  # ticks up the counter for the image name

                                        if len(proteins) == 1:

                                            # gets the position of the cell to crop the image
                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            # gets the outline
                                            outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                            outline_image_cell[outline_image_cell > 255] = 255

                                            # adds the neighborhood radius
                                            if show_neighborhood == True:

                                                outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                # adds the neighborhood coordinates to the blanck image
                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                            path_name_1 = 'PathName_' + proteins[0]
                                            file_name_1 = 'FileName_' + proteins[0]

                                            path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])

                                            image_1 = get_pixel_values(path_1)

                                            image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value

                                            image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast

                                            image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines

                                            image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == False:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                image_normalized_g = outline_image_cell.astype(np.uint8)
                                                image_normalized_b = outline_image_cell.astype(np.uint8)

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == True:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                                image_g[image_g > 255] = 255
                                                image_normalized_g = image_g.astype(np.uint8)

                                                image_b = np.add(outline_image_cell, outline_image_neighborhood)
                                                image_b[image_b > 255] = 255
                                                image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                            croped_final_image.save(export_file_path_full)  # Save to a file

                                            if show_neighboring_cells == True:

                                                outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                for pixel in cell_a_og:
                                                    outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                for index, cell_n in Cells_n.iterrows():

                                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                        cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                        cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                        for pixel in cell_a:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                        if len(cell_p) == 2 and select_neighbors[1] == cell_p[0][0]:

                                                            cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                            for pixel in cell_p_a:
                                                                outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                image_n = Image.fromarray(rgb_array, mode='RGB')

                                                croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                        if len(proteins) == 2:

                                            # gets the position of the cell to crop the image
                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            # gets the outline
                                            outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                            outline_image_cell[outline_image_cell > 255] = 255

                                            # adds the neighborhood radius
                                            if show_neighborhood == True:

                                                outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                # adds the neighborhood coordinates to the blanck image
                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                            path_name_1 = 'PathName_' + proteins[0]
                                            file_name_1 = 'FileName_' + proteins[0]
                                            path_name_2 = 'PathName_' + proteins[1]
                                            file_name_2 = 'FileName_' + proteins[1]

                                            path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                            path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])

                                            image_1 = get_pixel_values(path_1)
                                            image_2 = get_pixel_values(path_2)

                                            image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                            image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value

                                            image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                            image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast

                                            image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                            image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines

                                            image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                            image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == False:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                image_normalized_g = outline_image_cell.astype(np.uint8)
                                                image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == True:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                image_g = np.add(outline_image_cell, outline_image_neighborhood)
                                                image_g[image_g > 255] = 255
                                                image_normalized_g = image_g.astype(np.uint8)

                                                image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                image_b[image_b > 255] = 255
                                                image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                            croped_final_image.save(export_file_path_full)  # Save to a file

                                            if show_neighboring_cells == True:

                                                outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                for pixel in cell_a_og:
                                                    outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                for index, cell_n in Cells_n.iterrows():

                                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                        cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                        cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                        for pixel in cell_a:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                        if len(cell_p) == 2 and select_neighbors[1] == cell_p[0][0]:

                                                            cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                            for pixel in cell_p_a:
                                                                outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                image_n = Image.fromarray(rgb_array, mode='RGB')

                                                croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                        if len(proteins) == 3:

                                            # gets the position of the cell to crop the image
                                            postition_x = round(cell['Location_Center_X'])
                                            postition_y = round(cell['Location_Center_Y'])

                                            # gets the outline
                                            outline_image_cell = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                            outline_image_cell[outline_image_cell > 255] = 255

                                            # adds the neighborhood radius
                                            if show_neighborhood == True:

                                                outline_image_neighborhood = np.multiply(outline_image_cell, 0)  # creates a blank image
                                                outline_image_neighborhood_cell = ast.literal_eval(cell['Neighborhood_coordinates'])  # makes a list out of the cell outlines

                                                # adds the neighborhood coordinates to the blanck image
                                                for pixel in outline_image_neighborhood_cell:

                                                    pixel = list(pixel)

                                                    if pixel[0] < 0:
                                                        pixel[0] = 0
                                                    if pixel[1] < 0:
                                                        pixel[1] = 0
                                                    if pixel[0] >= (outline_image_neighborhood.shape[1] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[0] = (sub_list[1] - 2)
                                                    if pixel[1] >= (outline_image_neighborhood.shape[0] - 1):
                                                        sub_list = list(outline_image_cell.shape)
                                                        pixel[1] = (sub_list[0] - 2)

                                                    outline_image_neighborhood[pixel[1], pixel[0]] = 255

                                            path_name_1 = 'PathName_' + proteins[0]
                                            file_name_1 = 'FileName_' + proteins[0]
                                            path_name_2 = 'PathName_' + proteins[1]
                                            file_name_2 = 'FileName_' + proteins[1]
                                            path_name_3 = 'PathName_' + proteins[2]
                                            file_name_3 = 'FileName_' + proteins[2]

                                            path_1 = str(cell[path_name_1]) + '/' + str(cell[file_name_1])
                                            path_2 = str(cell[path_name_2]) + '/' + str(cell[file_name_2])
                                            path_3 = str(cell[path_name_3]) + '/' + str(cell[file_name_3])

                                            image_1 = get_pixel_values(path_1)
                                            image_2 = get_pixel_values(path_2)
                                            image_3 = get_pixel_values(path_3)

                                            image_normalized_1 = np.divide(image_1, max_pixel[0])  # normalizes the image to the choosen max pixel value
                                            image_normalized_2 = np.divide(image_2, max_pixel[1])  # normalizes the image to the choosen max pixel value
                                            image_normalized_3 = np.divide(image_3, max_pixel[2])  # normalizes the image to the choosen max pixel value

                                            image_normalized_contrasted_1 = np.multiply(image_normalized_1, contrast_multiplier[0])  # creates the contrast
                                            image_normalized_contrasted_2 = np.multiply(image_normalized_2, contrast_multiplier[1])  # creates the contrast
                                            image_normalized_contrasted_3 = np.multiply(image_normalized_3, contrast_multiplier[2])  # creates the contrast

                                            image_normalized_contrasted_add_outline_1 = np.add(image_normalized_contrasted_1, outline_image_cell)  # adds the outlines
                                            image_normalized_contrasted_add_outline_2 = np.add(image_normalized_contrasted_2, outline_image_cell)  # adds the outlines
                                            image_normalized_contrasted_add_outline_3 = np.add(image_normalized_contrasted_3, outline_image_cell)  # adds the outlines

                                            image_normalized_contrasted_add_outline_1[image_normalized_contrasted_add_outline_1 > 255] = 255  # sets everything above 255 to 255
                                            image_normalized_contrasted_add_outline_2[image_normalized_contrasted_add_outline_2 > 255] = 255  # sets everything above 255 to 255
                                            image_normalized_contrasted_add_outline_3[image_normalized_contrasted_add_outline_3 > 255] = 255  # sets everything above 255 to 255

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == False:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)
                                                image_normalized_g = image_normalized_contrasted_add_outline_3.astype(np.uint8)
                                                image_normalized_b = image_normalized_contrasted_add_outline_2.astype(np.uint8)

                                            # generates the rgb image without neighborhood
                                            if show_neighborhood == True:
                                                image_normalized_r = image_normalized_contrasted_add_outline_1.astype(np.uint8)

                                                image_g = np.add(image_normalized_contrasted_add_outline_3, outline_image_neighborhood)
                                                image_g[image_g > 255] = 255
                                                image_normalized_g = image_g.astype(np.uint8)

                                                image_b = np.add(image_normalized_contrasted_add_outline_2, outline_image_neighborhood)
                                                image_b[image_b > 255] = 255
                                                image_normalized_b = image_b.astype(np.uint8)

                                            rgb_array = np.stack((image_normalized_r, image_normalized_g, image_normalized_b), axis=-1)
                                            image = Image.fromarray(rgb_array, mode='RGB')
                                            croped_final_image = crop_image(image, postition_x, postition_y, crop_size)
                                            export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'  # creates export file path
                                            croped_final_image.save(export_file_path_full)  # Save to a file

                                            if show_neighboring_cells == True:

                                                outline_image_cell_a_og = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                cell_a_og = ast.literal_eval(cell['Cell_pixel_area'])

                                                for pixel in cell_a_og:
                                                    outline_image_cell_a_og[pixel[1], pixel[0]] = 255

                                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])
                                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']]
                                                outline_image_cell_neigboorhood = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])
                                                outline_image_cell_cell_types = get_pixel_values(file_dir_list_cell[cell['ImageNumber'] - 1])

                                                for index, cell_n in Cells_n.iterrows():

                                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                                        cell_a = ast.literal_eval(cell_n['Cell_pixel_area'])
                                                        cell_p = ast.literal_eval(cell_n['cell_types_and_states'])

                                                        for pixel in cell_a:
                                                            outline_image_cell_neigboorhood[pixel[1], pixel[0]] = 255

                                                        if len(cell_p) == 2 and select_neighbors[1] == cell_p[0][0]:

                                                            cell_p_a = ast.literal_eval(cell_n['Cell_pixel_area'])

                                                            for pixel in cell_p_a:
                                                                outline_image_cell_cell_types[pixel[1], pixel[0]] = 255

                                                outline_image_cell_neigboorhood[outline_image_cell_neigboorhood > 255] = 255
                                                outline_image_cell_a_og[outline_image_cell_a_og > 255] = 255
                                                outline_image_cell_cell_types[outline_image_cell_cell_types > 255] = 255

                                                rgb_array = np.stack((outline_image_cell_a_og.astype(np.uint8), outline_image_cell_cell_types.astype(np.uint8), outline_image_cell_neigboorhood.astype(np.uint8)), axis=-1)

                                                image_n = Image.fromarray(rgb_array, mode='RGB')

                                                croped_final_image_neigboorhood = crop_image(image_n, postition_x, postition_y, crop_size)

                                                export_file_path_full = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder_phenotypes + '/' + 'image_full_show_neighboring_cells_' + str(proteins[0]) + '_' + str(cell['ImageNumber']) + '_' + str(image_counter) + '.png'

                                                croped_final_image_neigboorhood.save(export_file_path_full)  # Save to a file

                                    else:
                                        continue
                                    break

#tests the ISO
def test_ISO(Cells, outline_cell_image_dir, output_folder, dir_name_string, folder_name, crop_size = 30):

    file_dir_list_cell = []  # list for the cell image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".tiff"):  # returns all files that end with .tiff
                file_dir_list_cell.append(filedir)  # appends the found .tiff files into the list

    output_folder = dir_name_string[:- len(str(folder_name)) - 17] + 'images' + '/' + output_folder  # generates the output folder

    # checks if the folder exists
    if os.path.isdir(output_folder) == True:
        pass
    else:
        os.makedirs(output_folder)

    try:
        for file in os.listdir(output_folder):
            if file.endswith('.png'):
                os.remove(str(output_folder) + '/' + str(file))

    finally:

        #goes throw all the cells
        image_counter = -1
        for index, cell in Cells.iterrows():

            if isinstance(cell['ISO_visualisation_and_information'],str): #checks if a cell was corrected

                image_counter = image_counter + 1

                cell_type_ISO = ast.literal_eval(cell['cell_types_and_states_ISO'])# gets the final cell type

                #gets the position of the cell
                postition_x = round(cell['Location_Center_X'])
                postition_y = round(cell['Location_Center_Y'])

                ISO_visualisation = ast.literal_eval(cell['ISO_visualisation_and_information']) #gets the information necessary for the image

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
def show_ISO(Cells, new_folder_name, outline_cell_image_dir):

    UniqueNames_Full_cell = Cells.ImageNumber.unique()  # creates a list of all unique images
    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  # creates a dictionary of all unique images

    # resets indexes for conected images and they identefied objects
    for key in DataFrameDict_Full_cell.keys():
        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

    file_dir_list_cell = []  # list for the cell image outline directories
    for paths, dirs, files in sd.walk(outline_cell_image_dir):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".tiff"):  # returns all files that end with .tiff
                file_dir_list_cell.append(filedir)  # appends the found .tiff files into the list

    # goes throw each image
    for key, value in DataFrameDict_Full_cell.items():

        outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])  # gets the cell outline image associated with the image key
        outline_image_cell[outline_image_cell > 255] = 255

        # gives you the file name
        filename = os.path.basename(file_dir_list_cell[key - 1])
        filename_string = str(filename)

        # gives you the directory name
        dir_name = os.path.dirname(file_dir_list_cell[key - 1])
        dir_name_string = str(dir_name)

        # gives you the folder name
        folder_name = os.path.basename(dir_name)
        folder_name_string = str(folder_name)

        filedir_test_image = dir_name_string[:-len(folder_name_string) - 17] + 'images/' + new_folder_name  # creates a new directory for the neighborhood analysis

        # checks if the new directory already exists
        if os.path.isdir(filedir_test_image) == True:
            pass
        else:
            os.makedirs(filedir_test_image)  # creates a new directory if necessary

        export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'  # creates a export file path

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
def area_manipulation_test(Cells, new_folder_name, outline_cell_image_dir, method = 'select_areas'):

    if method == 'select_areas':

        UniqueNames_Full_cell = Cells.ImageNumber.unique()  # creates a list of all unique images

        print(UniqueNames_Full_cell)

        DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueNames_Full_cell}  # creates a dictionary of all unique images

        # resets indexes for conected images and they identefied objects
        for key in DataFrameDict_Full_cell.keys():
            DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

        file_dir_list_cell = []  # list for the cell image outline directories
        for paths, dirs, files in sd.walk(outline_cell_image_dir):  # goes throw all files and folders in given directory

            for file in os.listdir(paths):  # goes throw all files in a folder
                filedir = os.path.join(paths, file)  # returns full file directory

                if filedir.endswith(".tiff"):  # returns all files that end with .tiff
                    file_dir_list_cell.append(filedir)  # appends the found .tiff files into the list

        # goes throw each image
        for key, value in DataFrameDict_Full_cell.items():

            outline_image_cell = get_pixel_values(file_dir_list_cell[key - 1])  # gets the cell outline image associated with the image key
            outline_image_cell[outline_image_cell > 255] = 255

            # gives you the file name
            filename = os.path.basename(file_dir_list_cell[key - 1])
            filename_string = str(filename)

            # gives you the directory name
            dir_name = os.path.dirname(file_dir_list_cell[key - 1])
            dir_name_string = str(dir_name)

            # gives you the folder name
            folder_name = os.path.basename(dir_name)
            folder_name_string = str(folder_name)

            filedir_test_image = dir_name_string[:-len(folder_name_string) - 17] + 'images/' + new_folder_name  # creates a new directory for the neighborhood analysis

            # checks if the new directory already exists
            if os.path.isdir(filedir_test_image) == True:
                pass
            else:
                os.makedirs(filedir_test_image)  # creates a new directory if necessary

            export_file_path = filedir_test_image + '/' + filename_string[:-5] + '_' + str(new_folder_name).replace(' ', '_') + '.png'  # creates a export file path

            cell_pixel_area = DataFrameDict_Full_cell[key]['Cell_pixel_area'].to_list()  # makes a list out of the cell outlines

            # assigns the cell area a gray scale value
            for cell_a in cell_pixel_area:
                cell_a = ast.literal_eval(cell_a)

                for pixel in cell_a:
                    outline_image_cell[pixel[1], pixel[0]] = 100

            data = Image.fromarray(outline_image_cell.astype(np.uint8))
            data.save(export_file_path)