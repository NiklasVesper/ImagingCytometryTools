import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from itertools import permutations, groupby, chain
import random
import os

'''
Function that sets the pixel positive area and returns them as connected spots. 

Parameters that can be tweaked are the intensity, the distance allowed between identified pixels, 
the amount of neighboring pixels that need to be found as well as the size of the identified area.
'''

def set_pixel_positive_area(Cells, proteins, new_folder_name, Intensity = [1,4294967295], Distance = 2, Neighbor_count = 1, Area = 5, apply_data = False):

    # gives you the file name
    filename = os.path.basename(Cells)
    filename_string = str(filename)

    # gives you the directory name
    dir_name = os.path.dirname(Cells)
    dir_name_string = str(dir_name)

    # gives you the folder name
    folder_name = os.path.basename(dir_name)
    folder_name_string = str(folder_name)

    #goes three steps in the directory
    dir_name_string_2 = dir_name_string[:-len(str(folder_name))-1]
    folder_name_2 = os.path.basename(dir_name_string_2)
    dir_name_string_3 = dir_name_string_2[:-len(str(folder_name_2))-1]

    Cells = pd.read_csv(Cells) #imports the .csv file

    print('I am looking at: ' + str(proteins) + ' in: ' + str(os.path.join(dir_name, filename)))

    protein_counter = -1 #counter for the different proteins
    for protein in proteins: #goes throw all the different proteins

        protein_counter = protein_counter +1 #ticks up the protein counter

        UniqueImageNumber = Cells.ImageNumber.unique() #finds all unique image numbers

        #finds all unique path and file names
        path_name = 'PathName_' + protein
        file_name = 'FileName_' + protein
        UniquePathName_1 = Cells[path_name].unique()
        UniqueFileName_1 = Cells[file_name].unique()

        unique_datas = [] #empty list for unique pixel positive areas

        for number in UniqueImageNumber: #goes through all unique images

            img = Image.open(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))  # imports image
            pixel_value_list = list(img.getdata())  # gets pixel data of image
            width, height = img.size  # shapes the list to the image size
            pixel_values_og = np.array(pixel_value_list).reshape((height, width))  # turns the image into an array
            pixel_values_analysis = np.array(pixel_value_list).reshape((height, width))  # image array for analysis

            # image array for visualizing the analysis
            pixel_values_visualisation = np.array(pixel_value_list).reshape((height, width))  # image array for visualising the analysis
            pixel_values_visualisation[pixel_values_visualisation > 0] = 0 #sets all values in the array to 0

            for (pixel_x, pixel_y), pixel_value in np.ndenumerate(pixel_values_analysis):  # goes throw each pixel value

                if pixel_value <= Intensity[protein_counter][1] and pixel_value >= Intensity[protein_counter][0]:  # checks if the pixel value is within the intensity boundries
                    pixel_values_analysis[pixel_x, pixel_y] = 255  # assigns the value 255 if the value is within the intensity boundries

                else:
                    pixel_values_analysis[pixel_x, pixel_y] = 0  # assigns the value 0 if the value is outside the intensity boundries

            #creates a number range depending on the input distance
            numbers = [0]
            for x in range(-Distance[protein_counter], Distance[protein_counter] + 1):
                if x < 0 or x > 0:
                    numbers.append(x)
                    numbers.append(x)

            directions = np.array(list(set(permutations(numbers, 2))))  # creates an array of the possible directions

            pixel_positive_areas = []  #emty list for the areas

            for (pixel_x, pixel_y), pixel_value in np.ndenumerate(pixel_values_analysis):  # goes throw each pixel

                if pixel_value == 0 or (pixel_x, pixel_y) in chain(*pixel_positive_areas): #checks of the pixel is 0
                    continue

                elif pixel_value == 255: #checks of the pixel is 255

                    connected_pixels = [(pixel_x, pixel_y)]  #list of connected pixels
                    stack = [(pixel_x, pixel_y)]  # stack for the seach

                    visited = set()  # set for searched pixels
                    visited.add((pixel_x, pixel_y))

                    while stack:  # while there is a stack run the following
                        x, y = stack.pop()  # removes the last element of the stack

                        neighbouring_pixels = []  # list for neighbouring_pixels

                        # creates all the potential neighboring pixels to be checked
                        for direction in directions:
                            neighbouring_pixels.append((x + direction[0], y + direction[1]))

                        neighbouring_pixel_sublist = [] #empty list for the neighboring pixels
                        for neighbouring_pixel in neighbouring_pixels:

                            #if the neighboring pixel is within the image constrains and also positive it is appended in to the list
                            if 0 <= neighbouring_pixel[0] < height and 0 <= neighbouring_pixel[1] < width and pixel_values_analysis[neighbouring_pixel[0], neighbouring_pixel[1]] == 255:
                                neighbouring_pixel_sublist.append(neighbouring_pixel)

                        #checks if enough pixels are found to be added into the pixel positive area
                        if len(neighbouring_pixel_sublist) >= Neighbor_count[protein_counter]:

                            #adds the pixels to the lists and sets if there have not been visited
                            for sub_pixel in neighbouring_pixel_sublist:

                                if sub_pixel not in visited:
                                    visited.add(sub_pixel)
                                    stack.append(sub_pixel)
                                    connected_pixels.append(sub_pixel)

                    #checks if the area of connected pixels is above the threshold
                    if len(connected_pixels) >= Area[protein_counter]:
                        pixel_positive_areas.append(connected_pixels)

            unique_data = list(map(list, {tuple(sorted(sublist)) for sublist in pixel_positive_areas})) #removes all multiples evern tow there shoulb be none

            #visualisation for the areas
            for ppa in unique_data:

                #choses a random value for the area color
                color = round(random.uniform(25, 240))
                for element in ppa:
                    pixel_values_visualisation[element[0], element[1]] = color

            #creates an output folder based input settings
            output_folder = (dir_name_string_3 + '/' + new_folder_name + '/' + protein
                             + '_Intensity-' + str(Intensity[protein_counter])
                             + '_Distance-' + str(Distance[protein_counter])
                             + '_Neighbor_count-' + str(Neighbor_count[protein_counter])
                             + '_Area-' + str(Area[protein_counter])
                             )

            #checks if the folder already exists
            if os.path.isdir(output_folder) == True:
                pass
            else:
                os.makedirs(output_folder)

            #creates an export path and saves the image
            export_file_path_full = (output_folder + '/' + str(UniqueFileName_1[number - 1])[:-5] + '.png')
            plt.imsave(export_file_path_full, pixel_values_visualisation, cmap='nipy_spectral')

            #appends the pixel positive ares equal to the number of rows
            row_count = np.count_nonzero(Cells['ImageNumber'] == number)
            for x in range(0, row_count):
                unique_datas.append(unique_data)

        new_column = pd.DataFrame({'pixel_positive_area_' + protein: unique_datas}) #creates a new dataframe column name based on the protein name
        Cells = pd.concat([Cells, new_column], axis=1) #merges the dataframes

    if apply_data == True: #exports the file if wanted

        filedir_pixel_positive_area = dir_name_string[:-len(folder_name_string)] + new_folder_name  # creates a new directory

        # checks if the new directory already exists
        if os.path.isdir(filedir_pixel_positive_area) == True:
            pass
        else:
            os.makedirs(filedir_pixel_positive_area)

        export_file_path = filedir_pixel_positive_area + '/' + filename_string[: - 4] + '_' + 'pixel_positive_area' + '.csv'  # creates a export file path

        #checks if the export file path already exists
        if os.path.isfile(export_file_path) == True:
            print('The file: ' + export_file_path + ' already exists!')

        else:
            Cells.to_csv(export_file_path) #exports