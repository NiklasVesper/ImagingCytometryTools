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

def set_pixel_positive_area(Cells, proteins, filestring, new_folder_name, Intensity = [1,4294967295], Distance = 2, Neighbor_count = 1, Area = 5, apply_data = False):

    filename = os.path.basename(Cells)
    filename_string = str(filename)

    dir_name = os.path.dirname(Cells)
    dir_name_string = str(dir_name)

    folder_name = os.path.basename(dir_name)
    folder_name_string = str(folder_name)

    dir_name_string_2 = dir_name_string[:-len(str(folder_name))-1]
    folder_name_2 = os.path.basename(dir_name_string_2)
    dir_name_string_3 = dir_name_string_2[:-len(str(folder_name_2))-1]

    Cells = pd.read_csv(Cells)

    print('I am looking at: ' + str(proteins) + ' in: ' + str(os.path.join(dir_name, filename)))

    protein_counter = -1
    for protein in proteins:

        protein_counter = protein_counter + 1

        UniqueImageNumber = Cells.ImageNumber.unique()

        path_name = 'PathName_' + protein
        file_name = 'FileName_' + protein
        UniquePathName_1 = Cells[path_name].unique()
        UniqueFileName_1 = Cells[file_name].unique()

        unique_datas = []

        for number in UniqueImageNumber:

            img = Image.open(str(UniquePathName_1[number - 1]) + '/' + str(UniqueFileName_1[number - 1]))
            pixel_value_list = list(img.getdata())
            width, height = img.size
            pixel_values_analysis = np.array(pixel_value_list).reshape((height, width))

            pixel_values_visualisation = np.array(pixel_value_list).reshape((height, width))
            pixel_values_visualisation[pixel_values_visualisation > 0] = 0

            for (pixel_x, pixel_y), pixel_value in np.ndenumerate(pixel_values_analysis):

                if pixel_value <= Intensity[protein_counter][1] and pixel_value >= Intensity[protein_counter][0]:
                    pixel_values_analysis[pixel_x, pixel_y] = 255

                else:
                    pixel_values_analysis[pixel_x, pixel_y] = 0

            numbers = [0]
            for x in range(-Distance[protein_counter], Distance[protein_counter] + 1):
                if x < 0 or x > 0:
                    numbers.append(x)
                    numbers.append(x)

            directions = np.array(list(set(permutations(numbers, 2))))

            pixel_positive_areas = []

            for (pixel_x, pixel_y), pixel_value in np.ndenumerate(pixel_values_analysis):

                if pixel_value == 0 or (pixel_x, pixel_y) in chain(*pixel_positive_areas):
                    continue

                elif pixel_value == 255:

                    connected_pixels = [(pixel_x, pixel_y)]
                    stack = [(pixel_x, pixel_y)]

                    visited = set()
                    visited.add((pixel_x, pixel_y))

                    while stack:
                        x, y = stack.pop()

                        neighbouring_pixels = []

                        for direction in directions:
                            neighbouring_pixels.append((x + direction[0], y + direction[1]))

                        neighbouring_pixel_sublist = []
                        for neighbouring_pixel in neighbouring_pixels:

                            if 0 <= neighbouring_pixel[0] < height and 0 <= neighbouring_pixel[1] < width and pixel_values_analysis[neighbouring_pixel[0], neighbouring_pixel[1]] == 255:
                                neighbouring_pixel_sublist.append(neighbouring_pixel)

                        if len(neighbouring_pixel_sublist) >= Neighbor_count[protein_counter]:

                            for sub_pixel in neighbouring_pixel_sublist:

                                if sub_pixel not in visited:
                                    visited.add(sub_pixel)
                                    stack.append(sub_pixel)
                                    connected_pixels.append(sub_pixel)

                    if len(connected_pixels) >= Area[protein_counter]:
                        pixel_positive_areas.append(connected_pixels)

            unique_data = list(map(list, {tuple(sorted(sublist)) for sublist in pixel_positive_areas}))

            for ppa in unique_data:

                color = round(random.uniform(25, 240))
                for element in ppa:
                    pixel_values_visualisation[element[0], element[1]] = color

            output_folder = (dir_name_string_3 + '/' + 'images' + '/' + new_folder_name + '/' + protein
                             + '_Intensity-' + str(Intensity[protein_counter])
                             + '_Distance-' + str(Distance[protein_counter])
                             + '_Neighbor_count-' + str(Neighbor_count[protein_counter])
                             + '_Area-' + str(Area[protein_counter]))

            if os.path.isdir(output_folder) == True:
                pass
            else:
                os.makedirs(output_folder)

            export_file_path_full = (output_folder + '/' + str(UniqueFileName_1[number - 1])[:-5] + '.png')
            plt.imsave(export_file_path_full, pixel_values_visualisation, cmap='nipy_spectral')

            row_count = np.count_nonzero(Cells['ImageNumber'] == number)
            for x in range(0, row_count):
                unique_datas.append(unique_data)

        new_column = pd.DataFrame({'pixel_positive_area_' + protein: unique_datas})
        Cells = pd.concat([Cells, new_column], axis=1)

    if apply_data == True:

        filedir_pixel_positive_area = dir_name_string[:-len(folder_name_string)] + new_folder_name

        if os.path.isdir(filedir_pixel_positive_area) == True:
            pass
        else:
            os.makedirs(filedir_pixel_positive_area)

        export_file_path = filedir_pixel_positive_area + '/' + filename_string[:-len(filestring) - 4] + '_' +str(new_folder_name).replace(' ', '_') + '.csv'  # creates a export file path

        if os.path.isfile(export_file_path) == True:
            print('The file: ' + export_file_path + ' already exists!')

        else:
            Cells.to_csv(export_file_path) 