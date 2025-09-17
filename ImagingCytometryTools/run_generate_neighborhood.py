import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.generate_neighborhood import generate_neighborhood
from ImagingCytometryTools.generate_neighborhood import generate_neighborhood_add_outline

'''
This function runs the functions that analyze the surrounding neighborhood of cells identified with the segmentation.

Here you can either choose specific radii or use the average cell size as an orientation point for the neighborhood radius.
This can also be combined with or without adding the cellular outlines.

You can also combine this with the analysis of the subcellular localization.
'''

def run_generate_neighborhood(directory, filestring, new_folder_name, add_cellular_information_and_outline = True, use_own_neighborhood_radius = [True, 0]):

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string:

                    filedir_neighborhood = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    if os.path.isdir(filedir_neighborhood) == True:
                        pass
                    else:
                        os.makedirs(filedir_neighborhood)

                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(str(os.path.basename(file))) - 4] + r'outline cell'

                    export_file_path = filedir_neighborhood + '/' + filename_string[:-len(filestring)-4] + str(new_folder_name).replace(' ', '_') + '.csv'

                    if add_cellular_information_and_outline == False:

                        if use_own_neighborhood_radius[0] == False:

                            if os.path.isfile(export_file_path) == True:
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                print('Running neighborhood analysis for: ' + str(filedir))
                                file = pd.read_csv(filedir)
                                neigboorhood_all = generate_neighborhood(file,
                                                                         file['AreaShape_MaxFeretDiameter'].mean())
                                neigboorhood_all.to_csv(export_file_path)

                        if use_own_neighborhood_radius[0] == True:

                            if os.path.isfile(export_file_path) == True:
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                if use_own_neighborhood_radius[1] == 0:
                                    print('You need to provide a value larger than 0 for the neighborhood radius!')
                                else:
                                    print('Running neighborhood analysis for: ' + str(filedir))
                                    file = pd.read_csv(filedir)
                                    neigboorhood_all = generate_neighborhood(file,
                                                                             use_own_neighborhood_radius[1])
                                    neigboorhood_all.to_csv(export_file_path)

                    elif add_cellular_information_and_outline == True:

                        if use_own_neighborhood_radius[0] == False:

                            if os.path.isfile(export_file_path) == True:
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                print('Running neighborhood analysis for: ' + str(filedir))
                                file = pd.read_csv(filedir)
                                neigboorhood_all = generate_neighborhood_add_outline(file,
                                                                                     file['AreaShape_MaxFeretDiameter'].mean(),
                                                                                     dir_name_string,
                                                                                     filename_string,
                                                                                     outline_cell)
                                neigboorhood_all.to_csv(export_file_path)

                        if use_own_neighborhood_radius[0] == True:

                            if os.path.isfile(export_file_path) == True:
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                if use_own_neighborhood_radius[1] == 0:
                                    print('You need to provide a value larger than 0 for the neighborhood radius!')
                                else:
                                    print('Running neighborhood analysis for: ' + str(filedir))
                                    file = pd.read_csv(filedir)
                                    neigboorhood_all = generate_neighborhood_add_outline(file,
                                                                                         use_own_neighborhood_radius[1],
                                                                                         dir_name_string,
                                                                                         filename_string,
                                                                                         outline_cell)
                                    neigboorhood_all.to_csv(export_file_path)

    print('I looked for my friends. ----------------------------------------------------------------------------------')