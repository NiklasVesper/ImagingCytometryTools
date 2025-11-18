import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.generate_subcellular_localization import generate_subcellular_localization_basic
from ImagingCytometryTools.generate_subcellular_localization import generate_subcellular_localization_advanced

'''
This function runs the functions that analyze the subcellular localization.
By selection of the number of nuclei dividing cells can be found.

Two options are avalible:
The first extracts the full cellular area and check whether the center of the nucleus is within that area.
The second checks whether the full nuclear area is within the cell.
'''

def run_generate_subcellular_localization(directory, filestrings, new_folder_name, mode = 'basic', nucleus_count = 1):

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestrings[0] in filename_string:

                    dir_list = []

                    dir_list.append(filedir)

                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    folder_name = os.path.basename(dir_name)
                    folder_name_string = str(folder_name)

                    filedir_subcellular = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    if os.path.isdir(filedir_subcellular) == True:
                        pass
                    else:
                        os.makedirs(filedir_subcellular)

                    export_file_path = filedir_subcellular + '/' + filename_string[:-len(filestrings[0]) - 4] + str(new_folder_name).replace(' ', '_') + '.csv' #creates a export file path

                    if os.path.isfile(export_file_path) == True:
                        print('The file: ' + export_file_path + ' already exists!')
                        continue

                    else:
                        for paths, dirs, files in sd.walk(dir_name_string):

                            for file in os.listdir(paths):

                                filedir = os.path.join(paths, file)
                                filedir_string = str(filedir)

                                if filestrings[1] in filedir_string:
                                    dir_list.append(filedir)
                                if filestrings[2] in filedir_string:
                                    dir_list.append(filedir)

                                if len(dir_list) == 3:

                                    print('Matching localization for: ' + str(dir_list))

                                    outline_cell = filedir_string[:-len(str(os.path.basename(file)))-4] + r'outline cell'
                                    outline_nucleus = filedir_string[:-len(str(os.path.basename(file)))-4] + r'outline nucleus'

                                    Full_cell = pd.read_csv(dir_list[0])
                                    Cytoplasm = pd.read_csv(dir_list[1])
                                    Nucleus = pd.read_csv(dir_list[2])

                                    if mode == 'basic':
                                        single_cells_and_organells = generate_subcellular_localization_basic(Full_cell,
                                                                                                             Cytoplasm,
                                                                                                             Nucleus,
                                                                                                             dir_name_string,
                                                                                                             filename_string,
                                                                                                             outline_cell,
                                                                                                             outline_nucleus,
                                                                                                             nucleus_count)
                                        single_cells_and_organells.to_csv(export_file_path)

                                    elif mode == 'advanced':
                                        single_cells_and_organells = generate_subcellular_localization_advanced(Full_cell,
                                                                                                                Cytoplasm,
                                                                                                                Nucleus,
                                                                                                                dir_name_string,
                                                                                                                filename_string,
                                                                                                                outline_cell,
                                                                                                                outline_nucleus)
                                        single_cells_and_organells.to_csv(export_file_path)

    print('That was harder than finding then finding two matching socks. ---------------------------------------------')