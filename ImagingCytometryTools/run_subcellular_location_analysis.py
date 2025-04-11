import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.cell_to_organell import cell_to_organell_basic
from ImagingCytometryTools.cell_to_organell import cell_to_organell_advanced

'''
This function runs the functions that analyze the subcellular localization.
By selection of the number of nuclei dividing cells can be found.

Two options are avalible:
The first extracts the full cellular area and check whether the center of the nucleus is within that area.
The second checks whether the full nuclear area is within the cell.
'''

#runs the different types of subcellular localization analysis
def run_subcellular_location_analysis(directory, filestrings, new_folder_name,mode = 'basic', nucleus_count = 1):

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestrings[0] in filename_string: #returns all files that contain the provided file string

                    dir_list = [] #empty list for all the file directories

                    dir_list.append(filedir)

                    #gives you the directory name
                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    #gives you the folder name
                    folder_name = os.path.basename(dir_name)
                    folder_name_string = str(folder_name)

                    filedir_subcellular = dir_name_string[:-len(folder_name_string)] + new_folder_name #creates a new directory

                    #checks if the new directory already exists
                    if os.path.isdir(filedir_subcellular) == True:
                        pass
                    else:
                        os.makedirs(filedir_subcellular)

                    export_file_path = filedir_subcellular + '/' + filename_string[:-len(filestrings[0]) - 4] + str(new_folder_name).replace(' ', '_') + '.csv' #creates a export file path

                    #checks if the new file directory already exists
                    if os.path.isfile(export_file_path) == True:
                        print('The file: ' + export_file_path + ' already exists!')
                        continue

                    else:
                        for paths, dirs, files in sd.walk(dir_name_string):  #goes throw all files and folders in given directory

                            for file in os.listdir(paths):  #goes throw all files in a folder

                                #gives you the file directory
                                filedir = os.path.join(paths, file)
                                filedir_string = str(filedir)

                                #looks for the other files containing the information of the organelles
                                if filestrings[1] in filedir_string:
                                    dir_list.append(filedir)
                                if filestrings[2] in filedir_string:
                                    dir_list.append(filedir)

                                if len(dir_list) == 3: #continues if all three files are found

                                    print('Matching localization for: ' + str(dir_list))

                                    #gets the outline folder directories
                                    outline_cell = filedir_string[:-len(str(os.path.basename(file)))-4] + r'outline cell'
                                    outline_nucleus = filedir_string[:-len(str(os.path.basename(file)))-4] + r'outline nucleus'

                                    #assigns the found files
                                    Full_cell = pd.read_csv(dir_list[0])
                                    Cytoplasm = pd.read_csv(dir_list[1])
                                    Nucleus = pd.read_csv(dir_list[2])

                                    if mode == 'basic':
                                        single_cells_and_organells = cell_to_organell_basic(Full_cell, Cytoplasm,Nucleus, dir_name_string, filename_string, outline_cell, outline_nucleus, nucleus_count)
                                        single_cells_and_organells.to_csv(export_file_path)
                                    elif mode == 'advanced':
                                        single_cells_and_organells = cell_to_organell_advanced(Full_cell, Cytoplasm, Nucleus, dir_name_string, filename_string, outline_cell, outline_nucleus)
                                        single_cells_and_organells.to_csv(export_file_path)

    print('That was harder than finding then finding two matching socks. ---------------------------------------------')