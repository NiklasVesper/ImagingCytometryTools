import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.neigboorhood import neigboorhood
from ImagingCytometryTools.neigboorhood import neigboorhood_add_outline

'''
This function runs the functions that analyze the surrounding neighborhood of cells identified with the segmentation.

Here you can either choose specific radii or use the average cell size as an orientation point for the neighborhood radius.
This can also be combined with or without adding the cellular outlines.

You can also combine this with the analysis of the subcellular localization.
'''

#runs the different types of neighborhood analysis
def run_neighborhood_analysis(directory, filestring, new_folder_name, add_cellular_information_and_outline = True, use_own_neighborhood_radius = [True, 0]):

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                #gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                #gives you the folder name
                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string: #returns all files that contain the provided file string

                    filedir_neighborhood = dir_name_string[:-len(folder_name_string)] + new_folder_name #creates a new directory for the neighborhood analysis

                    #checks if the new directory already exists
                    if os.path.isdir(filedir_neighborhood) == True:
                        pass
                    else:
                        os.makedirs(filedir_neighborhood) #creates a new directory if necessary

                    #gets you the directory for the outline images
                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(str(os.path.basename(file))) - 4] + r'outline cell' #directory for cellular outline

                    export_file_path = filedir_neighborhood + '/' + filename_string[:-len(filestring)-4] + str(new_folder_name).replace(' ', '_') + '.csv' #creates a export file path

                    #runs neighborhood analysis for cells where subcellular information is added
                    if add_cellular_information_and_outline == False:

                        if use_own_neighborhood_radius[0] == False:

                            if os.path.isfile(export_file_path) == True: #checks if the new file directory already exists
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                print('Running neighborhood analysis for: ' + str(filedir))
                                file = pd.read_csv(filedir) #imports the file
                                neigboorhood_all = neigboorhood(file,file['AreaShape_MaxFeretDiameter'].mean()) #runs neighborhood analysis with average cell diameter as radius
                                neigboorhood_all.to_csv(export_file_path) #exports file

                        if use_own_neighborhood_radius[0] == True:

                            if os.path.isfile(export_file_path) == True: #checks if the new file directory already exists
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                if use_own_neighborhood_radius[1] == 0:
                                    print('You need to provide a value larger than 0 for the neighborhood radius!')
                                else:
                                    print('Running neighborhood analysis for: ' + str(filedir))
                                    file = pd.read_csv(filedir) #imports the file
                                    neigboorhood_all = neigboorhood(file,use_own_neighborhood_radius[1]) #runs neighborhood analysis with a fixed value as radius
                                    neigboorhood_all.to_csv(export_file_path) #exports file

                    #runs neighborhood analysis for full cells
                    elif add_cellular_information_and_outline == True:

                        if use_own_neighborhood_radius[0] == False:

                            if os.path.isfile(export_file_path) == True: #checks if the new file directory already exists
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                print('Running neighborhood analysis for: ' + str(filedir))
                                file = pd.read_csv(filedir) #imports the file
                                neigboorhood_all = neigboorhood_add_outline(file, file['AreaShape_MaxFeretDiameter'].mean(),dir_name_string,filename_string, outline_cell) #runs neighborhood analysis with average cell diameter as radius
                                neigboorhood_all.to_csv(export_file_path) #exports file

                        if use_own_neighborhood_radius[0] == True:

                            if os.path.isfile(export_file_path) == True: #checks if the new file directory already exists
                                print('The file: ' + export_file_path + ' already exists!')
                                continue
                            else:
                                if use_own_neighborhood_radius[1] == 0:
                                    print('You need to provide a value larger than 0 for the neighborhood radius!')
                                else:
                                    print('Running neighborhood analysis for: ' + str(filedir))
                                    file = pd.read_csv(filedir) #imports the file
                                    neigboorhood_all = neigboorhood_add_outline(file, use_own_neighborhood_radius[1], dir_name_string, filename_string, outline_cell) #runs neighborhood analysis with a fixed value as radius
                                    neigboorhood_all.to_csv(export_file_path) #exports file

    print('I looked for my friends. ----------------------------------------------------------------------------------')