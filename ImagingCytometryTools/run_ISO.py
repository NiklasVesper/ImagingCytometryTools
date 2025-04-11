import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.ISO import ISO

'''
Runs the ISO function.
'''

def run_ISO(directory, filestring, new_folder_name, allowed_distance = 15, allowed_iterations = 3):

    for paths, dirs, files in sd.walk(directory):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".csv"):  # returns all files that end with .csv

                # gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                # gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                # gives you the folder name
                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string:  # returns all files that contain the provided file string

                    filedir_ISO = dir_name_string[:-len(folder_name_string)] + new_folder_name  # creates a new directory for the neighborhood analysis

                    # checks if the new directory already exists
                    if os.path.isdir(filedir_ISO) == True:
                        pass
                    else:
                        os.makedirs(filedir_ISO)  # creates a new directory if necessary

                    export_file_path = filedir_ISO + '/' + filename_string[:-4] + '_iso' + '.csv'  # creates a export file path

                    if os.path.isfile(export_file_path) == True:  # checks if the new file directory already exists
                        print('The file: ' + export_file_path + ' already exists!')
                        continue
                    
                    else:

                        print('Running ISO for: ' + filedir)

                        file = pd.read_csv(filedir)
                        cells_after_ISO = ISO(file,allowed_distance,allowed_iterations)
                        cells_after_ISO.to_csv(export_file_path)