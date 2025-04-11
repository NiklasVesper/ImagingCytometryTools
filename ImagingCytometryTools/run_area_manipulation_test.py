import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.visualize_images import area_manipulation_test

'''
Visualizes the selected cells so it can be evaluated if the correct cells were chosen.
'''

def run_area_manipulation_test(directory, filestring, new_folder_name, method='select_areas'):

    for paths, dirs, files in sd.walk(directory):  # goes throw all files and folders in given directory

        for file in os.listdir(paths):  # goes throw all files in a folder
            filedir = os.path.join(paths, file)  # returns full file directory

            if filedir.endswith(".csv"):  # returns all files that end with .csv

                # gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string:  # returns all files that contain the provided file string

                    filedir = os.path.join(paths, file)
                    filedir_string = str(filedir)

                    # gives you the directory name
                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    # gives you the folder name
                    folder_name = os.path.basename(dir_name)
                    folder_name_str = str(folder_name)

                    filedir_outline_image = dir_name_string[:-len(folder_name_str)] + r'outline cell'  # gets you the directory for the outline images
                    filedir_pd = pd.read_csv(filedir)  # imports the data frame

                    print('Testing: ' + str(filedir))

                    area_manipulation_test(filedir_pd, new_folder_name, filedir_outline_image, method=method)  # does the test