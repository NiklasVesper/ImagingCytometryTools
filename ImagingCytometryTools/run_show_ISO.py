import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.visualize_images import show_ISO

'''
Runs the function that shows the cells that were corrected by ISO.
'''

def run_show_ISO(directory, filestring, output_folder):

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

                    # gets you the directory for the outline images
                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(filename_string) - 1 - len(folder_name_string) - 1] + '/' + r'outline cell'  # directory for cellular outline

                    print('Showing ISO for: ' + filedir)

                    file = pd.read_csv(filedir) #imports the file
                    cells_after_ISO = show_ISO(file, output_folder, outline_cell)