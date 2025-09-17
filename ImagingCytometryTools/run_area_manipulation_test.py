import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.visualize_data import area_manipulation_test

'''
Visualizes the selected cells so it can be evaluated if the correct cells were chosen.
'''

def run_area_manipulation_test(directory, filestring, new_folder_name):

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string:

                    filedir = os.path.join(paths, file)
                    filedir_string = str(filedir)

                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    folder_name = os.path.basename(dir_name)
                    folder_name_str = str(folder_name)

                    filedir_outline_image = dir_name_string[:-len(folder_name_str)] + r'outline cell'
                    filedir_pd = pd.read_csv(filedir)

                    print('Testing: ' + str(filedir))

                    area_manipulation_test(filedir_pd, new_folder_name, filedir_outline_image)