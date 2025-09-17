import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.visualize_data import test_image_analysis_of_spatial_overlap

'''
Runs the test ISO function.
'''

def run_test_image_analysis_of_spatial_overlap(directory, filestring, output_folder, crop_size = 30):

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

                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(filename_string) - 1 - len(folder_name_string) - 1] + '/' + r'outline cell'

                    print('Testing ISO for: ' + filedir)

                    file = pd.read_csv(filedir)
                    test_image_analysis_of_spatial_overlap(file,
                                                           outline_cell,
                                                           output_folder,
                                                           dir_name_string,
                                                           folder_name,
                                                           crop_size=crop_size)