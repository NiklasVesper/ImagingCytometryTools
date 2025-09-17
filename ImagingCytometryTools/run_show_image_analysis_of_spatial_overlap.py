import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.visualize_data import show_image_analysis_of_spatial_overlap

'''
Runs the function that shows the cells that were corrected by ISO.
'''

def run_show_image_analysis_of_spatial_overlap(directory, filestring, output_folder):

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

                    print('Showing ISO for: ' + filedir)

                    file = pd.read_csv(filedir)
                    cells_after_ISO = show_image_analysis_of_spatial_overlap(file,
                                                                             output_folder,
                                                                             outline_cell)