import os
import scandir as sd
import pandas as pd

from ImagingCytometryTools.image_analysis_of_spatial_overlap import image_analysis_of_spatial_overlap

'''
Runs the image analysis of spatial overlap function.
'''

def run_image_analysis_of_spatial_overlap(directory, filestring, new_folder_name, allowed_distance = 15, allowed_iterations = 3):

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

                    filedir_ISO = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    if os.path.isdir(filedir_ISO) == True:
                        pass
                    else:
                        os.makedirs(filedir_ISO)

                    export_file_path = filedir_ISO + '/' + filename_string[:-len(filestring) - 4] + '_' + str(new_folder_name).replace(' ', '_') + '.csv'  # creates a export file path

                    if os.path.isfile(export_file_path) == True:
                        print('The file: ' + export_file_path + ' already exists!')
                        continue
                    
                    else:
                        print('Running ISO for: ' + filedir)

                        file = pd.read_csv(filedir)
                        cells_after_ISO = image_analysis_of_spatial_overlap(file,
                                                                            allowed_distance,
                                                                            allowed_iterations)
                        cells_after_ISO.to_csv(export_file_path)