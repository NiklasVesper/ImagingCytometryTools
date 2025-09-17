import scandir as sd
import os

from ImagingCytometryTools.set_pixel_positive_area import set_pixel_positive_area

'''
Sets the pixel positive area for all files.
'''

def run_set_pixel_positive_area(directory, filestring, new_folder_name, proteins, Intensity = [0,4294967295], Distance = 2, Neighbor_count = 1, Area = 5, apply_data = False):

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

                    filedir_pixel_positive_area = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    export_file_path = filedir_pixel_positive_area + '/' + filename_string[:-len(filestring) - 4] + '_' + str(new_folder_name).replace(' ', '_') + '.csv'  # creates a export file path

                    if os.path.isfile(export_file_path) == True:
                        print('The file: ' + export_file_path + ' already exists!')

                    else:
                        set_pixel_positive_area(filedir_string,
                                                proteins,
                                                filestring,
                                                new_folder_name,
                                                Intensity,
                                                Distance,
                                                Neighbor_count,
                                                Area,
                                                apply_data)