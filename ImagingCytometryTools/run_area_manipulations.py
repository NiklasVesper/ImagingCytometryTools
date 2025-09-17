import os
import scandir as sd
import pandas as pd
import pickle

from ImagingCytometryTools.interactive_tools import interactive_area_selection

from ImagingCytometryTools.area_manipulation import select_areas_from_dictionary
from ImagingCytometryTools.area_manipulation import select_areas_by_pixel_positive_area
from ImagingCytometryTools.area_manipulation import select_areas_from_both

'''
Runs the selection of cells intersecting with specially chosen areas.
'''

def run_interactive_area_selection(directory, filestring, export_file_path, protein, additional_images=None):

    area_selections = {}

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string:

                    print(filedir)

                    Cells = pd.read_csv(filedir)

                    UniqueImageNumber = Cells.ImageNumber.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                    for number in UniqueImageNumber:
                        Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                        image_name = Cells_on_image['ImageName'].unique()
                        image_number = Cells_on_image['ImageNumber'].unique()

                        gating_path_name = Cells_on_image['PathName_' + protein].unique()
                        gating_file_name = Cells_on_image['FileName_' + protein].unique()
                        gating_file_directory = gating_path_name[0] + '/' +gating_file_name[0]

                        if additional_images != None:

                            file_paths = []
                            for additional_protein in additional_images:

                                gating_path_name_add = Cells_on_image['PathName_' + additional_protein].unique()
                                gating_file_name_add = Cells_on_image['FileName_' + additional_protein].unique()
                                gating_file_directory_add = gating_path_name_add[0] + '/' + gating_file_name_add[0]

                                file_paths.append(gating_file_directory_add)

                            pixel_positions = interactive_area_selection(gating_file_directory, protein, file_paths=file_paths, additional_images=additional_images)
                            area_selections.update({image_name[0] + '_' + str(image_number[0]): pixel_positions})

                        else:
                            pixel_positions = interactive_area_selection(gating_file_directory,protein)
                            area_selections.update({image_name[0] + '_' + str(image_number[0]): pixel_positions})
    
    with open(export_file_path, 'wb') as d:
        pickle.dump(area_selections, d)

def run_select_areas(directory, filestring, new_folder_name, method):

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

                    filedir_area = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    if os.path.isdir(filedir_area) == True:
                        pass
                    else:
                        os.makedirs(filedir_area)

                    export_file_path = filedir_area + '/' + filename_string[:-len(filestring) - 4] + str(new_folder_name).replace(' ', '_') + '.csv'

                    if os.path.isfile(export_file_path) == True:
                        print('The file: ' + export_file_path + ' already exists!')
                        continue

                    else:

                        print('Select area for: ' + filedir)

                        if method[0] == 'from_dictionary':

                            file = pd.read_csv(filedir)
                            select_areas = select_areas_from_dictionary(file, method[1])
                            select_areas.to_csv(export_file_path)

                        if method[0] == 'by_pixel_positive_area':

                            file = pd.read_csv(filedir)
                            select_areas = select_areas_by_pixel_positive_area(file, method[1])
                            select_areas.to_csv(export_file_path)

                        if method[0] == 'both':
                            file = pd.read_csv(filedir)
                            select_areas = select_areas_from_both(file, method[1], method[2])
                            select_areas.to_csv(export_file_path)