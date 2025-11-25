import scandir as sd
import os

from ImagingCytometryTools.visualize_data import generate_image_galaries
from ImagingCytometryTools.get_data_from_files import get_max_pixel

'''
Function that generates image galaries for all files.
'''

def run_generate_image_galleries(directory, filestring, df_column, new_folder_name, max_pixel_images, proteins, normalisation_value = 'std', contrast_multiplier = [1,1,1], select_cell_type_and_state = [False], add_mixed_cells = True, generate_crops = False, crop_size = 40, show_neighborhood_radius = True, select_neighboring_cell_type_and_state = [False], show_neighboring_cells = False):

    max_pixels = []
    for potein in max_pixel_images:
        max_pixels.append(get_max_pixel(directory,potein,normalisation_value))

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)

                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string:

                    print('Generating image galleries for: ' + str(filedir))

                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(filename_string) - 1 - len(folder_name_string) -1] + '/' + r'outline cell'

                    generate_image_galaries(filedir,
                                            df_column,
                                            new_folder_name,
                                            outline_cell,
                                            proteins,
                                            max_pixels,
                                            contrast_multiplier=contrast_multiplier,
                                            select_cell_type_and_state=select_cell_type_and_state,
                                            add_mixed_cells=add_mixed_cells, 
                                            generate_crops=generate_crops, 
                                            crop_size=crop_size,
                                            show_neighborhood_radius=show_neighborhood_radius, 
                                            select_neighboring_cell_type_and_state=select_neighboring_cell_type_and_state,
                                            show_neighboring_cells=show_neighboring_cells)