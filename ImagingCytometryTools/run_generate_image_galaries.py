import scandir as sd
import os

from ImagingCytometryTools.visualize_images import generate_image_galaries
from ImagingCytometryTools.get_stuff import get_max_pixel

'''
Function that generates image galaries for all files.
'''

def run_generate_image_galaries(directory, filestring, new_folder_name, cell_types_and_states, proteins, normalisation_value = 'std', contrast_multiplier = [1,1,1], crop_size = 40, inclue_all = True, show_neighborhood = True, select_neighbors = [False], show_neighboring_cells = False, generate_crops = True):

    #gets the max pixels for each protein
    max_pixels = []
    for potein in proteins:
        max_pixels.append(get_max_pixel(directory,potein,normalisation_value))

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                #gives you the directory name
                dir_name = os.path.dirname(filedir)

                #gives you the folder name
                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string: #returns all files that contain the provided file string

                    # gets you the directory for the outline images
                    filedir_string = str(filedir)
                    outline_cell = filedir_string[:-len(filename_string) - 1 - len(folder_name_string) -1] + '/' + r'outline cell'  # directory for cellular outline

                    generate_image_galaries(filedir,
                                            new_folder_name,
                                            outline_cell,
                                            proteins,
                                            cell_types_and_states,
                                            max_pixels,
                                            contrast_multiplier=contrast_multiplier,
                                            crop_size=crop_size,
                                            add_mixed_cells=inclue_all,
                                            show_neighborhood=show_neighborhood,
                                            select_neighbors=select_neighbors,
                                            show_neighboring_cells=show_neighboring_cells,
                                            generate_crops=generate_crops)