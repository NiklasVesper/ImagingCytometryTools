import scandir as sd
import os

from ImagingCytometryTools.set_pixel_positive_area import set_pixel_positive_area

'''
Sets the pixel positive area for all files.
'''

def run_set_pixel_positive_area(directory, filestring, new_folder_name, proteins, Intensity = [0,4294967295], Distance = 2, Neighbor_count = 1, Area = 5, apply_data = False):

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                # gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string: #returns all files that contain the provided file string

                    #gives you the file directory
                    filedir = os.path.join(paths, file)
                    filedir_string = str(filedir)

                    #sets the pixel positive areas
                    set_pixel_positive_area(filedir_string, proteins, new_folder_name, Intensity, Distance, Neighbor_count, Area, apply_data)