import re
import pandas as pd
import numpy as np
from PIL import Image
import scandir as sd
import os

'''
A collection of functions to extract information from files.
'''

#extracts the used metal isotopes from a .txt file belonging to a .mcd file
def get_metals_from_IMC_data(txt_file):

    df = pd.read_csv(txt_file, sep="\t")

    MCD_file_first_column = df.columns.values.tolist()
    cleaned_metals = ['nan',]

    for element in MCD_file_first_column:
        metal = re.findall(r'\(.*?\)', element)

        if bool(metal) == True:

            for string in metal:
                string = string.replace('(', '').replace(')', '')
                cleaned_metals.append(string)

    return(cleaned_metals) 

#extracts the compensation matrix from a .txt file
def get_spillover_matrix(x, seperation = "\t"):
    comp_mat_df = pd.read_csv(x, sep=seperation)
    return(comp_mat_df)

#extracts the pixel values from a .txt file
def get_file_for_compensation(x, seperation = "\t"):
    data = pd.read_csv(x, sep = seperation)
    return(data)

#gets the protein marker names from a CellProfiler .csv file
def get_markers_from_segmentation(Cells):
    first_column = list(Cells.columns.values.tolist())
    clean_markers = []

    for element in first_column:
        marker_cond = re.findall('MeanIntensity_', element)

        if bool(marker_cond) == True:
            marker = element[element.find('MeanIntensity_'):]
            marker_str = str(marker)
            clean_markers.append(marker_str)

    return(clean_markers)

#gets the markers from the segmentation dataset
def get_markers_from_csv(Cells):
    first_column = list(Cells.columns.values.tolist())
    clean_markers = []

    for element in first_column:
        marker_cond = re.findall('Intensity_MeanIntensity_', element)

        if bool(marker_cond) == True:
            marker = element[element.find('Intensity_MeanIntensity_'):]
            marker_str = str(marker)
            clean_markers.append(marker_str)

    return(clean_markers)

#opens an image and extracts the pixel information as a numpy array
def get_pixel_values(x):
    img = Image.open(x)
    pixel_value_list = list(img.getdata())
    width, height = img.size
    pixel_values = np.array(pixel_value_list).reshape((height, width))

    return(pixel_values)

#gets different types of maximum intensity pixel
def get_max_pixel(directory, protein_string, normalisation_value = 'mean'):

    max_pixel_list = []

    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".tiff"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                if protein_string in filename_string and 'compensated images' in dir_name_string:

                    max_pixel = get_pixel_values(filedir).max()
                    max_pixel_list.append(max_pixel)

    if normalisation_value == 'mean':
        max_pixel_norm = np.mean(max_pixel_list)

    if normalisation_value == 'std':
        max_pixel_norm = np.mean(max_pixel_list) + np.std(max_pixel_list)

    if normalisation_value == 'max':
        max_pixel_norm = np.max(max_pixel_list)

    if normalisation_value == 'min':
        max_pixel_norm = np.min(max_pixel_list)
    
    return(max_pixel_norm)