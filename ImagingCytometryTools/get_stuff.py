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

    df = pd.read_csv(txt_file, sep="\t") #reads the .txt file provided 

    MCD_file_first_column = df.columns.values.tolist() #gets the column names that contain the metal isotope
    cleaned_metals = ['nan',] #empty list for metal isotopes. nan is required because the spillover matrix contains one empty filed.

    for element in MCD_file_first_column: #goes thow each column name
        metal = re.findall(r'\(.*?\)', element) #looks for brackets and takes the element within the brackets

        if bool(metal) == True: #checks if a metal isotope is found

            for string in metal: #goes trow the metal isotope
                string = string.replace('(', '').replace(')', '') #replaces the brackets
                cleaned_metals.append(string) #appends the finished extracted metal isotope into the list

    return(cleaned_metals) 

#extracts the compensation matrix from a .txt file
def get_spillover_matrix(x, seperation = "\t"):
    comp_mat_df = pd.read_csv(x, sep=seperation) #reads the .txt file provided
    return(comp_mat_df)

#extracts the pixel values from a .txt file
def get_file_for_compensation(x, seperation = "\t"):
    data = pd.read_csv(x, sep = seperation) #reads the .txt file provided
    return(data)

#gets the protein marker names from a CellProfiler .csv file
def get_markers_from_segmentation(Cells):
    first_column = list(Cells.columns.values.tolist()) # gets you the colum names from a data frame as list
    clean_markers = [] #empty list for protein markers

    for element in first_column:
        marker_cond = re.findall('MeanIntensity_', element) #condition to find the protein markers

        if bool(marker_cond) == True: #checks if a protein marker is found
            marker = element[element.find('MeanIntensity_'):] #find all the symbols after mean intensity
            marker_str = str(marker) #turn the marker into a string
            clean_markers.append(marker_str) # appends the protein marker into the list

    return(clean_markers)

#gets the markers from the segmentation dataset
def get_markers_from_csv(Cells):
    first_column = list(Cells.columns.values.tolist()) # gets you the colum names from a data frame as list
    clean_markers = [] #empty list for protein markers

    for element in first_column:
        marker_cond = re.findall('Intensity_MeanIntensity_', element) #condition to find the protein markers

        if bool(marker_cond) == True: #checks if a protein marker is found
            marker = element[element.find('Intensity_MeanIntensity_'):] #find all the symbols after mean intensity
            marker_str = str(marker)
            clean_markers.append(marker_str) # appends the protein marker into the list

    return(clean_markers)

#opens an image and extracts the pixel information as a numpy array
def get_pixel_values(x):
    img = Image.open(x) #imports image
    pixel_value_list = list(img.getdata()) #gets pixel data of image
    width, height = img.size #shapes the list to the image size
    pixel_values = np.array(pixel_value_list).reshape((height, width)) #turns the image into an array

    return(pixel_values)

#gets different types of maximum intensity pixel
def get_max_pixel(directory, protein_string, normalisation_value = 'mean'):

    max_pixel_list = []

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".tiff"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                if protein_string in filename_string and 'compensated images' in dir_name_string: #returns all files that contain the provided file string

                    #gets the max pixel for each found image
                    max_pixel = get_pixel_values(filedir).max()
                    max_pixel_list.append(max_pixel)

    #returns the mean of all max pixels
    if normalisation_value == 'mean':
        max_pixel_norm = np.mean(max_pixel_list)

    # returns the mean of all max pixels
    if normalisation_value == 'std':
        max_pixel_norm = np.mean(max_pixel_list) + np.std(max_pixel_list)

    # returns the mean of all max pixels
    if normalisation_value == 'max':
        max_pixel_norm = np.max(max_pixel_list)

    # returns the mean of all max pixels
    if normalisation_value == 'min':
        max_pixel_norm = np.min(max_pixel_list)
    
    return(max_pixel_norm)