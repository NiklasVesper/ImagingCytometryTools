import re
import pandas as pd
import scipy as sc
import numpy as np

# Gets the metals from a MCD dataset
def get_metals_from_data(x):

    df = pd.read_csv(x, sep="\t") #reads the MCD file provided

    MCD_file_first_column = df.columns.values.tolist() # Gets the keys from the pandas Dataframe as list
    cleaned_metals = ['nan',] # nan is reqired because the spillover matrix contains one empty filed

    for element in MCD_file_first_column:
        metal = re.findall(r'\(.*?\)', element) # looks for brackets and takes the element within the brackets

        if bool(metal) == True: # if there is an element that fits the description do the following

            for string in metal: # go throw each symbol in the metal
                string = string.replace('(', '').replace(')', '')  # remove bracket
                cleaned_metals.append(string)  # append to metal list

    return(cleaned_metals) # returns metal list

# Gets the markers from a MCD dataset
def get_markers_from_data(x): # Helper funktion to get the markers from a MCD dataset

    df = pd.read_csv(x, sep="\t")
    MCD_file_first_column = df.columns.values.tolist() # gets you the keys from a pd df as list
    clean_markers = ['nan',] # empty list for markers

    for element in MCD_file_first_column: # gets you each element in the key list
        metal_cond = re.findall(r'\(.*?\)', element) # Here only a condition to find the metals

        if bool(metal_cond) == True: # if there is an element that fits the description do the following
                marker = element[:element.find("(")] # find all the symbols infront of '('

                marker_str = str(marker) # turns marker into strings for the dictionary
                clean_markers.append(marker_str) # append to marker list

    return(clean_markers) # returns marker list

#
def get_comp_mat(x):

    comp_mat_df = pd.read_csv(x, sep="\t")

    return(comp_mat_df)

#
def remove_non_recorded_channels_and_invert(mat,chan):
    channel_indexes = [mat.columns.get_loc(c) - 1 for c in chan if c in mat]
    comp_mat_final = mat.iloc[: , 1:]

    for c in chan:
        comp_mat_final = comp_mat_final.drop(c, axis = 1)

    comp_mat_final = comp_mat_final.drop(index = channel_indexes)
    comp_mat_final_np = comp_mat_final.to_numpy()
    comp_mat_final_inv = np.linalg.inv(np.float64(comp_mat_final_np))


    '''
    comp_mat_temp_np = np.array([get_metals_from_data(pd.read_csv(filedir, sep="\t"))[1:]])  # returns metals as a np array
    for element in x:  # goes throw each element in metals list
        selected_row = pd.DataFrame(sm[x][sm[x].nan == element])
        comp_mat_temp_np = np.append(comp_mat_temp_np, selected_row.iloc[:, 1:].to_numpy(), axis=0)
    comp_mat_inv = np.linalg.inv(np.float64(comp_mat_temp_np[1:]))
    '''
    return(comp_mat_final_inv)

def get_data_for_comp(x):
    data = pd.read_csv(x, sep="\t")
    return(data)

#----------------------





def get_mat_from_data_for_comp(x,dic,lis):
    data = pd.read_csv(x, sep="\t")
    data = data.rename(columns = dic)[lis[1:]]
    #data_np = pd.DataFrame(data).to_numpy()
    return(data)

def get_mat_from_data_for_XYZ(x):
    data = pd.read_csv(x, sep="\t")
    return(data)
'''
def get_comp_mat_from_data(dir,x):
    comp_mat_temp_np = np.array(x)  # returns metals as a np array
    for element in x:  # goes throw each element in metals list
        selected_row = pd.DataFrame(dir[x][dir[x].nan == element])
        comp_mat_temp_np = np.append(comp_mat_temp_np, selected_row.iloc[:, 1:].to_numpy(), axis=0)
    comp_mat_inv = sc.linalg.inv(comp_mat_temp_np[1:]) # np.float64()
    comp_mat_inv_pd = pd.DataFrame(comp_mat_inv, columns = metals[1:])
    return(comp_mat_inv_pd)
'''
def get_comp_mat_from_data(x):
    comp_mat_temp_np = np.array([get_metals_from_data(pd.read_csv(filedir, sep="\t"))[1:]])  # returns metals as a np array
    for element in x:  # goes throw each element in metals list
        selected_row = pd.DataFrame(sm[x][sm[x].nan == element])
        comp_mat_temp_np = np.append(comp_mat_temp_np, selected_row.iloc[:, 1:].to_numpy(), axis=0)
    comp_mat_inv = np.linalg.inv(np.float64(comp_mat_temp_np[1:]))
    return(comp_mat_inv)

def get_mat_from_data(x):
    data = pd.read_csv(x, sep="\t")
    data = data.rename(columns = Metals)[metals[1:]]
    data_np = pd.DataFrame(data).to_numpy()
    return(data_np)

def get_comp_mat_from_data(x):
    comp_mat_temp_np = np.array([get_metals_from_data(pd.read_csv(filedir, sep="\t"))[1:]])  # returns metals as a np array
    for element in x:  # goes throw each element in metals list
        selected_row = pd.DataFrame(sm[x][sm[x].nan == element])
        comp_mat_temp_np = np.append(comp_mat_temp_np, selected_row.iloc[:, 1:].to_numpy(), axis=0)
    comp_mat_inv = np.linalg.inv(np.float64(comp_mat_temp_np[1:]))
    return(comp_mat_inv)

#
def get_markers_from_segmentation(Cells): # Helper funktion to get the markers from the segmentation dataset
    first_column = list(Cells.columns.values.tolist()) # gets you the keys from a pd df as list
    clean_markers = [] # empty list for markers
    for element in first_column: # gets you each element in the key list
        marker_cond = re.findall('MeanIntensity_', element) # Here only a condition to find the markers
        if bool(marker_cond) == True: # if there is an element that fits the description do the following
            marker = element[element.find('MeanIntensity_'):] # find all the symbols after mean intensity
            marker_str = str(marker) # turns marker into strings for the list
            clean_markers.append(marker_str) # append to marker list
    return(clean_markers) # returns marker list