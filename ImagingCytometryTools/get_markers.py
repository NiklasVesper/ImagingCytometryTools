import re
import pandas as pd


#gets protein markers from a cellprofiler file
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
