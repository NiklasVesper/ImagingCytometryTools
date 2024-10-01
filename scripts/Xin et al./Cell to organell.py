import pandas as pd
import os
import scandir as sd

from ImagingCytometryToolsIMC.cell_to_organell import cell_to_organell_basic
from ImagingCytometryToolsIMC.neigboorhood import neigboorhood


folder_dir = r'D:\Basel'# folder directory


for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith(".csv"):# returns all files that end with txt

            filename = os.path.basename(file) # gives you the file name
            filename_string = str(filename)

            if 'RunCellpose_C' in filename_string: # checks for a condition in the string
                print(filedir)

                filedir_string = str(filedir)[:-len(filename)-5] # creates a file directory string

                filedir_images_string = filedir_string + r'\neigboorhood' # creates a new directory

                if os.path.isdir(filedir_images_string) == True: #checks if the new directory allready exists
                    pass
                else:
                    os.makedirs(filedir_images_string) # creates the new directory

                neigboorhood_all = neigboorhood(pd.read_csv(filedir)) # calculates the neighborhood
                neigboorhood_all.to_csv(filedir_images_string + '/' + filename_string[:-18] + '_neighborhood.csv') # exports the analysis


for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith(".csv"): # returns all files that end with txt

            filename = os.path.basename(file) # gives you the file name
            filename_string = str(filename)

            if 'RunCellpose_C' in filename_string: # checks for a condition in the string
                dir_list = [] # creates a list to hold the other directories

                dir_list.append(filedir)
                filedir_string = str(filedir)[:-len(filename)-5] # creates a file directory string

                filedir_images_string = filedir_string + r'\subcellular localisation' # creates a new directory

                if os.path.isdir(filedir_images_string) == True: #checks if the new directory allready exists
                    pass
                else:
                    os.makedirs(filedir_images_string) # creates the new directory

                for paths, dirs, files in sd.walk(filedir_string + '/csv'):  # goes throw all files and folders in given directory

                    for file in os.listdir(paths):  # goes throw all files in a folder
                        filedir = os.path.join(paths, file)
                        filedir_string = str(filedir) # creates a file directory string

                        if 'RunCellpose_N' in filedir_string: # checks for a condition in the other files
                            dir_list.append(filedir)
                        if 'Cytoplasm' in filedir_string: # checks for a condition in the other files
                            dir_list.append(filedir)

                        if len(dir_list) == 3: # if three files are found the analisis is done

                            print(dir_list)

                            Full_cell = pd.read_csv(dir_list[0]) # reads the found files
                            Cytoplasm = pd.read_csv(dir_list[1])
                            Nucleus = pd.read_csv(dir_list[2])

                            single_cells_and_organells = cell_to_organell_basic(Full_cell, Cytoplasm, Nucleus,1) # matches subcellular locations
                            single_cells_and_organells.to_csv(filedir_images_string + '/' + filename_string[:-18] +'_single_cells_and_organells.csv') # exports the file


for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        file_list = []
        if filedir.endswith("single_cells_and_organells.csv"): # returns all files that end with txt

            file_list.append(filedir)

            filename = os.path.basename(file)
            filename_string = str(filename)
            filedir_string = str(filedir)[:-len(filename)-1]

            filedir_images_string = filedir_string[:-len('subcellular localisation/')] + r'\subcellular localisation and neigboorhood'

            if os.path.isdir(filedir_images_string) == True:
                pass
            else:
                os.makedirs(filedir_images_string)

        for file in file_list:
            print(file)

            filename = os.path.basename(file)
            filename_string = str(filename)
            filedir_string = str(filedir)[:-len(filename) - 1]
            filedir_images_string = filedir_string[:-len('subcellular localisation/')] + r'\subcellular localisation and neigboorhood'

            df = pd.read_csv(filedir)
            df = df.drop('Unnamed: 0', axis=1)

            neigboorhood_sub = neigboorhood(df)
            neigboorhood_sub.to_csv(filedir_images_string + '/' + filename_string[:-len('single_cells_and_organells')-5] + '_subcell_neighborhood.csv')

'''
for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        dir_list = []
        if filedir.endswith("RunCellpose_C.csv"): # returns all files that end with txt

            filename = os.path.basename(file)
            filename_string = str(filename)

            filedir_string = str(filedir)[:-len(filename) - 5]
            filedir_images_string = filedir_string + r'\subcellular localisation'

            dir_list.append(filedir)

            if os.path.isdir(filedir_images_string) == True:
                pass
            else:
                os.makedirs(filedir_images_string)

            for paths, dirs, files in sd.walk(filedir_string + "/csv"):  # goes throw all files and folders in given directory

                for file in os.listdir(paths):  # goes throw all files in a folder
                    filedir = os.path.join(paths, file)  # returns full file directory

                    if filedir.endswith("RunCellpose_N.csv"):
                        dir_list.append(filedir)
                    if filedir.endswith("Cytoplasm.csv"):
                        dir_list.append(filedir)
                    if filedir.endswith("Outline_cell_1.tiff"):

                        filename = os.path.basename(file)
                        filename_string = str(filename)
                        filedir_string = str(filedir)[:-len(filename)-1]

                        dir_list.append(filedir_string)

                    if filedir.endswith("Outline_nucleus_1.tiff"):

                        filename = os.path.basename(file)
                        filename_string = str(filename)
                        filedir_string = str(filedir)[:-len(filename)]

                        dir_list.append(filedir_string)

            if len(dir_list) == 5:
                print(dir_list)

                Full_cell = pd.read_csv(dir_list[0])
                Cytoplasm = pd.read_csv(dir_list[1])
                Nucleus = pd.read_csv(dir_list[2])
                Cell_outlines = dir_list[3]
                Nucleus_outlines = dir_list[4]

                single_cells_and_organells = cell_to_organell_advanced(Full_cell, Cytoplasm, Nucleus,1,Cell_outlines,Nucleus_outlines)
                #single_cells_and_organells.to_csv(filedir_images_string + '/' + filename_string[:-18] +'_single_cells_and_organells.csv')
'''