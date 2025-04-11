import pandas as pd
import scandir as sd
import os

from assign_stuff import phenotype_row

'''
This function runs the functions that assigns each cell a phenotype based on the mean intensity.
'''

def run_assign_cell_types(directory, filestring, new_folder_name, overwrite_existing_files = False):

    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns fuull file directory

            if filedir.endswith(".csv"): # returns all files that end with txt

                # gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                # gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                # gives you the folder name
                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string:

                    filedir_phenotyping = dir_name_string[:-len(folder_name_string)] + new_folder_name #creates a new directory

                    # checks if the new directory already exists
                    if os.path.isdir(filedir_phenotyping) == True:
                        pass
                    else:
                        os.makedirs(filedir_phenotyping)

                    export_file_path = filedir_phenotyping + '/' + filename_string[:-len(filestring) - 4] + str(new_folder_name).replace(' ', '_') + '.csv' #creates a export file path

                    if overwrite_existing_files == False:

                        # checks if the new file directory already exists
                        if os.path.isfile(export_file_path) == True:
                            print('The file: ' + export_file_path + ' already exists!')
                            continue

                        else:
                            print('Running phenotypings for: ' + str(filedir))

                            phenotypes = pd.read_csv(filedir) #imports the data frame
                            phenotype_row(phenotypes) #assigns each row a phenotype
                            phenotypes.to_csv(export_file_path) #exports the file

                    else:
                        print('Running phenotypings for: ' + str(filedir))

                        phenotypes = pd.read_csv(filedir)  # imports the data frame
                        phenotype_row(phenotypes)  # assigns each row a phenotype
                        phenotypes.to_csv(export_file_path)  # exports the file

    print('I phenotyped all the cells. -------------------------------------------------------------------------------')