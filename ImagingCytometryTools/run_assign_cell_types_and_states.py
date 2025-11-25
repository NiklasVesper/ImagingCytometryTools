import pandas as pd
import scandir as sd
import os

from assign_phenotypes_and_metadata import phenotype_row

'''
This function runs the functions that assigns each cell a phenotype based on the mean intensity.
'''

def run_assign_cell_types_and_states(directory, filestring, new_folder_name, overwrite_existing_files = False):

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

                    filedir_phenotyping = dir_name_string[:-len(folder_name_string)] + new_folder_name

                    if os.path.isdir(filedir_phenotyping) == True:
                        pass
                    else:
                        os.makedirs(filedir_phenotyping)

                    export_file_path = filedir_phenotyping + '/' + filename_string[:-len(filestring) - 4] + str(new_folder_name).replace(' ', '_') + '.csv'

                    if overwrite_existing_files == False:

                        if os.path.isfile(export_file_path) == True:
                            print('The file: ' + export_file_path + ' already exists!')
                            continue

                        else:
                            print('Running phenotypings for: ' + str(filedir))

                            phenotypes = pd.read_csv(filedir)
                            phenotype_row(phenotypes)
                            phenotypes.to_csv(export_file_path)

                    else:
                        print('Running phenotypings for: ' + str(filedir))

                        phenotypes = pd.read_csv(filedir)
                        phenotype_row(phenotypes)
                        phenotypes.to_csv(export_file_path)

    print('I phenotyped all the cells. -------------------------------------------------------------------------------')