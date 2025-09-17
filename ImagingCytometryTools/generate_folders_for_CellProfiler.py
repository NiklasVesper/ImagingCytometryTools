import os
import scandir as sd

'''
This function automatically generates folders for the CellProfiler pipeline.

Here two variations are possible depending in how you want to analyse your data:

Generate the folders per .mcd file to group together multiple images into one output .csv file in the end,
or one filer per image group to analyze each image separately.
'''

def generate_folders_for_CellProfiler(data_direction, generate_folders_per_mcd_file=True):

    for paths, dirs, files in sd.walk(data_direction):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".mcd"):

                if generate_folders_per_mcd_file == True:

                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    filedir_CP_csv = dir_name_string + '/' + 'csv and outlines' + '/' + 'csv'
                    filedir_CP_outline_cell = dir_name_string + '/' + 'csv and outlines' + '/' + 'outline cell'
                    filedir_CP_outline_nucleus = dir_name_string + '/' + 'csv and outlines' + '/' + 'outline nucleus'

                    if os.path.isdir(filedir_CP_csv) == True:
                        pass
                    else:
                        os.makedirs(filedir_CP_csv)

                    if os.path.isdir(filedir_CP_outline_cell) == True:
                        pass
                    else:
                        os.makedirs(filedir_CP_outline_cell)

                    if os.path.isdir(filedir_CP_outline_nucleus) == True:
                        pass
                    else:
                        os.makedirs(filedir_CP_outline_nucleus)

                elif generate_folders_per_mcd_file == False:

                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    for paths, dirs, files in sd.walk(dir_name_string):

                        for file in os.listdir(paths):
                            filedir = os.path.join(paths, file)

                            if filedir.endswith("_compensated"):

                                filedir_CP_csv = filedir + '/' + 'csv and outlines' + '/' + 'csv'
                                filedir_CP_outline_cell = filedir + '/' + 'csv and outlines' + '/' + 'outline cell'
                                filedir_CP_outline_nucleus = filedir + '/' + 'csv and outlines' + '/' + 'outline nucleus'

                                if os.path.isdir(filedir_CP_csv) == True:
                                    pass
                                else:
                                    os.makedirs(filedir_CP_csv)

                                if os.path.isdir(filedir_CP_outline_cell) == True:
                                    pass
                                else:
                                    os.makedirs(filedir_CP_outline_cell)

                                if os.path.isdir(filedir_CP_outline_nucleus) == True:
                                    pass
                                else:
                                    os.makedirs(filedir_CP_outline_nucleus)

    print('I made some folders for you. ------------------------------------------------------------------------------')