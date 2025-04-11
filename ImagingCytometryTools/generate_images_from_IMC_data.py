from imctools.io.txt.txtparser import TxtParser
import os
import scandir as sd

'''
This function generates images out of IMC data (.mcd format).

Here you can either generate single .tiff files per sample or in addition to that .ometiffs.

In addition to that if you change the name back to the original .txt file so you can view the compensated .txt file in the MCD viewer.
Otherwise, just view them in Napari (1). But you have to be careful if you rerun the script because it will be recognized as an uncompensated file.

(1) Chiu, C. L., & Clack, N. (2022).
Napari: a Python multi-dimensional image viewer platform for the research community. Microscopy and Microanalysis, 28(S1), 1576-1577.
'''

#generates .tiff / .ometiff files out of IMC data
def generate_images_from_IMC_data(data_direction, generate_ome_tiff = False, change_name_for_MCD_viewer = False):

    for paths, dirs, files in sd.walk(data_direction): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith("_comp.txt"): #returns all files that end with _comp.txt

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                #creates a new file name
                filename_new = os.path.basename(file).strip('_comp.txt')
                filename_new_string = str(filename_new)

                #gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                filedir_Images = dir_name_string + '/' + 'compensated images' #creates a new directory for .tiff files
                filedir_Images_ome = dir_name_string + '/' + 'compensated images ome tiff' #creates a new directory for .ometiff files

                #checks if the new directory already exists
                if os.path.isdir(filedir_Images) == True:
                    pass
                else:
                    os.makedirs(filedir_Images) #creates a new directory if necessary

                #checks if the new directory already exists
                if generate_ome_tiff == True:
                    if os.path.isdir(filedir_Images_ome) == True:
                        pass
                    else:
                        os.makedirs(filedir_Images_ome) #creates a new directory if necessary

                old_file_name = os.path.join(dir_name_string, filename_string) #recreates old file directory
                new_file_name = os.path.join(dir_name_string, filename_string.replace('_comp', '')) #creates new file name in that directory

                print('Generating images for: ' + str(filedir))

                #checks if the new directory already exists
                if os.path.exists(new_file_name):
                    os.remove(new_file_name)

                os.rename(old_file_name, new_file_name) #renames the old file in the directory in the new one

                parser = TxtParser(str(new_file_name)) #creates a parser
                ac_data = parser.get_acquisition_data() #gets the data for the parser

                if generate_ome_tiff == True:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0) #generates .tiff files
                    ac_data.save_ome_tiff(filedir_Images_ome + '/' + filename_string[:-4] + '.ome.tiff') #generates .ometiff files

                elif generate_ome_tiff == False:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0) #generates .tiff files

                if change_name_for_MCD_viewer == False:
                    os.rename(new_file_name, old_file_name) #renames the new file in the directory in the old one

                elif change_name_for_MCD_viewer == True:
                    print('If you rerun the compensation now in the same folder this file will be compensated again, because it has the same name as the original file!')

            elif filedir.endswith(".txt"): #returns all files that end with .txt

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                #creates a new file name
                filename_new = os.path.basename(file).strip('.txt')
                filename_new_string = str(filename_new)

                # gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                filedir_Images = dir_name_string + '/' + str(filename_string[:-4]) + '/' + 'images' #creates a new directory for .tiff files
                filedir_Images_ome = dir_name_string + '/' + str(filename_string[:-4]) + '/' + 'images ome tiff' #creates a new directory for .ometiff files

                #checks if the new directory already exists
                if os.path.isdir(filedir_Images) == True:
                    pass
                else:
                    os.makedirs(filedir_Images) #creates a new directory if necessary

                #checks if the new directory already exists
                if generate_ome_tiff == True:
                    if os.path.isdir(filedir_Images_ome) == True:
                        pass
                    else:
                        os.makedirs(filedir_Images_ome) #creates a new directory if necessary

                file_name = os.path.join(dir_name_string, filename_string) #recreates file directory

                print('Generating images for: ' + str(filedir))

                parser = TxtParser(str(file_name)) #creates a parser
                ac_data = parser.get_acquisition_data() #gets the data for the parser

                if generate_ome_tiff == True:
                    ac_data.save_tiffs(filedir_Images, basename=filename_new_string) #generates .tiff files
                    ac_data.save_ome_tiff(filedir_Images_ome + '/' + filename_string[:-4] + '.ome.tiff') #generates .ometiff files

                elif generate_ome_tiff == False:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0) #generates .tiff files

    print('I made pretty images for you. -----------------------------------------------------------------------------')