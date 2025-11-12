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

def generate_images_from_IMC_data(data_direction, generate_ome_tiff = False, change_name_for_MCD_viewer = False):

    for paths, dirs, files in sd.walk(data_direction):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith("_comp.txt"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                filename_new = os.path.basename(file).strip('_comp.txt')
                filename_new_string = str(filename_new)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                filedir_Images = dir_name_string + '/' + 'compensated images'
                filedir_Images_ome = dir_name_string + '/' + 'compensated images ome tiff'

                if os.path.isdir(filedir_Images) == True:
                    pass
                else:
                    os.makedirs(filedir_Images)

                if generate_ome_tiff == True:
                    if os.path.isdir(filedir_Images_ome) == True:
                        pass
                    else:
                        os.makedirs(filedir_Images_ome)

                old_file_name = os.path.join(dir_name_string, filename_string)
                new_file_name = os.path.join(dir_name_string, filename_string.replace('_comp', ''))

                print('Generating images for: ' + str(filedir))

                if os.path.exists(new_file_name):
                    os.remove(new_file_name)

                os.rename(old_file_name, new_file_name)

                parser = TxtParser(str(new_file_name))
                ac_data = parser.get_acquisition_data()

                if generate_ome_tiff == True:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0)
                    ac_data.save_ome_tiff(filedir_Images_ome + '/' + filename_string[:-4] + '.ome.tiff')

                elif generate_ome_tiff == False:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0)

                if change_name_for_MCD_viewer == False:
                    os.rename(new_file_name, old_file_name)

                elif change_name_for_MCD_viewer == True:
                    print('If you rerun the compensation now in the same folder this file will be compensated again, because it has the same name as the original file!')

            elif filedir.endswith(".txt"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                filename_new = os.path.basename(file).strip('.txt')
                filename_new_string = str(filename_new)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                filedir_Images = dir_name_string + '/' + str(filename_string[:-4]) + '/' + 'images'
                filedir_Images_ome = dir_name_string + '/' + str(filename_string[:-4]) + '/' + 'images ome tiff'

                if os.path.isdir(filedir_Images) == True:
                    pass
                else:
                    os.makedirs(filedir_Images)

                if generate_ome_tiff == True:

                    if os.path.isdir(filedir_Images_ome) == True:
                        pass
                    else:
                        os.makedirs(filedir_Images_ome)

                file_name = os.path.join(dir_name_string, filename_string)

                print('Generating images for: ' + str(filedir))

                parser = TxtParser(str(file_name))
                ac_data = parser.get_acquisition_data()

                if generate_ome_tiff == True:
                    ac_data.save_tiffs(filedir_Images, basename=filename_new_string)
                    ac_data.save_ome_tiff(filedir_Images_ome + '/' + filename_string[:-4] + '.ome.tiff')

                elif generate_ome_tiff == False:
                    ac_data.save_tiffs(filedir_Images,basename = filename_new_string, compression=0)

    print('I made pretty images for you. -----------------------------------------------------------------------------')