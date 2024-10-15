from imctools.io.txt.txtparser import TxtParser
import os
import scandir as sd

'''
The "generate images" script makes images out of compensated IMC data (.mcd format).
The script takes only compensated images.

Here the images were compensated with the script available at https://github.com/BodenmillerGroup/cyTOFcompensation/tree/master.
The corresponding paper you find under: https://doi.org/10.1016/j.cels.2018.02.010
'''

folder_dir = r'D:\ATF6' #folder directory

for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns full file directory

        if filedir.endswith("_comp.txt"): #returns all files that end with "_comp.txt"

            print(filedir) #prints file directory

            filename = os.path.basename(file) #gives you the file name
            filename_string = str(filename) #turns the filename into a string
            filename_new = os.path.basename(file).strip('_comp.txt') #strips '_comp.txt' from the file name and creates a new one
            filename_new_string = str(filename_new) #turns the new filename into a string

            filedir_string = str(filedir)[:-len(filename)] #gives you the file directory as a string

            filedir_images_string = filedir_string + 'Images' #creates a new directory

            if os.path.isdir(filedir_images_string) == True: #checks if the new directory allready exists
                pass
            else: #creates a new directory if necessary
                os.makedirs(filedir_images_string)

            filedir_string_full = str(filedir) #turns the directory into a string

            old_file_name = os.path.join(filedir_string[:-1], filename_string) # creates the old file name
            new_file_name = os.path.join(filedir_string[:-1], filename_string.replace('_comp','')) # creates new file name

            if os.path.exists(new_file_name): #removes a duplicate file if necessary
                os.remove(new_file_name)

            os.rename(old_file_name, new_file_name) #renames the file so it is taken by the .txt parser

            parser = TxtParser(str(new_file_name)) # calls parser
            ac_data = parser.get_acquisition_data() # gets data
            ac_data.save_tiffs(filedir_images_string, basename = filename_new_string, compression=0) # saves tiffs
