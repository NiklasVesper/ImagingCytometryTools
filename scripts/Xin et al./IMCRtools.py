from imctools.io.txt.txtparser import TxtParser
import os
import scandir as sd

folder_dir = r'D:\Basel' # folder directory

#generates images
for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith("_comp.txt"): # returns all files that end with txt

            print(filedir)

            filename = os.path.basename(file) # gives you the file name
            filename_string = str(filename)
            filename_new = os.path.basename(file).strip('_comp.txt') #.strip('.txt') # returns filename without '.txt'
            filename_new_string = str(filename_new)

            filedir_string = str(filedir)[:-len(filename)] # gives you the file directory

            filedir_images_string = filedir_string + 'Images' # creates a new directory

            if os.path.isdir(filedir_images_string) == True: #checks if the new directory allready exists
                pass
            else:
                os.makedirs(filedir_images_string) # creates a new directory

            filedir_string_full = str(filedir)

            old_file_name = os.path.join(filedir_string[:-1], filename_string)
            new_file_name = os.path.join(filedir_string[:-1], filename_string.replace('_comp',''))

            if os.path.exists(new_file_name):
                os.remove(new_file_name)

            os.rename(old_file_name, new_file_name)

            parser = TxtParser(str(new_file_name)) # calls parser
            ac_data = parser.get_acquisition_data() # gets data
            ac_data.save_tiffs(filedir_images_string,basename = filename_new_string, compression=0) # saves tiffs 