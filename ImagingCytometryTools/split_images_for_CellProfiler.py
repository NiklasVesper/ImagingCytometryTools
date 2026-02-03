import os
import scandir as sd
import numpy as np
from PIL import Image, ImageOps

def split_images_for_CellProfiler(data_direction, crop_count=4, image_size_threshold=600000):

    if crop_count == 4:

        for paths, dirs, files in sd.walk(data_direction):

            for file in os.listdir(paths):
                filedir = os.path.join(paths, file)

                filename = os.path.basename(file)
                filename_string = str(filename)

                if filename.endswith(".tiff"):

                    filename_new = os.path.basename(file).strip('.tiff')
                    filename_new_string = str(filename_new)

                    dir_name = os.path.dirname(filedir)
                    dir_name_string = str(dir_name)

                    image = Image.open(filedir)
                    width, height = image.size

                    if width * height > image_size_threshold:
                        
                        print('Generating image crops for: ' + str(filedir))

                        for x in range(crop_count):

                            if x == 0:
                                image_croped = image.crop(((0, 0, width/2, height/2)))

                            if x == 1:
                                image_croped = image.crop((((width / 2) + 1, 0, width, height / 2)))

                            if x == 2:
                                image_croped = image.crop(((0, (height / 2) + 1, width / 2, height)))

                            if x == 3:
                                image_croped = image.crop((((width / 2) + 1, (height / 2) + 1, width, height)))

                            filedir_Images = dir_name_string + '/' + 'image crops ' + str(x+1)

                            if os.path.isdir(filedir_Images) == True:
                                pass

                            else:
                                os.makedirs(filedir_Images)

                            export_file_path_full = filedir_Images + '/' + filename_new_string + '_' + str(x+1) + '.tiff'
                            image_croped.save(export_file_path_full)

                        os.remove(filedir)

    else:
        print('Currently, you can only choose 4 crops because I am lazy.')


#ome tiff
'''
from PIL import Image
import tifffile
import numpy as np
img = Image.open("input_image.ome.tiff")
box = (100, 100, 400, 400)
cropped = img.crop(box)
cropped_array = np.array(cropped)
tifffile.imwrite("cropped_output.ome.tiff",cropped_array,photometric='minisblack',
metadata={'axes': 'YX'})
'''