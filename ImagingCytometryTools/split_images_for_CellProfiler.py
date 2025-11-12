import os
import scandir as sd
import numpy as np
from PIL import Image, ImageOps

def split_images_for_CellProfiler(data_direction, crop_count=4):

    if crop_count != 4:

        print('not added yet')

    else:

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

                    filedir_Images = dir_name_string + '/' + 'image crops'

                    if os.path.isdir(filedir_Images) == True:
                        pass

                    else:
                        os.makedirs(filedir_Images)

                    print('Generating image crops for: ' + str(filedir))

                    image = Image.open(filedir)
                    width, height = image.size

                    for x in range(crop_count):

                        if x == 0:
                            image_croped = image.crop(((0, 0, width/2, height/2)))
                            export_file_path_full = filedir_Images + '/' + filename_new_string + '_' + str(1) + '.tiff'
                            image_croped.save(export_file_path_full)

                        if x == 1:
                            image_croped = image.crop((((width/2)+1, 0, width, height/2)))
                            export_file_path_full = filedir_Images + '/' + filename_new_string + '_' + str(2) + '.tiff'
                            image_croped.save(export_file_path_full)

                        if x == 2:
                            image_croped = image.crop(((0, (height/2)+1, width/2, height)))
                            export_file_path_full = filedir_Images + '/' + filename_new_string + '_' + str(3) + '.tiff'
                            image_croped.save(export_file_path_full)

                        if x == 3:
                            image_croped = image.crop((((width/2)+1, (height/2) + 1, width, height)))
                            export_file_path_full = filedir_Images + '/' + filename_new_string + '_' + str(4) + '.tiff'
                            image_croped.save(export_file_path_full)


#ome tiff
'''
from PIL import Image
import tifffile

img = Image.open("input_image.ome.tiff")
box = (100, 100, 400, 400)
cropped = img.crop(box)

import numpy as np
cropped_array = np.array(cropped)

tifffile.imwrite(
    "cropped_output.ome.tiff",
    cropped_array,
    photometric='minisblack',
    metadata={'axes': 'YX'}  # adjust axes as needed (e.g., 'CYX', 'ZYX')
)
'''