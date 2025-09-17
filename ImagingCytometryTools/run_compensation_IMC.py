import os
import numpy as np
import pandas as pd
import scandir as sd
from tqdm import tqdm

from ImagingCytometryTools.get_data_from_files import get_metals_from_IMC_data
from ImagingCytometryTools.get_data_from_files import get_spillover_matrix
from ImagingCytometryTools.get_data_from_files import get_file_for_compensation

from ImagingCytometryTools.helper_functions import remove_non_recorded_channels
from ImagingCytometryTools.helper_functions import invert_spillover_matrix

from ImagingCytometryTools.compensation import compensate_row_no_GPU
from ImagingCytometryTools.compensation import compensate_row_GPU

'''
Compensation is performed by generating the inverse spillover matrix to generate the compensation matrix.
Then, each data point row is multiplied with the compensation matrix to generate the compensated data point (1).
The compensation of images is based on the theoretical foundations of the compensation of flow cytometry data (2,3).

If you want to use your GPU, you need to find and install the CUDA version that is compatible with your Nvidia graphics card.
CUDA: https://www.nvidia.com/en-us/

(1) Roederer, M. (2002).
Compensation in flow cytometry. Current protocols in cytometry, 22(1), 1-14.

(2) Chevrier, S., Crowell, H. L., Zanotelli, V. R., Engler, S., Robinson, M. D., & Bodenmiller, B. (2018). 
Compensation of signal spillover in suspension and imaging mass cytometry. Cell systems, 6(5), 612-620.

(3)Novo, D., Grégori, G., & Rajwa, B. (2013). 
Generalized unmixing model for multispectral flow cytometry utilizing nonsquare compensation matrices. Cytometry Part A, 83(5), 508-520.
'''

def run_compensation_IMC(data_direction, spillover_matrix, use_GPU = False):

    for paths, dirs, files in sd.walk(data_direction):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith("comp.txt"):
                continue

            elif filedir.endswith(".txt"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                filedir_comp = dir_name_string + '/' + 'compensated data' + '/' + filename_string[:-4] + '_compensated' #creates a new directory

                if os.path.isdir(filedir_comp) == True:
                    pass
                else:
                    os.makedirs(filedir_comp)

                export_file_path = filedir_comp + '/' + filename_string[:-4] + '_comp.txt'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                    continue

                else:
                    print('I am compensating: ' + str(filedir))

                    channels = get_metals_from_IMC_data(filedir)

                    spill_mat = get_spillover_matrix(spillover_matrix)
                    comp_mat_channels = spill_mat.columns.values.tolist()
                    non_recorded_channels = list(set(channels).symmetric_difference(set(comp_mat_channels)))
                    comp_mat_for_use = remove_non_recorded_channels(spill_mat,non_recorded_channels)
                    comp_mat = invert_spillover_matrix(comp_mat_for_use)

                    data_for_comp = get_file_for_compensation(filedir)
                    pixel_positions = data_for_comp.iloc[:, :6]
                    pixel_values = data_for_comp.iloc[:, 6:]
                    pixel_values_np = pixel_values.to_numpy()

                    comp_data_np = []

                    for row in tqdm(pixel_values_np):

                        if use_GPU == False:
                            compensated_row = compensate_row_no_GPU(row, comp_mat)
                            comp_data_np.append(compensated_row)

                        elif use_GPU == True:
                            compensated_row = compensate_row_GPU(row, comp_mat)
                            comp_data_np.append(compensated_row)

                    comp_data = pd.DataFrame(comp_data_np, columns = pixel_values.columns.values.tolist())
                    comp_data_final = pd.concat([pixel_positions,comp_data], axis=1)
                    comp_data_final.to_csv(export_file_path, sep='\t', index=False)

    print('I am done compensating. -----------------------------------------------------------------------------------')