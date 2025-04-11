from numba import njit
import numpy as np

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

#for GPU accelerated compensation
@njit
def compensate_row_GPU(data_row, final_compensation_matrix):
    compenstated_row = []  # empty list to store the compensated row

    datapoint_counter = 0  # counter to know which part of the compensation matrix needs to be called to compensate the data point
    for x in range(len(data_row)):  # goes throw each data point in the row
        compenstated_datapoint = np.multiply(data_row,final_compensation_matrix[:, datapoint_counter])
        compenstated_row.append(np.sum(compenstated_datapoint))  # compensates a single data point
        datapoint_counter = datapoint_counter + 1  # ticks up the data point counter to move to the next row in the compensation matrix

    compenstated_row = [max(0, x) for x in compenstated_row]

    comp_val_counter = 0
    for comp_val in compenstated_row:
        if comp_val > data_row[comp_val_counter]:
            compenstated_row[comp_val_counter] = data_row[comp_val_counter]
        comp_val_counter = comp_val_counter + 1

    return(compenstated_row)

#for "normal" compensation
def compensate_row_no_GPU(data_row, final_compensation_matrix):
    compenstated_row = []  # empty list to store the compensated row

    datapoint_counter = 0  # counter to know which part of the compensation matrix needs to be called to compensate the data point
    for x in range(len(data_row)):  # goes throw each data point in the row
        compenstated_datapoint = np.multiply(data_row,final_compensation_matrix[:, datapoint_counter])
        compenstated_row.append(np.sum(compenstated_datapoint))  # compensates a single data point
        datapoint_counter = datapoint_counter + 1  # ticks up the data point counter to move to the next row in the compensation matrix

    compenstated_row = [max(0, x) for x in compenstated_row]

    comp_val_counter = 0
    for comp_val in compenstated_row:
        if comp_val > data_row[comp_val_counter]:
            compenstated_row[comp_val_counter] = data_row[comp_val_counter]
        comp_val_counter = comp_val_counter + 1

    return(compenstated_row)