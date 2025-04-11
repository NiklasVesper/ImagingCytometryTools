import pandas
import numpy as np

'''
A collection of functions for the analysis.
'''

#removes the non recorded channels from the spillover matrix
def remove_non_recorded_channels(mat,chan):

    channel_indexes = [mat.columns.get_loc(c) - 1 for c in chan if c in mat] #gets chanel indexes from the spillover matrix
    spill_mat_final = mat.iloc[: , 1:]

    for c in chan:
        spill_mat_final = spill_mat_final.drop(c, axis = 1) #drops the spillover matrix rows that are not recorded

    spill_mat_final = spill_mat_final.drop(index = channel_indexes) #drops the spillover matrix columns that are not recorded

    return(spill_mat_final)

#inverts the spillover matrix into the compensation matrix
def invert_spillover_matrix(spillover_matrix):
    spill_mat_final_np = spillover_matrix.to_numpy() #converts into numpy array
    comp_mat = np.linalg.inv(spill_mat_final_np) #inverts the spillover matrix

    return(comp_mat)