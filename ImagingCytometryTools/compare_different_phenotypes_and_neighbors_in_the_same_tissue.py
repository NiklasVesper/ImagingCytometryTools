import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import scandir as sd
import os
import pandas as pd
import ast
import seaborn as sns
import statistics
import numpy as np
from scipy.stats import mannwhitneyu

'''
Runs the statistics on the different phenotypes and cell states.

Here, one can select between calculating the overall abundance or cell category (immune/tissue cell) of a phenotype.
In addition to that, the overall cellular state can be accessed, with or without specific neighboring cells.
'''

def compare_different_phenotypes_and_neighbors_in_the_same_tissue(directory, filestring, df_column, output_folder, select_cell_type_and_state, select_files = [False], add_mixed_cells=False, analyse_neighboring_cells=False, select_neighboring_cell_type_and_state=[False, ['', [], []],['', [], []]], compare_means=[False,[]], split_by='file'):

    file_counter = 0
    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string and file_counter == 0:

                    file_counter = file_counter + 1

                    Cells = pd.read_csv(filedir)

                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]

                    cell_type_state_and_neighbor_1 = Cells.iloc[0:0]
                    cell_type_state_and_neighbor_2 = Cells.iloc[0:0]
                    neighboring_cell_type_state_and_neighbor_1 = Cells.iloc[0:0]
                    neighboring_cell_type_state_and_neighbor_2 = Cells.iloc[0:0]

                    cell_type_state_and_neighbor_1['Selected_cell_counter'] = ''
                    cell_type_state_and_neighbor_2['Selected_cell_counter'] = ''
                    neighboring_cell_type_state_and_neighbor_1['Selected_cell_counter'] = ''
                    neighboring_cell_type_state_and_neighbor_2['Selected_cell_counter'] = ''

                    if split_by == 'image':

                        cell_type_state_and_neighbor_1['Individual_image_counter'] = ''
                        cell_type_state_and_neighbor_2['Individual_image_counter'] = ''
                        neighboring_cell_type_state_and_neighbor_1['Individual_image_counter'] = ''
                        neighboring_cell_type_state_and_neighbor_2['Individual_image_counter'] = ''

                elif file_counter > 0:
                    continue

    Selected_cell_counter_1 = -1
    Selected_cell_counter_2 = -1
    for paths, dirs, files in sd.walk(directory):

        for file in os.listdir(paths):
            filedir = os.path.join(paths, file)

            if filedir.endswith(".csv"):

                filename = os.path.basename(file)
                filename_string = str(filename)

                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string and select_files[0] == False:

                    print('Processing: ' + filedir)

                    Cells = pd.read_csv(filedir)

                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]

                    UniqueImageNumber = Cells.ImageNumber.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                    for number in UniqueImageNumber:

                        Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                        for index, cell in Cells_on_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == select_cell_type_and_state[0][0] and set(select_cell_type_and_state[0][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[0][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Selected_cell_counter_1 = Selected_cell_counter_1 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                if split_by == 'image':
                                                    neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:
    
                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):
    
                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])
    
                                                for neighboring_type_and_state in neighboring_cell_types_and_states:
    
                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary
                                                        
                                        if select_neighboring_cell_type_and_state[0] == False:
                                            
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                    if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Selected_cell_counter_2 = Selected_cell_counter_2 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):
                                                    
                                                neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                if split_by == 'image':
                                                    neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[2][0] and set(select_neighboring_cell_type_and_state[2][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[2][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False:
                                            
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[0][0] and set(select_cell_type_and_state[0][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[0][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        Selected_cell_counter_1 = Selected_cell_counter_1 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                if len(ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])) == 2:

                                                    neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()                                                   
                                                    neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                    if split_by == 'image':
                                                        neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                    neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary
                                                    
                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:
                                            
                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        
                                                        sub_counter = sub_counter + 1
                                                        
                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary
                                                        
                                                        if sub_counter == 0: 
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == False:
                                            
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                        
                                        Selected_cell_counter_2 = Selected_cell_counter_2 + 1
                                        
                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                if len(ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])) == 2:

                                                    neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()                                                   
                                                    neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                    if split_by == 'image':
                                                        neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                    neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary
                                        
                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[2][0] and set(select_neighboring_cell_type_and_state[2][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[2][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        
                                                        sub_counter = sub_counter + 1
                                                        
                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary
                                                        
                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == False:
                                            
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                if filestring in filename_string and select_files[0] == True and filename in select_files[1]:

                    print('Processing: ' + filedir)

                    Cells = pd.read_csv(filedir)

                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]

                    UniqueImageNumber = Cells.ImageNumber.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImageNumber}

                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                    for number in UniqueImageNumber:

                        Cells_on_image = pd.DataFrame(DataFrameDict_Full_cell[number])

                        for index, cell in Cells_on_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == select_cell_type_and_state[0][0] and set(select_cell_type_and_state[0][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[0][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Selected_cell_counter_1 = Selected_cell_counter_1 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):
                                                neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                if split_by == 'image':
                                                    neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False:
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                    if type_and_state[0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(type_and_state[1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Selected_cell_counter_2 = Selected_cell_counter_2 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):
                                                neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                if split_by == 'image':
                                                    neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[2][0] and set(select_neighboring_cell_type_and_state[2][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[2][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False:
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[0][0] and set(select_cell_type_and_state[0][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[0][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        Selected_cell_counter_1 = Selected_cell_counter_1 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                if len(ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])) == 2:
                                                    neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                    neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                    if split_by == 'image':
                                                        neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                    neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_1.loc[len(neighboring_cell_type_state_and_neighbor_1)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == False:
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_1

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_1.loc[len(cell_type_state_and_neighbor_1)] = cell_dictionary

                                    if cell_types_and_states[0][0] == select_cell_type_and_state[1][0] and set(select_cell_type_and_state[1][1]).issubset(set(cell_types_and_states[0][1:])) is True and set(select_cell_type_and_state[1][2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        Selected_cell_counter_2 = Selected_cell_counter_2 + 1

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == False:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                if len(ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])) == 2:
                                                    neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                    neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                    if split_by == 'image':
                                                        neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                    neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                        if analyse_neighboring_cells == True and select_neighboring_cell_type_and_state[0] == True:

                                            sub_counter = -1
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_on_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[2][0] and set(select_neighboring_cell_type_and_state[2][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[2][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                        sub_counter = sub_counter + 1

                                                        neighboring_cell_dictionary = Cells_on_image.iloc[neighboring_cell].to_dict()
                                                        neighboring_cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                        if split_by == 'image':
                                                            neighboring_cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                        neighboring_cell_type_state_and_neighbor_2.loc[len(neighboring_cell_type_state_and_neighbor_2)] = neighboring_cell_dictionary

                                                        if sub_counter == 0:
                                                            cell_dictionary = cell.to_dict()
                                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                                            if split_by == 'image':
                                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == False:
                                            cell_dictionary = cell.to_dict()
                                            cell_dictionary['Selected_cell_counter'] = Selected_cell_counter_2

                                            if split_by == 'image':
                                                cell_dictionary['Individual_image_counter'] = filename_string + '_image_' + str(number)

                                            cell_type_state_and_neighbor_2.loc[len(cell_type_state_and_neighbor_2)] = cell_dictionary

    final_output_folder = output_folder + '/' + 'Comparing ' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + ' and ' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2])

    if split_by == 'file':

        if os.path.isdir(final_output_folder.replace('+', '') + '/' + 'Percentages per file') == True:
            pass
        else:
            os.makedirs(final_output_folder.replace('+', '') + '/' + 'Percentages per file')

    if split_by == 'image':

        if os.path.isdir(final_output_folder.replace('+', '') + '/' + 'Percentages per image') == True:
            pass
        else:
            os.makedirs(final_output_folder.replace('+', '') + '/' + 'Percentages per image')

    if split_by == 'file' and compare_means[0] == True:

        if os.path.isdir(final_output_folder.replace('+', '') + '/' + 'Mean per file') == True:
            pass
        else:
            os.makedirs(final_output_folder.replace('+', '') + '/' + 'Mean per file')

    if split_by == 'image' and compare_means[0] == True:

        if os.path.isdir(final_output_folder.replace('+', '') + '/' + 'Mean per image') == True:
            pass
        else:
            os.makedirs(final_output_folder.replace('+', '') + '/' + 'Mean per image')

    if split_by == 'cell':

        if os.path.isdir(final_output_folder.replace('+', '') + '/' + 'Mean per cell') == True:
            pass
        else:
            os.makedirs(final_output_folder.replace('+', '') + '/' + 'Mean per cell')

    if analyse_neighboring_cells == True:

        if split_by == 'file':

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file')

        if split_by == 'image':

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image')

        if split_by == 'cell':

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell')

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell')

        if split_by == 'file' and compare_means[0] == True:

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file')

        if split_by == 'image' and compare_means[0] == True:

            if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image') == True:
                pass
            else:
                os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image')

    export_file_path_1 = final_output_folder.replace('+', '') + '/' + str(select_cell_type_and_state[0][0]) + ' cells positive for ' + str(select_cell_type_and_state[0][1]) + ' negative for ' + str(select_cell_type_and_state[0][2]) + '.csv'
    export_file_path_2 = final_output_folder.replace('+', '') + '/' + str(select_cell_type_and_state[1][0]) + ' cells positive for ' + str(select_cell_type_and_state[1][1]) + ' negative for ' + str(select_cell_type_and_state[1][2]) + '.csv'
    neighboring_export_file_path_1 = final_output_folder.replace('+', '') + '/' + 'Cells neighboring ' + str(select_cell_type_and_state[0][0]) + ' cells positive for ' + str(select_cell_type_and_state[0][1]) + ' negative for ' + str(select_cell_type_and_state[0][2]) + '.csv'
    neighboring_export_file_path_2 = final_output_folder.replace('+', '') + '/' + 'Cells neighboring ' + str(select_cell_type_and_state[1][0]) + ' cells positive for ' + str(select_cell_type_and_state[1][1]) + ' negative for ' + str(select_cell_type_and_state[1][2]) + '.csv'

    if analyse_neighboring_cells == False:

        cell_type_state_and_neighbor_1.to_csv(export_file_path_1)
        cell_type_state_and_neighbor_2.to_csv(export_file_path_2)

    if analyse_neighboring_cells == True:

        cell_type_state_and_neighbor_1.to_csv(export_file_path_1)
        cell_type_state_and_neighbor_2.to_csv(export_file_path_2)
        neighboring_cell_type_state_and_neighbor_1.to_csv(neighboring_export_file_path_1)
        neighboring_cell_type_state_and_neighbor_2.to_csv(neighboring_export_file_path_2)

    if split_by == 'file':

        state_list = []

        for row in cell_type_state_and_neighbor_1[df_column]:
            cell_types_and_states = ast.literal_eval(row)

            for type_and_state in cell_types_and_states:

                if type_and_state[0] == select_cell_type_and_state[0][0]:

                    for state in type_and_state[1:]:
                        state_list.append(state)

        for row in cell_type_state_and_neighbor_2[df_column]:
            cell_types_and_states = ast.literal_eval(row)

            for type_and_state in cell_types_and_states:

                if type_and_state[0] == select_cell_type_and_state[1][0]:

                    for state in type_and_state[1:]:
                        state_list.append(state)

        state_list = list(set(state_list))

        UniqueFileName = cell_type_state_and_neighbor_1.ImageName.unique()
        DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
        for key in DataFrameDict_Full_cell.keys():
            DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_1[:][cell_type_state_and_neighbor_1.ImageName == key].reset_index()

        state_list_percent_1_final = []

        for file in UniqueFileName:

            Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

            state_list_counter_cells_1 = [0 for x in range(0, len(state_list))]

            for row in Cells_in_file[df_column]:

                cell_types_and_states = ast.literal_eval(row)
                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                for type_and_state in cell_types_and_states:

                    if type_and_state[0] == select_cell_type_and_state[0][0]:

                        for state in type_and_state[1:]:
                            state_list_counter_cell[state_list.index(state)] = 1

                state_list_counter_cells_1 = [state_list_counter_cells_1[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

            state_list_percent_1 = [(state_count / Cells_in_file.shape[0]) * 100 for state_count in state_list_counter_cells_1]
            state_list_percent_1_final.append(state_list_percent_1)

        UniqueFileName = cell_type_state_and_neighbor_2.ImageName.unique()
        DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
        for key in DataFrameDict_Full_cell.keys():
            DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_2[:][cell_type_state_and_neighbor_2.ImageName == key].reset_index()

        state_list_percent_2_final = []

        for file in UniqueFileName:

            Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

            state_list_counter_cells_2 = [0 for x in range(0, len(state_list))]

            for row in Cells_in_file[df_column]:

                cell_types_and_states = ast.literal_eval(row)
                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                for type_and_state in cell_types_and_states:

                    if type_and_state[0] == select_cell_type_and_state[0][0]:

                        for state in type_and_state[1:]:
                            state_list_counter_cell[state_list.index(state)] = 1

                state_list_counter_cells_2 = [state_list_counter_cells_2[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

            state_list_percent_2 = [(state_count / Cells_in_file.shape[0]) * 100 for state_count in state_list_counter_cells_2]

            state_list_percent_2_final.append(state_list_percent_2)

        for x in range(0, len(state_list)):

            if x == 0:

                plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                plt.xlim(0, 0.5)
                plt.ylim(0, 0.5)

                export_file_path = final_output_folder.replace('+', '') + '/' + 'Percentages per file/' + 'p_value legend.png'
                plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

            list_1 = [element[x] for element in state_list_percent_1_final]
            list_2 = [element[x] for element in state_list_percent_2_final]

            plot_list = [list_1, list_2]

            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list, widths=0.4)

            if 0.005 < p < 0.05:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
            elif 0.0005 < p < 0.005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
            elif 0.00005 < p < 0.0005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
            elif p < 0.00005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
            axs.set_title(str(state_list[x]))
            axs.set_ylabel('Percent ' + str(state_list[x]))
            axs.set_xticks([1, 2],
                           [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                            str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                           rotation=-90)

            export_file_path = final_output_folder.replace('+', '') + '/' + 'Percentages per file/' + str(state_list[x]) + '.png'
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
            plt.close()

        if compare_means[0] == True:

            if not compare_means[1]:

                mean_cols = [col for col in cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

            else:

                mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

            plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
            plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
            plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
            plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
            plt.xlim(0, 0.5)
            plt.ylim(0, 0.5)

            export_file_path = final_output_folder.replace('+', '') + '/' + 'Mean per file/' + 'p_value legend.png'
            plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

            for mean_col in mean_cols:

                mean_list_1 = []
                mean_list_2 = []

                UniqueFileName = cell_type_state_and_neighbor_1.ImageName.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_1[:][cell_type_state_and_neighbor_1.ImageName == key].reset_index()

                for file in UniqueFileName:
                    mean_per_file = pd.DataFrame(DataFrameDict_Full_cell[file]).loc[:, mean_col].mean()
                    mean_list_1.append(mean_per_file)

                UniqueFileName = cell_type_state_and_neighbor_2.ImageName.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_2[:][cell_type_state_and_neighbor_2.ImageName == key].reset_index()

                for file in UniqueFileName:
                    mean_per_file = pd.DataFrame(DataFrameDict_Full_cell[file]).loc[:, mean_col].mean()
                    mean_list_2.append(mean_per_file)

                plot_list = [mean_list_1, mean_list_2]

                U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                max_plot = max([max(plot) for plot in plot_list])
                min_plot = min([min(plot) for plot in plot_list])

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                axs.boxplot(plot_list, widths=0.4)

                if 0.005 < p < 0.05:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                elif 0.0005 < p < 0.005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                elif 0.00005 < p < 0.0005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                elif p < 0.00005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                axs.set_ylabel(str(mean_col))
                axs.set_title(str(mean_col))
                axs.set_xticks([1, 2],
                               [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                               rotation=-90)

                export_file_path = final_output_folder.replace('+', '') + '/' + 'Mean per file/' + str(mean_col) + '.png'
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                plt.close()

        if analyse_neighboring_cells == True:

            if select_neighboring_cell_type_and_state[1][0] != select_neighboring_cell_type_and_state[2][0]:

                print('Different neighboring cell lineages are selected. Since they can not be compared directly only the primary selected cell lines are compared.')

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            if select_neighboring_cell_type_and_state[0] == False:

                neighboring_state_list_heat_map = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                neighboring_state_list_heat_map = list(set(neighboring_state_list_heat_map))

                data_1 = {}
                data_2 = {}

                for state in neighboring_state_list_heat_map:
                    data_1[state] = []
                    data_2[state] = []

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + 'selected neighbors/' + str(lineage))

                neighboring_state_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                neighboring_state_list = list(set(neighboring_state_list))

                if select_neighboring_cell_type_and_state[0] == False:

                    non_contained_markers = list(set(neighboring_state_list).symmetric_difference(set(neighboring_state_list_heat_map)))

                    for non_contained_marker in non_contained_markers:
                        data_1[non_contained_marker].append(np.nan)
                        data_2[non_contained_marker].append(np.nan)

                UniqueFileName = neighboring_cell_type_state_and_neighbor_1.ImageName.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.ImageName == key].reset_index()

                state_list_percent_1_final = []
                cell_lineage_percent_1 = []

                for file in UniqueFileName:

                    Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

                    state_list_counter_cells_1 = [0 for x in range(0, len(neighboring_state_list))]

                    lineage_counter = 0
                    for row in Cells_in_file[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1

                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_1 = [state_list_counter_cells_1[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_1 = [(state_count / Cells_in_file.shape[0]) * 100 for state_count in state_list_counter_cells_1]

                    state_list_percent_1_final.append(state_list_percent_1)
                    cell_lineage_percent_1.append((lineage_counter / Cells_in_file.shape[0]) * 100)

                UniqueFileName = neighboring_cell_type_state_and_neighbor_2.ImageName.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.ImageName == key].reset_index()

                state_list_percent_2_final = []
                cell_lineage_percent_2 = []

                for file in UniqueFileName:

                    Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

                    state_list_counter_cells_2 = [0 for x in range(0, len(neighboring_state_list))]

                    lineage_counter = 0
                    for row in Cells_in_file[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1

                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_2 = [state_list_counter_cells_2[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_2 = [(state_count / Cells_in_file.shape[0]) * 100 for state_count in state_list_counter_cells_2]

                    state_list_percent_2_final.append(state_list_percent_2)
                    cell_lineage_percent_2.append((lineage_counter / Cells_in_file.shape[0]) * 100)

                if select_neighboring_cell_type_and_state[0] == False:

                    plot_list = [cell_lineage_percent_1, cell_lineage_percent_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage))
                    axs.set_ylabel('Percent ' + str(lineage) + ' cells')
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(lineage) + '/' + str(lineage) + '.png'
                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                for x in range(0, len(neighboring_state_list)):

                    if x == 0:
                        plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                        plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                        plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                        plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                        plt.xlim(0, 0.5)
                        plt.ylim(0, 0.5)

                        if select_neighboring_cell_type_and_state[0] == False:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(lineage) + '/' + 'p_value legend.png'

                        if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + 'selected neighbors' + '/' + str(lineage) + '/' + 'p_value legend.png'

                        plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                    list_1 = [element[x] for element in state_list_percent_1_final]
                    list_2 = [element[x] for element in state_list_percent_2_final]

                    plot_list = [list_1, list_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage) + str(neighboring_state_list[x]))
                    axs.set_ylabel('Percent ' + str(neighboring_state_list[x]))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + 'selected neighbors' + '/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                    if select_neighboring_cell_type_and_state[0] == False:

                        data_1[neighboring_state_list[x]].append(statistics.mean(list_1))
                        data_2[neighboring_state_list[x]].append(statistics.mean(list_2))

            if select_neighboring_cell_type_and_state[0] == False:

                df_1 = pd.DataFrame(data_1, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.csv'
                df_1.to_csv(export_file_path)

                df_2 = pd.DataFrame(data_2, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.csv'
                df_2.to_csv(export_file_path)

                nan_mask_1 = df_1.isna()
                nan_mask_2 = df_2.isna()

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_1, vmax=100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_1, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_2, vmax=100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_2, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per file/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')
                plt.close()

        if compare_means[0] == True:

            if not compare_means[1]:
                mean_cols = [col for col in neighboring_cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

            else:
                mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + 'selected neighbors/' + str(lineage))

                plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                plt.xlim(0, 0.5)
                plt.ylim(0, 0.5)

                if select_neighboring_cell_type_and_state[0] == False:

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + str(lineage) + '/' + 'p_value legend.png'

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + 'selected neighbors/' + str(lineage) + '/' + 'p_value legend.png'

                plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                for mean_col in mean_cols:

                    mean_list_1_final = []
                    mean_list_2_final = []

                    UniqueFileName = neighboring_cell_type_state_and_neighbor_1.ImageName.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.ImageName == key].reset_index()

                    for file in UniqueFileName:

                        Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

                        mean_list_1 = []
                        for index, cell in Cells_in_file.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:
                                    
                                    mean_list_1.append(cell[mean_col])

                        if len(mean_list_1) == 0:
                            continue
                        else:
                            mean_list_1_final.append(statistics.mean(mean_list_1))

                    UniqueFileName = neighboring_cell_type_state_and_neighbor_2.ImageName.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueFileName}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.ImageName == key].reset_index()

                    for file in UniqueFileName:

                        Cells_in_file = pd.DataFrame(DataFrameDict_Full_cell[file])

                        mean_list_2 = []
                        for index, cell in Cells_in_file.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:
                                    mean_list_2.append(cell[mean_col])

                        if len(mean_list_2) == 0:
                            continue
                        else:
                            mean_list_2_final.append(statistics.mean(mean_list_2))

                    plot_list = [mean_list_1_final, mean_list_2_final]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_ylabel(str(mean_col))
                    axs.set_title(str(mean_col))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + str(lineage) + '/' + str(mean_col) + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per file/' + 'selected neighbors/' + str(lineage) + '/' + str(mean_col) + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

    if split_by == 'image':

        state_list = []

        for row in cell_type_state_and_neighbor_1[df_column]:
            cell_types_and_states = ast.literal_eval(row)

            for type_and_state in cell_types_and_states:

                if type_and_state[0] == select_cell_type_and_state[0][0]:

                    for state in type_and_state[1:]:
                        state_list.append(state)

        for row in cell_type_state_and_neighbor_2[df_column]:
            cell_types_and_states = ast.literal_eval(row)

            for type_and_state in cell_types_and_states:

                if type_and_state[0] == select_cell_type_and_state[1][0]:

                    for state in type_and_state[1:]:
                        state_list.append(state)

        state_list = list(set(state_list))

        UniqueImage = cell_type_state_and_neighbor_1.Individual_image_counter.unique()
        DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
        for key in DataFrameDict_Full_cell.keys():
            DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_1[:][cell_type_state_and_neighbor_1.Individual_image_counter == key].reset_index()

        state_list_percent_1_final = []

        for Image in UniqueImage:

            Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

            state_list_counter_cells_1 = [0 for x in range(0, len(state_list))]

            for row in Cells_in_image[df_column]:

                cell_types_and_states = ast.literal_eval(row)
                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                for type_and_state in cell_types_and_states:

                    if type_and_state[0] == select_cell_type_and_state[0][0]:

                        for state in type_and_state[1:]:
                            state_list_counter_cell[state_list.index(state)] = 1

                state_list_counter_cells_1 = [state_list_counter_cells_1[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

            state_list_percent_1 = [(state_count / Cells_in_image.shape[0]) * 100 for state_count in state_list_counter_cells_1]

            state_list_percent_1_final.append(state_list_percent_1)

        UniqueImage = cell_type_state_and_neighbor_2.Individual_image_counter.unique()
        DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
        for key in DataFrameDict_Full_cell.keys():
            DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_2[:][cell_type_state_and_neighbor_2.Individual_image_counter == key].reset_index()

        state_list_percent_2_final = []

        for Image in UniqueImage:

            Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

            state_list_counter_cells_2 = [0 for x in range(0, len(state_list))]

            for row in Cells_in_image[df_column]:

                cell_types_and_states = ast.literal_eval(row)
                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                for type_and_state in cell_types_and_states:

                    if type_and_state[0] == select_cell_type_and_state[0][0]:

                        for state in type_and_state[1:]:
                            state_list_counter_cell[state_list.index(state)] = 1

                state_list_counter_cells_2 = [state_list_counter_cells_2[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

            state_list_percent_2 = [(state_count / Cells_in_image.shape[0]) * 100 for state_count in state_list_counter_cells_2]

            state_list_percent_2_final.append(state_list_percent_2)

        for x in range(0, len(state_list)):

            if x == 0:

                plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                plt.xlim(0, 0.5)
                plt.ylim(0, 0.5)

                export_file_path = final_output_folder.replace('+', '') + '/' + 'Percentages per image/' + 'p_value legend.png'
                plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

            list_1 = [element[x] for element in state_list_percent_1_final]
            list_2 = [element[x] for element in state_list_percent_2_final]

            plot_list = [list_1, list_2]

            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list, widths=0.4)

            if 0.005 < p < 0.05:
                axs.plot([1.25, 1.75],[max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
            elif 0.0005 < p < 0.005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
            elif 0.00005 < p < 0.0005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
            elif p < 0.00005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
            axs.set_title(str(state_list[x]))
            axs.set_ylabel('Percent ' + str(state_list[x]))
            axs.set_xticks([1, 2],
                           [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                            str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                           rotation=-90)

            export_file_path = final_output_folder.replace('+', '') + '/' + 'Percentages per image/' + str(state_list[x]) + '.png'
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
            plt.close()

        if compare_means[0] == True:

            if not compare_means[1]:

                mean_cols = [col for col in cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

            else:

                mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

            plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
            plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
            plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
            plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
            plt.xlim(0, 0.5)
            plt.ylim(0, 0.5)

            export_file_path = final_output_folder.replace('+', '') + '/' + 'Mean per image/' + 'p_value legend.png'
            plt.savefig(export_file_path, bbox_inches='tight', dpi=600)
                
            for mean_col in mean_cols:

                mean_list_1 = []
                mean_list_2 = []

                UniqueImage = cell_type_state_and_neighbor_1.Individual_image_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_1[:][cell_type_state_and_neighbor_1.Individual_image_counter == key].reset_index()

                for Image in UniqueImage:

                    mean_per_image = pd.DataFrame(DataFrameDict_Full_cell[Image]).loc[:, mean_col].mean()
                    mean_list_1.append(mean_per_image)

                UniqueImage = cell_type_state_and_neighbor_2.Individual_image_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = cell_type_state_and_neighbor_2[:][cell_type_state_and_neighbor_2.Individual_image_counter == key].reset_index()

                for Image in UniqueImage:
                    mean_per_image = pd.DataFrame(DataFrameDict_Full_cell[Image]).loc[:, mean_col].mean()
                    mean_list_2.append(mean_per_image)

                plot_list = [mean_list_1,mean_list_2]

                U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                max_plot = max([max(plot) for plot in plot_list])
                min_plot = min([min(plot) for plot in plot_list])

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                axs.boxplot(plot_list, widths=0.4)

                if 0.005 < p < 0.05:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                elif 0.0005 < p < 0.005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                elif 0.00005 < p < 0.0005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                elif p < 0.00005:
                    axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                    axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                axs.set_ylabel(str(mean_col))
                axs.set_title(str(mean_col))
                axs.set_xticks([1, 2],
                               [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                               rotation=-90)

                export_file_path = final_output_folder.replace('+', '') + '/' + 'Mean per image/' + str(mean_col) + '.png'
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                plt.close()

        if analyse_neighboring_cells == True:

            if select_neighboring_cell_type_and_state[1][0] != select_neighboring_cell_type_and_state[2][0]:
                print('Different neighboring cell lineages are selected. Since they can not be compared directly only the primary selected cell lines are compared.')

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            if select_neighboring_cell_type_and_state[0] == False:

                neighboring_state_list_heat_map = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                neighboring_state_list_heat_map = list(set(neighboring_state_list_heat_map))

                data_1 = {}
                data_2 = {}

                for state in neighboring_state_list_heat_map:
                    data_1[state] = []
                    data_2[state] = []

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + 'selected neighbors/' + str(lineage))

                neighboring_state_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                neighboring_state_list = list(set(neighboring_state_list))

                if select_neighboring_cell_type_and_state[0] == False:

                    non_contained_markers = list(set(neighboring_state_list).symmetric_difference(set(neighboring_state_list_heat_map)))

                    for non_contained_marker in non_contained_markers:

                        data_1[non_contained_marker].append(np.nan)
                        data_2[non_contained_marker].append(np.nan)

                UniqueImage = neighboring_cell_type_state_and_neighbor_1.Individual_image_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.Individual_image_counter == key].reset_index()

                state_list_percent_1_final = []
                cell_lineage_percent_1 = []

                for Image in UniqueImage:

                    Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

                    state_list_counter_cells_1 = [0 for x in range(0, len(neighboring_state_list))]
                    
                    lineage_counter = 0
                    for row in Cells_in_image[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1
                                
                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_1 = [state_list_counter_cells_1[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_1 = [(state_count / Cells_in_image.shape[0]) * 100 for state_count in state_list_counter_cells_1]

                    state_list_percent_1_final.append(state_list_percent_1)
                    cell_lineage_percent_1.append((lineage_counter / Cells_in_image.shape[0]) * 100)

                UniqueImage = neighboring_cell_type_state_and_neighbor_2.Individual_image_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.Individual_image_counter == key].reset_index()

                state_list_percent_2_final = []
                cell_lineage_percent_2 = []

                for Image in UniqueImage:

                    Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

                    state_list_counter_cells_2 = [0 for x in range(0, len(neighboring_state_list))]

                    lineage_counter = 0
                    for row in Cells_in_image[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1

                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_2 = [state_list_counter_cells_2[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_2 = [(state_count / Cells_in_image.shape[0]) * 100 for state_count in state_list_counter_cells_2]

                    state_list_percent_2_final.append(state_list_percent_2)
                    cell_lineage_percent_2.append((lineage_counter / Cells_in_image.shape[0]) * 100)

                if select_neighboring_cell_type_and_state[0] == False:

                    plot_list = [cell_lineage_percent_1, cell_lineage_percent_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage))
                    axs.set_ylabel('Percent ' + str(lineage) + ' cells')
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(lineage) + '/' + str(lineage) + '.png'
                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()
                
                for x in range(0, len(neighboring_state_list)):

                    if x == 0:
                        plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                        plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                        plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                        plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                        plt.xlim(0, 0.5)
                        plt.ylim(0, 0.5)

                        if select_neighboring_cell_type_and_state[0] == False:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(lineage) + '/' + 'p_value legend.png'

                        if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + 'selected neighbors' + '/' + str(lineage) + '/' + 'p_value legend.png'

                        plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                    list_1 = [element[x] for element in state_list_percent_1_final]
                    list_2 = [element[x] for element in state_list_percent_2_final]

                    plot_list = [list_1, list_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage) + str(neighboring_state_list[x]))
                    axs.set_ylabel('Percent ' + str(neighboring_state_list[x]))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + 'selected neighbors' + '/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                    if select_neighboring_cell_type_and_state[0] == False:

                        data_1[neighboring_state_list[x]].append(statistics.mean(list_1))
                        data_2[neighboring_state_list[x]].append(statistics.mean(list_2))

            if select_neighboring_cell_type_and_state[0] == False:

                df_1 = pd.DataFrame(data_1, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.csv'
                df_1.to_csv(export_file_path)

                df_2 = pd.DataFrame(data_2, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.csv'
                df_2.to_csv(export_file_path)

                nan_mask_1 = df_1.isna()
                nan_mask_2 = df_2.isna()

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_1, vmax= 100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_1, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_2, vmax= 100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_2, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percentages per image/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')
                plt.close()

        if compare_means[0] == True:

            if not compare_means[1]:
                mean_cols = [col for col in neighboring_cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

            else:
                mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + 'selected neighbors/' + str(lineage))

                plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                plt.xlim(0, 0.5)
                plt.ylim(0, 0.5)

                if select_neighboring_cell_type_and_state[0] == False:
                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + str(lineage) + '/' + 'p_value legend.png'

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:
                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + 'selected neighbors/' + str(lineage) + '/' + 'p_value legend.png'

                plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                for mean_col in mean_cols:

                    mean_list_1_final = []
                    mean_list_2_final = []

                    UniqueImage = neighboring_cell_type_state_and_neighbor_1.Individual_image_counter.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.Individual_image_counter == key].reset_index()

                    for Image in UniqueImage:

                        Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

                        mean_list_1 = []
                        for index, cell in Cells_in_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:

                                    mean_list_1.append(cell[mean_col])

                        if len(mean_list_1) == 0:
                            continue
                        else:
                            mean_list_1_final.append(statistics.mean(mean_list_1))

                    UniqueImage = neighboring_cell_type_state_and_neighbor_2.Individual_image_counter.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueImage}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.Individual_image_counter == key].reset_index()

                    for Image in UniqueImage:

                        Cells_in_image = pd.DataFrame(DataFrameDict_Full_cell[Image])

                        mean_list_2 = []
                        for index, cell in Cells_in_image.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:

                                    mean_list_2.append(cell[mean_col])

                        if len(mean_list_2) == 0:
                            continue
                        else:
                            mean_list_2_final.append(statistics.mean(mean_list_2))

                    plot_list = [mean_list_1_final, mean_list_2_final]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_ylabel(str(mean_col))
                    axs.set_title(str(mean_col))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + str(lineage) + '/' + str(mean_col) + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per image/' + 'selected neighbors/' + str(lineage) + '/' + str(mean_col) + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

    if split_by == 'cell':

        if not compare_means[1]:
            mean_cols = [col for col in cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

        else:
            mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

        for mean_col in mean_cols:

            plot_list = [cell_type_state_and_neighbor_1[mean_col].to_list(), cell_type_state_and_neighbor_2[mean_col].to_list()]

            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list, widths=0.4)

            if 0.005 < p < 0.05:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
            elif 0.0005 < p < 0.005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
            elif 0.00005 < p < 0.0005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
            elif p < 0.00005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
            axs.set_ylabel(str(mean_col))
            axs.set_title(str(mean_col))
            axs.set_xticks([1, 2],
                           [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                            str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                           rotation=-90)

            export_file_path = final_output_folder.replace('+', '') + '/' + 'Mean per cell/' + str(mean_col) + '.png'
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
            plt.close()

        if analyse_neighboring_cells == True:

            if select_neighboring_cell_type_and_state[1][0] != select_neighboring_cell_type_and_state[2][0]:
                print('Different neighboring cell lineages are selected. Since they can not be compared directly only the primary selected cell lines are compared.')

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            if select_neighboring_cell_type_and_state[0] == False:

                neighboring_state_list_heat_map = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        for state in type_and_state[1:]:
                            neighboring_state_list_heat_map.append(state)

                neighboring_state_list_heat_map = list(set(neighboring_state_list_heat_map))

                data_1 = {}
                data_2 = {}

                for state in neighboring_state_list_heat_map:
                    data_1[state] = []
                    data_2[state] = []

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + 'selected neighbors/' + str(lineage))

                neighboring_state_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == lineage:

                            for state in type_and_state[1:]:
                                neighboring_state_list.append(state)

                neighboring_state_list = list(set(neighboring_state_list))

                if select_neighboring_cell_type_and_state[0] == False:

                    non_contained_markers = list(set(neighboring_state_list).symmetric_difference(set(neighboring_state_list_heat_map)))

                    for non_contained_marker in non_contained_markers:
                        data_1[non_contained_marker].append(np.nan)
                        data_2[non_contained_marker].append(np.nan)

                UniqueCell = neighboring_cell_type_state_and_neighbor_1.Selected_cell_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueCell}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.Selected_cell_counter == key].reset_index()

                state_list_percent_1_final = []
                cell_lineage_percent_1 = []

                for unique_cells in UniqueCell:

                    Cells_in_cells = pd.DataFrame(DataFrameDict_Full_cell[unique_cells])

                    state_list_counter_cells_1 = [0 for x in range(0, len(neighboring_state_list))]

                    lineage_counter = 0
                    for row in Cells_in_cells[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1

                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_1 = [state_list_counter_cells_1[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_1 = [(state_count / Cells_in_cells.shape[0]) * 100 for state_count in state_list_counter_cells_1]

                    state_list_percent_1_final.append(state_list_percent_1)
                    cell_lineage_percent_1.append((lineage_counter / Cells_in_cells.shape[0]) * 100)

                UniqueCell = neighboring_cell_type_state_and_neighbor_2.Selected_cell_counter.unique()
                DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueCell}
                for key in DataFrameDict_Full_cell.keys():
                    DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.Selected_cell_counter == key].reset_index()

                state_list_percent_2_final = []
                cell_lineage_percent_2 = []

                for unique_cells in UniqueCell:

                    Cells_in_cells = pd.DataFrame(DataFrameDict_Full_cell[unique_cells])

                    state_list_counter_cells_2 = [0 for x in range(0, len(neighboring_state_list))]

                    lineage_counter = 0
                    for row in Cells_in_cells[df_column]:

                        cell_types_and_states = ast.literal_eval(row)
                        state_list_counter_cell = [0 for x in range(0, len(neighboring_state_list))]
                        for type_and_state in cell_types_and_states:

                            if type_and_state[0] == lineage:
                                lineage_counter = lineage_counter + 1

                                for state in type_and_state[1:]:
                                    state_list_counter_cell[neighboring_state_list.index(state)] = 1

                        state_list_counter_cells_2 = [state_list_counter_cells_2[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    state_list_percent_2 = [(state_count / Cells_in_cells.shape[0]) * 100 for state_count in state_list_counter_cells_2]

                    state_list_percent_2_final.append(state_list_percent_2)
                    cell_lineage_percent_2.append((lineage_counter / Cells_in_cells.shape[0]) * 100)

                if select_neighboring_cell_type_and_state[0] == False:

                    plot_list = [cell_lineage_percent_1, cell_lineage_percent_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage))
                    axs.set_ylabel('Percent ' + str(lineage) + ' cells')
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(lineage) + '/' + str(lineage) + '.png'
                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()
                
                for x in range(0, len(neighboring_state_list)):

                    if x == 0:
                        plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                        plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                        plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                        plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                        plt.xlim(0, 0.5)
                        plt.ylim(0, 0.5)

                        if select_neighboring_cell_type_and_state[0] == False:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(lineage) + '/' + 'p_value legend.png'

                        if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                            export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + 'selected neighbors' + '/' + str(lineage) + '/' + 'p_value legend.png'

                        plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                    list_1 = [element[x] for element in state_list_percent_1_final]
                    list_2 = [element[x] for element in state_list_percent_2_final]

                    plot_list = [list_1, list_2]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(lineage) + str(neighboring_state_list[x]))
                    axs.set_ylabel('Percent ' + str(neighboring_state_list[x]))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + 'selected neighbors' + '/' + str(lineage) + '/' + neighboring_state_list[x] + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                    if select_neighboring_cell_type_and_state[0] == False:

                        data_1[neighboring_state_list[x]].append(statistics.mean(list_1))
                        data_2[neighboring_state_list[x]].append(statistics.mean(list_2))

            if select_neighboring_cell_type_and_state[0] == False:

                df_1 = pd.DataFrame(data_1, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.csv'
                df_1.to_csv(export_file_path)

                df_2 = pd.DataFrame(data_2, index=lineage_list)
                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.csv'
                df_2.to_csv(export_file_path)

                nan_mask_1 = df_1.isna()
                nan_mask_2 = df_2.isna()

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_1, vmax=100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_1, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(select_cell_type_and_state[0][0]) + ' positive ' + str(select_cell_type_and_state[0][1]) + ' negative ' + str(select_cell_type_and_state[0][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')

                plt.figure(figsize=(22, 10))
                ax = sns.heatmap(df_2, vmax=100, annot=True, cmap="coolwarm", linewidths=0.5, mask=nan_mask_2, linewidth=0.5, norm=LogNorm())
                ax.set(xlabel="", ylabel="")
                ax.xaxis.tick_top()
                plt.title(str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2]))

                export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Percent per cell/' + str(select_cell_type_and_state[1][0]) + ' positive ' + str(select_cell_type_and_state[1][1]) + ' negative ' + str(select_cell_type_and_state[1][2]) + '.png'
                plt.savefig(export_file_path, dpi=300, bbox_inches='tight')
                plt.close()

        if compare_means[0] == True:

            if not compare_means[1]:
                mean_cols = [col for col in neighboring_cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

            else:
                mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]

            if select_neighboring_cell_type_and_state[0] == False:

                lineage_list = []

                for row in neighboring_cell_type_state_and_neighbor_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                for row in neighboring_cell_type_state_and_neighbor_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:
                        lineage_list.append(type_and_state[0])

                lineage_list = list(set(lineage_list))
                lineage_list.remove('CD45+')
                lineage_list.remove('CD45-')

            if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                lineage_list = [select_neighboring_cell_type_and_state[1][0]]

            for lineage in lineage_list:

                if select_neighboring_cell_type_and_state[0] == False:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + str(lineage))

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    if os.path.isdir(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + 'selected neighbors/' + str(lineage)) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + 'selected neighbors/' + str(lineage))

                plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
                plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
                plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
                plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
                plt.xlim(0, 0.5)
                plt.ylim(0, 0.5)

                if select_neighboring_cell_type_and_state[0] == False:

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + str(lineage) + '/' + 'p_value legend.png'

                if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                    export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + 'selected neighbors/' + str(lineage) + '/' + 'p_value legend.png'

                plt.savefig(export_file_path, bbox_inches='tight', dpi=600)

                for mean_col in mean_cols:

                    mean_list_1_final = []
                    mean_list_2_final = []

                    UniqueCell = neighboring_cell_type_state_and_neighbor_1.Selected_cell_counter.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueCell}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_1[:][neighboring_cell_type_state_and_neighbor_1.Selected_cell_counter == key].reset_index()

                    for unique_cell in UniqueCell:

                        Cells_in_cell = pd.DataFrame(DataFrameDict_Full_cell[unique_cell])

                        mean_list_1 = []
                        for index, cell in Cells_in_cell.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:

                                    mean_list_1.append(cell[mean_col])

                        if len(mean_list_1) == 0:
                            continue
                        else:
                            mean_list_1_final.append(statistics.mean(mean_list_1))

                    UniqueCell = neighboring_cell_type_state_and_neighbor_2.Selected_cell_counter.unique()
                    DataFrameDict_Full_cell = {elem: pd.DataFrame() for elem in UniqueCell}
                    for key in DataFrameDict_Full_cell.keys():
                        DataFrameDict_Full_cell[key] = neighboring_cell_type_state_and_neighbor_2[:][neighboring_cell_type_state_and_neighbor_2.Selected_cell_counter == key].reset_index()

                    for unique_cell in UniqueCell:

                        Cells_in_cell = pd.DataFrame(DataFrameDict_Full_cell[unique_cell])

                        mean_list_2 = []
                        for index, cell in Cells_in_cell.iterrows():

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == lineage:

                                    mean_list_2.append(cell[mean_col])

                        if len(mean_list_2) == 0:
                            continue
                        else:
                            mean_list_2_final.append(statistics.mean(mean_list_2))

                    plot_list = [mean_list_1_final, mean_list_2_final]

                    U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.46, max_plot * 1.25, '*', fontsize=14, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.43, max_plot * 1.25, '**', fontsize=14, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.4, max_plot * 1.25, '***', fontsize=14, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.25, '****', fontsize=14, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_ylabel(str(mean_col))
                    axs.set_title(str(mean_col))
                    axs.set_xticks([1, 2],
                                   [str(select_cell_type_and_state[0][0]) + ' positive: ' + str(select_cell_type_and_state[0][1]) + ' negative: ' + str(select_cell_type_and_state[0][2]),
                                    str(select_cell_type_and_state[1][0]) + ' positive: ' + str(select_cell_type_and_state[1][1]) + ' negative: ' + str(select_cell_type_and_state[1][2])],
                                   rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + str(lineage) + '/' + str(mean_col) + '.png'

                    if select_neighboring_cell_type_and_state[0] == True and select_neighboring_cell_type_and_state[1][0] == select_neighboring_cell_type_and_state[2][0]:

                        export_file_path = final_output_folder.replace('+', '') + '/Neighboring cells/' + 'Mean per cell/' + 'selected neighbors/' + str(lineage) + '/' + str(mean_col) + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()