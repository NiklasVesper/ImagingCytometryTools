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
from scipy.stats import wilcoxon

from assign_phenotypes_and_metadata import meta_data_list_1
from assign_phenotypes_and_metadata import meta_data_list_2

'''
Runs the statistics on the different phenotypes and cell states.

Here, one can select between calculating the overall abundance or cell category (immune/tissue cell) of a phenotype.
In addition to that, the overall cellular state can be accessed, with or without specific neighboring cells.
'''

def compare_the_same_phenotypes_and_neighbors_in_different_tissue(directory, filestring, df_column, output_folder, phenotypes, percent_of='all', add_mixed_cells=False, analyse_neighboring_cells=False, percent_of_neigborhood='all', surrounding_phenotypes = [], select_neighboring_cell_type_and_state=[False, ['', [], []]], analyse_states=False, compare_means=[False], split_by='file', statistics='mannwhitneyu', select_area=False, colors=['#D01515','#142EFB']):

    if statistics == 'wilcoxon':

        from assign_phenotypes_and_metadata import sample_parings

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

                    cell_type_state_and_neighbor_meta_data_1 = Cells.iloc[0:0]
                    cell_type_state_and_neighbor_meta_data_2 = Cells.iloc[0:0]

                    cell_type_state_and_neighbor_meta_data_1['Meta_data_category'] = ''
                    cell_type_state_and_neighbor_meta_data_2['Meta_data_category'] = ''

                elif file_counter > 0:
                    continue

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

                if filestring in filename_string and filename_string in meta_data_list_1:

                    print('Processing: ' + filedir)

                    Cells = pd.read_csv(filedir)

                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]

                    for index, cell in Cells.iterrows():

                        cell_dictionary = cell.to_dict()
                        cell_dictionary['Meta_data_category'] = 1
                        cell_type_state_and_neighbor_meta_data_1.loc[len(cell_type_state_and_neighbor_meta_data_1)] = cell_dictionary

                if filestring in filename_string and filename_string in meta_data_list_2:

                    print('Processing: ' + filedir)

                    Cells = pd.read_csv(filedir)

                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]

                    for index, cell in Cells.iterrows():

                        cell_dictionary = cell.to_dict()
                        cell_dictionary['Meta_data_category'] = 2
                        cell_type_state_and_neighbor_meta_data_2.loc[len(cell_type_state_and_neighbor_meta_data_2)] = cell_dictionary

    analysis_values_phenotyping = ['Percent_of_all_per_image',
                                   'Percent_of_immune_per_image',
                                   'Percent_of_tissue_per_image',
                                   'FileName',
                                   'FileName_Image',
                                   'MetaData']

    lineage_data = [lineage[0] for lineage in phenotypes]
    phenotype_data = [lineage[0] + ' pos ' + str(lineage[1]) + ' neg ' + str(lineage[2]) for lineage in phenotypes]

    for x in range(0, len(phenotypes)):
        analysis_values_phenotyping.append(lineage_data[x])
        analysis_values_phenotyping.append(phenotype_data[x])

    analysis_df_phenotypes = pd.DataFrame(columns=analysis_values_phenotyping)

    if analyse_neighboring_cells == True:

        analysis_values_neigborhood = ['FileName',
                                       'FileName_Image',
                                       'MetaData']

        phenotype_data = [lineage[0] + ' pos ' + str(lineage[1]) + ' neg ' + str(lineage[2]) for lineage in phenotypes]
        neighboring_phenotype_data = [lineage[0] + ' pos ' + str(lineage[1]) + ' neg ' + str(lineage[2]) for lineage in surrounding_phenotypes]

        for x in range(0, len(phenotypes)):
            analysis_values_neigborhood.append('surrounding_cells_per_image_counter ' + str(phenotype_data[x]))
            analysis_values_neigborhood.append('surrounding_immune_cells_per_image_counter ' + str(phenotype_data[x]))
            analysis_values_neigborhood.append('surrounding_tissue_cells_per_image_counter ' + str(phenotype_data[x]))
            for k in range(0, len(surrounding_phenotypes)):
                analysis_values_neigborhood.append(str(phenotype_data[x]) + ' neighbors: ' + str(neighboring_phenotype_data[k]))

        analysis_df_neigborhood = pd.DataFrame(columns=analysis_values_neigborhood)
        
    UniqueFileName = cell_type_state_and_neighbor_meta_data_1.ImageName.unique()
    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

    for key in DataFrameDict_file.keys():
        DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_1[:][cell_type_state_and_neighbor_meta_data_1.ImageName == key].reset_index()

    for file in UniqueFileName:

        Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

        UniqueImage = Cells_in_file.ImageNumber.unique()
        DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

        for key in DataFrameDict_image.keys():
            DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index()

        for image in UniqueImage:

            analysis_row = []
            analysis_row_neighborhood = []
            phenotype_counter = -1
            for phenotype in phenotypes:
                phenotype_counter = phenotype_counter + 1

                if phenotype_counter == 0:
                    Percent_of_all_per_image_counter = 0
                    Percent_of_immune_per_image_counter = 0
                    Percent_of_tissue_per_image_counter = 0

                percent_of_lineage_per_image_counter = 0
                phenotype_per_image_counter_counter = 0

                Cells_in_image = pd.DataFrame(DataFrameDict_image[image])

                surrounding_immune_per_image_counter = 0
                surrounding_tissue_per_image_counter = 0
                surrounding_phenotype_counter_list = [0 for x in range(0, len(surrounding_phenotypes))]

                for index, cell in Cells_in_image.iterrows():

                    if select_area == True:

                        if cell['Selected Areas'] == 'in the area':

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + (len(cell_types_and_states) - 1)

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + (len(cell_types_and_states) - 1)

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                        if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            neighborhood_counter = 0
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                        neighborhood_counter = neighborhood_counter + 1

                                            if neighborhood_counter == 0:
                                                phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                            if neighborhood_counter > 0:

                                                for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                    neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                    if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                        surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                    if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                        surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                    surrounding_phenotype_counter = -1
                                                    for surrounding_phenotype in surrounding_phenotypes:

                                                        surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                        for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                            if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                                surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + 1

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + 1

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    if len(neighboring_cell_types_and_states) == 2:

                                                        if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            neighborhood_counter = 0
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        neighborhood_counter = neighborhood_counter + 1

                                            if neighborhood_counter == 0:
                                                phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                            if neighborhood_counter > 0:

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                    neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                    surrounding_phenotype_counter = -1
                                                    for surrounding_phenotype in surrounding_phenotypes:

                                                        surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                        if len(neighboring_cell_types_and_states) == 2:

                                                            if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                                surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                        elif cell['Selected Areas'] == 'outside the area':

                            continue

                    if select_area == False:

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + (len(cell_types_and_states) - 1)

                            if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + (len(cell_types_and_states) - 1)

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == phenotype[0]:

                                    percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                    phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                    if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            surrounding_phenotype_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:

                                                surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:

                                                        surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter == 0:

                                            phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                        if neighborhood_counter > 0:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                        if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                        if add_mixed_cells == False:

                            if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + 1

                            if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + 1

                            if len(cell_types_and_states) == 2:

                                if cell_types_and_states[0][0] == phenotype[0]:

                                    percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                    phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                    if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                        if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                            surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                        if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                            surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            surrounding_phenotype_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:

                                                surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:

                                                        surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter == 0:

                                            phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                        if neighborhood_counter > 0:

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    if len(neighboring_cell_types_and_states) == 2:

                                                        if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                if phenotype_counter == 0:

                    Percent_of_all_per_image_counter = Percent_of_immune_per_image_counter + Percent_of_tissue_per_image_counter

                    analysis_row.append(Percent_of_all_per_image_counter)
                    analysis_row.append(Percent_of_immune_per_image_counter)
                    analysis_row.append(Percent_of_tissue_per_image_counter)
                    analysis_row.append(file)
                    analysis_row.append(file + '_' + str(image))
                    analysis_row.append(1)

                    if analyse_neighboring_cells == True:
                        analysis_row_neighborhood.append(file)
                        analysis_row_neighborhood.append(file + '_' + str(image))
                        analysis_row_neighborhood.append(1)

                analysis_row.append(percent_of_lineage_per_image_counter)
                analysis_row.append(phenotype_per_image_counter_counter)

                if analyse_neighboring_cells == True:
                    surrounding_cells_per_image_counter = surrounding_immune_per_image_counter + surrounding_tissue_per_image_counter

                    analysis_row_neighborhood.append(surrounding_cells_per_image_counter)
                    analysis_row_neighborhood.append(surrounding_immune_per_image_counter)
                    analysis_row_neighborhood.append(surrounding_tissue_per_image_counter)

                    for count in surrounding_phenotype_counter_list:
                        analysis_row_neighborhood.append(count)

            analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row

            if analyse_neighboring_cells == True:
                analysis_df_neigborhood.loc[len(analysis_df_neigborhood)] = analysis_row_neighborhood

    UniqueFileName = cell_type_state_and_neighbor_meta_data_2.ImageName.unique()
    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

    for key in DataFrameDict_file.keys():
        DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_2[:][cell_type_state_and_neighbor_meta_data_2.ImageName == key].reset_index()

    for file in UniqueFileName:

        Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

        UniqueImage = Cells_in_file.ImageNumber.unique()
        DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

        for key in DataFrameDict_image.keys():
            DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index()

        for image in UniqueImage:

            analysis_row = []
            analysis_row_neighborhood = []
            phenotype_counter = -1
            for phenotype in phenotypes:
                phenotype_counter = phenotype_counter + 1

                if phenotype_counter == 0:
                    Percent_of_all_per_image_counter = 0
                    Percent_of_immune_per_image_counter = 0
                    Percent_of_tissue_per_image_counter = 0

                percent_of_lineage_per_image_counter = 0
                phenotype_per_image_counter_counter = 0

                Cells_in_image = pd.DataFrame(DataFrameDict_image[image])

                surrounding_immune_per_image_counter = 0
                surrounding_tissue_per_image_counter = 0
                surrounding_phenotype_counter_list = [0 for x in range(0, len(surrounding_phenotypes))]

                for index, cell in Cells_in_image.iterrows():

                    if select_area == True:

                        if cell['Selected Areas'] == 'in the area':

                            cell_types_and_states = ast.literal_eval(cell[df_column])

                            if add_mixed_cells == True:

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + (len(cell_types_and_states) - 1)

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + (len(cell_types_and_states) - 1)

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                        if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            neighborhood_counter = 0
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                        neighborhood_counter = neighborhood_counter + 1

                                            if neighborhood_counter == 0:
                                                phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                            if neighborhood_counter > 0:

                                                for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                    neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                    if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                        surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                    if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                        surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                    surrounding_phenotype_counter = -1
                                                    for surrounding_phenotype in surrounding_phenotypes:

                                                        surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                        for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                            if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                                surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + 1

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + 1

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                        if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    if len(neighboring_cell_types_and_states) == 2:

                                                        if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                        if select_neighboring_cell_type_and_state[0] == True:

                                            neighborhood_counter = 0
                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        neighborhood_counter = neighborhood_counter + 1

                                            if neighborhood_counter == 0:
                                                phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                            if neighborhood_counter > 0:

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                    neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                    surrounding_phenotype_counter = -1
                                                    for surrounding_phenotype in surrounding_phenotypes:

                                                        surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                        if len(neighboring_cell_types_and_states) == 2:

                                                            if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                                surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                        elif cell['Selected Areas'] == 'outside the area':

                            continue
                            
                    if select_area == False:

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + (len(cell_types_and_states) - 1)

                            if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + (len(cell_types_and_states) - 1)

                            for type_and_state in cell_types_and_states:

                                if type_and_state[0] == phenotype[0]:
                                    percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                    phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                    if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            surrounding_phenotype_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:

                                                surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                    if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                        surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter == 0:
                                            phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                        if neighborhood_counter > 0:

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                    surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                if len(neighboring_cell_types_and_states) >= 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                    surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                        if neighboring_type_and_state[0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_type_and_state[1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                        if add_mixed_cells == False:

                            if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                Percent_of_immune_per_image_counter = Percent_of_immune_per_image_counter + 1

                            if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                Percent_of_tissue_per_image_counter = Percent_of_tissue_per_image_counter + 1

                            if len(cell_types_and_states) == 2:

                                if cell_types_and_states[0][0] == phenotype[0]:
                                    percent_of_lineage_per_image_counter = percent_of_lineage_per_image_counter + 1

                                if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                    phenotype_per_image_counter_counter = phenotype_per_image_counter_counter + 1

                                    if select_neighboring_cell_type_and_state[0] == False and analyse_neighboring_cells == True:

                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            surrounding_phenotype_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:

                                                surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                if len(neighboring_cell_types_and_states) == 2:

                                                    if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                        surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                                    if select_neighboring_cell_type_and_state[0] == True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter == 0:
                                            phenotype_per_image_counter_counter = phenotype_per_image_counter_counter - 1

                                        if neighborhood_counter > 0:

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45+':
                                                surrounding_immune_per_image_counter = surrounding_immune_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            if len(neighboring_cell_types_and_states) == 2 and neighboring_cell_types_and_states[-1][0] == 'CD45-':
                                                surrounding_tissue_per_image_counter = surrounding_tissue_per_image_counter + (len(neighboring_cell_types_and_states) - 1)

                                            for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                                neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                                surrounding_phenotype_counter = -1
                                                for surrounding_phenotype in surrounding_phenotypes:

                                                    surrounding_phenotype_counter = surrounding_phenotype_counter + 1

                                                    if len(neighboring_cell_types_and_states) == 2:

                                                        if neighboring_cell_types_and_states[0][0] == surrounding_phenotype[0] and set(surrounding_phenotype[1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(surrounding_phenotype[2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                            surrounding_phenotype_counter_list[surrounding_phenotype_counter] = surrounding_phenotype_counter_list[surrounding_phenotype_counter] + 1

                if phenotype_counter == 0:
                    Percent_of_all_per_image_counter = Percent_of_immune_per_image_counter + Percent_of_tissue_per_image_counter

                    analysis_row.append(Percent_of_all_per_image_counter)
                    analysis_row.append(Percent_of_immune_per_image_counter)
                    analysis_row.append(Percent_of_tissue_per_image_counter)
                    analysis_row.append(file)
                    analysis_row.append(file + '_' + str(image))
                    analysis_row.append(2)

                    if analyse_neighboring_cells == True:
                        analysis_row_neighborhood.append(file)
                        analysis_row_neighborhood.append(file + '_' + str(image))
                        analysis_row_neighborhood.append(2)

                analysis_row.append(percent_of_lineage_per_image_counter)
                analysis_row.append(phenotype_per_image_counter_counter)

                if analyse_neighboring_cells == True:

                    surrounding_cells_per_image_counter = surrounding_immune_per_image_counter + surrounding_tissue_per_image_counter

                    analysis_row_neighborhood.append(surrounding_cells_per_image_counter)
                    analysis_row_neighborhood.append(surrounding_immune_per_image_counter)
                    analysis_row_neighborhood.append(surrounding_tissue_per_image_counter)

                    for count in surrounding_phenotype_counter_list:
                        analysis_row_neighborhood.append(count)

            analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row

            if analyse_neighboring_cells == True:
                analysis_df_neigborhood.loc[len(analysis_df_neigborhood)] = analysis_row_neighborhood

    analysis_df_phenotypes = analysis_df_phenotypes.loc[:, ~analysis_df_phenotypes.columns.duplicated()].copy()

    remove_row_list = list(analysis_df_phenotypes.loc[~(analysis_df_phenotypes['Percent_of_all_per_image'] > 0)]['FileName_Image'])

    if select_area == True:
        
        analysis_df_phenotypes = analysis_df_phenotypes[~analysis_df_phenotypes['FileName_Image'].isin(remove_row_list)]

    columns_final_values = ['FileName', 'MetaData']
    df_final_values = pd.DataFrame(columns=columns_final_values)

    for phenotype in phenotypes:

        df_final_values_phenotype_column = []
        FileName = []
        MetaData = []

        plot_list_1 = []
        plot_list_2 = []

        if statistics == 'mannwhitneyu':

            UniqueMeta_data = analysis_df_phenotypes.MetaData.unique()
            DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_phenotypes.MetaData.unique()}

            for key in DataFrameDict_Metadata.keys():
                DataFrameDict_Metadata[key] = analysis_df_phenotypes[:][analysis_df_phenotypes.MetaData == key].reset_index()

            for metadata in UniqueMeta_data:

                Cells_in_metadata = DataFrameDict_Metadata[metadata]

                if split_by == 'file':

                    UniqueFile = Cells_in_metadata.FileName.unique()
                    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                    for key in DataFrameDict_file.keys():
                        DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                    for file in UniqueFile:

                        Cells_in_file = DataFrameDict_file[file]

                        if percent_of == 'all':
                            sum_of_cells = sum(Cells_in_file['Percent_of_all_per_image'].to_list())

                        if percent_of == 'immune':
                            sum_of_cells = sum(Cells_in_file['Percent_of_immune_per_image'].to_list())

                        if percent_of == 'tissue':
                            sum_of_cells = sum(Cells_in_file['Percent_of_tissue_per_image'].to_list())

                        if percent_of == 'lineage':
                            sum_of_cells = sum(Cells_in_file[phenotype[0]].to_list())

                        phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                        sum_of_selected_cells = sum(Cells_in_file[phenotype_df].to_list())

                        if metadata == 1:

                            MetaData.append(1)

                            if int(sum_of_selected_cells) > 0:
                                plot_list_1.append((sum_of_selected_cells/sum_of_cells)*100)
                                df_final_values_phenotype_column.append((sum_of_selected_cells/sum_of_cells)*100)
                            if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                plot_list_1.append(0)
                                df_final_values_phenotype_column.append(0)

                        if metadata == 2:

                            MetaData.append(2)

                            if int(sum_of_selected_cells) > 0:
                                plot_list_2.append((sum_of_selected_cells/sum_of_cells) * 100)
                                df_final_values_phenotype_column.append((sum_of_selected_cells/sum_of_cells) * 100)
                            if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                plot_list_2.append(0)
                                df_final_values_phenotype_column.append(0)

                        FileName.append(file)

                if split_by == 'image':

                    UniqueImage = Cells_in_metadata.FileName_Image.unique()
                    DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                    for key in DataFrameDict_image.keys():
                        DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                    for image in UniqueImage:

                        Cells_in_image = DataFrameDict_image[image]

                        if percent_of == 'all':
                            cells = Cells_in_image['Percent_of_all_per_image']

                        if percent_of == 'immune':
                            cells = Cells_in_image['Percent_of_immune_per_image']

                        if percent_of == 'tissue':
                            cells = Cells_in_image['Percent_of_tissue_per_image']

                        if percent_of == 'lineage':
                            cells = Cells_in_image[phenotype[0]]

                        phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                        selected_cells = Cells_in_image[phenotype_df]

                        if metadata == 1:

                            if int(selected_cells) > 0:
                                plot_list_1.append(float((selected_cells / cells) * 100))
                            if int(selected_cells) == 0 or int(cells) == 0:
                                plot_list_1.append(0)

                        if metadata == 2:

                            if int(selected_cells) > 0:
                                plot_list_2.append(float((selected_cells / cells) * 100))
                            if int(selected_cells) == 0 or int(cells) == 0:
                                plot_list_2.append(0)

        if statistics == 'wilcoxon':

            UniqueMeta_data = analysis_df_phenotypes.MetaData.unique()
            DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_phenotypes.MetaData.unique()}

            for key in DataFrameDict_Metadata.keys():
                DataFrameDict_Metadata[key] = analysis_df_phenotypes[:][analysis_df_phenotypes.MetaData == key].reset_index()

            for metadata in UniqueMeta_data:

                Cells_in_metadata = DataFrameDict_Metadata[metadata]

                if split_by == 'file':

                    if metadata == 1:

                        UniqueFile = Cells_in_metadata.FileName.unique()
                        DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                        for key in DataFrameDict_file.keys():
                            DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                        for file in UniqueFile:

                            FileName.append(file)
                            MetaData.append(1)

                            Cells_in_file = DataFrameDict_file[file]

                            if percent_of == 'all':
                                sum_of_cells = sum(Cells_in_file['Percent_of_all_per_image'].to_list())

                            if percent_of == 'immune':
                                sum_of_cells = sum(Cells_in_file['Percent_of_immune_per_image'].to_list())

                            if percent_of == 'tissue':
                                sum_of_cells = sum(Cells_in_file['Percent_of_tissue_per_image'].to_list())

                            if percent_of == 'lineage':
                                sum_of_cells = sum(Cells_in_file[phenotype[0]].to_list())

                            phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                            sum_of_selected_cells = sum(Cells_in_file[phenotype_df].to_list())

                            if int(sum_of_selected_cells) > 0:
                                plot_list_1.append((sum_of_selected_cells / sum_of_cells) * 100)
                                df_final_values_phenotype_column.append((sum_of_selected_cells / sum_of_cells) * 100)

                            if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                plot_list_1.append(0)
                                df_final_values_phenotype_column.append(0)

                            paired_file = sample_parings[file]
                            paired_file_values = analysis_df_phenotypes.loc[analysis_df_phenotypes['FileName'] == paired_file]

                            FileName.append(paired_file)
                            MetaData.append(2)

                            if percent_of == 'all':
                                sum_of_paired_cells = sum(paired_file_values['Percent_of_all_per_image'].to_list())

                            if percent_of == 'immune':
                                sum_of_paired_cells = sum(paired_file_values['Percent_of_immune_per_image'].to_list())

                            if percent_of == 'tissue':
                                sum_of_paired_cells = sum(paired_file_values['Percent_of_tissue_per_image'].to_list())

                            if percent_of == 'lineage':
                                sum_of_paired_cells = sum(paired_file_values[phenotype[0]].to_list())

                            phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                            sum_of_selected_paired_cells = sum(paired_file_values[phenotype_df].to_list())

                            if int(sum_of_selected_paired_cells) > 0:
                                plot_list_2.append((sum_of_selected_paired_cells / sum_of_paired_cells) * 100)
                                df_final_values_phenotype_column.append((sum_of_selected_paired_cells / sum_of_paired_cells) * 100)

                            if int(sum_of_selected_paired_cells) == 0 or int(sum_of_paired_cells) == 0:
                                plot_list_2.append(0)
                                df_final_values_phenotype_column.append(0)

                    else:
                        continue

                if split_by == 'image':

                    UniqueImage = Cells_in_metadata.FileName_Image.unique()
                    DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                    for key in DataFrameDict_image.keys():
                        DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                    for image in UniqueImage:

                        Cells_in_image = DataFrameDict_image[image]

                        if percent_of == 'all':
                            cells = Cells_in_image['Percent_of_all_per_image']

                        if percent_of == 'immune':
                            cells = Cells_in_image['Percent_of_immune_per_image']

                        if percent_of == 'tissue':
                            cells = Cells_in_image['Percent_of_tissue_per_image']

                        if percent_of == 'lineage':
                            cells = Cells_in_image[phenotype[0]]

                        phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                        selected_cells = Cells_in_image[phenotype_df]

                        if metadata == 1:

                            if int(selected_cells) > 0:
                                plot_list_1.append(float((selected_cells / cells) * 100))
                            if int(selected_cells) == 0 or int(cells) == 0:
                                plot_list_1.append(0)

                        if metadata == 2:

                            if int(selected_cells) > 0:
                                plot_list_2.append(float((selected_cells / cells) * 100))
                            if int(selected_cells) == 0 or int(cells) == 0:
                                plot_list_2.append(0)

        df_final_values[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])] = df_final_values_phenotype_column
        df_final_values['FileName'] = FileName
        df_final_values['MetaData'] = MetaData

        plot_list = [plot_list_1, plot_list_2]

        try:

            if statistics == 'mannwhitneyu':
                U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')

            if statistics == 'wilcoxon':
                result = wilcoxon(plot_list[0], plot_list[1])
                p = result.pvalue

            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

            if statistics == 'wilcoxon' and split_by == 'file':

                for x in range(0, len(plot_list[0])):
                    axs.plot([1, 2], [plot_list[0][x], plot_list[1][x]], color='gray', alpha=0.25)

            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c=colors[0], s=45)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c=colors[1], s=45)
            axs.boxplot(plot_list, widths=0.4)

            if 0.005 < p < 0.05:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.45, max_plot * 1.2, '*', fontsize=18, color='black')
            elif 0.0005 < p < 0.005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.41, max_plot * 1.2, '**', fontsize=18, color='black')
            elif 0.00005 < p < 0.0005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.37, max_plot * 1.2, '***', fontsize=18, color='black')
            elif p < 0.00005:
                axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                axs.text(1.34, max_plot * 1.2, '****', fontsize=18, color='black')

            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
            axs.set_xticks([1, 2], ['Metadata 1','Metadata 2'], rotation=-90)

            if select_neighboring_cell_type_and_state[0] == False:

                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per ' + split_by + '/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per ' + split_by)
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per ' + split_by + '/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per ' + split_by)
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '.png'

            if select_neighboring_cell_type_and_state[0] == True:

                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per ' + split_by + '/' + 'neighboring ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '\n neigboring: ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per ' + split_by)
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neighboring ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per ' + split_by + '/'  + 'neighboring ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '\n neigboring: ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per ' + split_by)
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neighboring ' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '.png'

            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
            plt.close()

        except ValueError:
            print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not exist')
            continue

    plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
    plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
    plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
    plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
    plt.xlim(0, 0.5)
    plt.ylim(0, 0.5)

    plt.savefig(output_folder + '/' + 'Phenotyping/' + 'p_value_legend.png', bbox_inches='tight', dpi=600)
    plt.close()

    if select_neighboring_cell_type_and_state[0] == False:

        analysis_df_phenotypes.to_csv(output_folder + '/' + 'Phenotyping/' + 'Phenotypes.csv')
        df_final_values.to_csv(final_output_folder + 'Phenotype_plots.csv')

    if select_neighboring_cell_type_and_state[0] == True:

        analysis_df_phenotypes.to_csv(output_folder + '/' + 'Phenotyping/' + 'Phenotypes neighboring' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '.csv')
        df_final_values.to_csv(final_output_folder + 'Phenotypes neighboring' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '_plots.csv')

#-----------------------------------------------------------------------------------------------------------------------------------

    if analyse_neighboring_cells == True:

        if select_area == True:

            analysis_df_neigborhood = analysis_df_neigborhood[~analysis_df_neigborhood['FileName_Image'].isin(remove_row_list)]

        data = analysis_df_neigborhood.iloc[:, :3]

        phenotype_counter = -1
        for phenotype in phenotypes:

            columns_final_values = ['FileName', 'MetaData']
            df_final_values_neigborhood = pd.DataFrame(columns=columns_final_values)

            phenotype_counter = phenotype_counter + 1

            analysis_df_neigborhood_sub = analysis_df_neigborhood.iloc[:,((3*(phenotype_counter+1))+(phenotype_counter*len(surrounding_phenotypes))):((3*(phenotype_counter+2))+((phenotype_counter+1)*len(surrounding_phenotypes)))]

            analysis_df_neigborhood_sub= pd.concat([data, analysis_df_neigborhood_sub], axis=1)

            for surrounding_phenotype in surrounding_phenotypes:

                df_final_values_phenotype_column = []
                FileName = []
                MetaData = []

                plot_list_1 = []
                plot_list_2 = []

                if statistics == 'mannwhitneyu':

                    UniqueMeta_data = analysis_df_neigborhood_sub.MetaData.unique()
                    DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_neigborhood_sub.MetaData.unique()}

                    for key in DataFrameDict_Metadata.keys():
                        DataFrameDict_Metadata[key] = analysis_df_neigborhood_sub[:][analysis_df_neigborhood_sub.MetaData == key].reset_index()

                    for metadata in UniqueMeta_data:

                        Cells_in_metadata = DataFrameDict_Metadata[metadata]

                        if split_by == 'file':

                            UniqueFile = Cells_in_metadata.FileName.unique()
                            DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                            for key in DataFrameDict_file.keys():
                                DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                            for file in UniqueFile:

                                Cells_in_file = DataFrameDict_file[file]

                                if percent_of_neigborhood == 'all':
                                    sum_of_cells = sum(Cells_in_file['surrounding_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                if percent_of_neigborhood == 'immune':
                                    sum_of_cells = sum(Cells_in_file['surrounding_immune_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                if percent_of_neigborhood == 'tissue':
                                    sum_of_cells = sum(Cells_in_file['surrounding_tissue_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                                surrounding_phenotype_df = str(surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))

                                sum_of_selected_cells = sum(Cells_in_file[phenotype_df + ' neighbors: ' + surrounding_phenotype_df].to_list())

                                if metadata == 1:

                                    MetaData.append(1)

                                    if int(sum_of_selected_cells) > 0:
                                        plot_list_1.append((sum_of_selected_cells / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_selected_cells / sum_of_cells) * 100)
                                    if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                        plot_list_1.append(0)
                                        df_final_values_phenotype_column.append(0)

                                if metadata == 2:

                                    MetaData.append(2)

                                    if int(sum_of_selected_cells) > 0:
                                        plot_list_2.append((sum_of_selected_cells / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_selected_cells / sum_of_cells) * 100)
                                    if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                        plot_list_2.append(0)
                                        df_final_values_phenotype_column.append(0)

                                FileName.append(file)

                        if split_by == 'image':

                            UniqueImage = Cells_in_metadata.FileName_Image.unique()
                            DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                            for key in DataFrameDict_image.keys():
                                DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                            for image in UniqueImage:

                                Cells_in_image = DataFrameDict_image[image]

                                if percent_of_neigborhood == 'all':
                                    cells = Cells_in_image['surrounding_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                if percent_of_neigborhood == 'immune':
                                    cells = Cells_in_image['surrounding_immune_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                if percent_of_neigborhood == 'tissue':
                                    cells = Cells_in_image['surrounding_tissue_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                                surrounding_phenotype_df = str(surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))

                                selected_cells = Cells_in_image[phenotype_df + ' neighbors: ' + surrounding_phenotype_df]

                                if metadata == 1:

                                    if int(selected_cells) > 0:
                                        plot_list_1.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_1.append(0)

                                if metadata == 2:

                                    if int(selected_cells) > 0:
                                        plot_list_2.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_2.append(0)

                if statistics == 'wilcoxon':

                    UniqueMeta_data = analysis_df_neigborhood_sub.MetaData.unique()
                    DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_neigborhood_sub.MetaData.unique()}

                    for key in DataFrameDict_Metadata.keys():
                        DataFrameDict_Metadata[key] = analysis_df_neigborhood_sub[:][analysis_df_neigborhood_sub.MetaData == key].reset_index()

                    for metadata in UniqueMeta_data:

                        Cells_in_metadata = DataFrameDict_Metadata[metadata]

                        if split_by == 'file':

                            if metadata == 1:

                                UniqueFile = Cells_in_metadata.FileName.unique()
                                DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                                for key in DataFrameDict_file.keys():
                                    DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                                for file in UniqueFile:

                                    FileName.append(file)
                                    MetaData.append(1)

                                    Cells_in_file = DataFrameDict_file[file]

                                    if percent_of_neigborhood == 'all':
                                        sum_of_cells = sum(Cells_in_file['surrounding_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    if percent_of_neigborhood == 'immune':
                                        sum_of_cells = sum(Cells_in_file['surrounding_immune_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    if percent_of_neigborhood == 'tissue':
                                        sum_of_cells = sum(Cells_in_file['surrounding_tissue_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                                    surrounding_phenotype_df = str(surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))

                                    sum_of_selected_cells = sum(Cells_in_file[phenotype_df + ' neighbors: ' + surrounding_phenotype_df].to_list())

                                    if int(sum_of_selected_cells) > 0:
                                        plot_list_1.append((sum_of_selected_cells / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_selected_cells / sum_of_cells) * 100)
                                    if int(sum_of_selected_cells) == 0 or int(sum_of_cells) == 0:
                                        plot_list_1.append(0)
                                        df_final_values_phenotype_column.append(0)

                                    paired_file = sample_parings[file]
                                    paired_file_values = analysis_df_neigborhood_sub.loc[analysis_df_neigborhood_sub['FileName'] == paired_file]

                                    FileName.append(paired_file)
                                    MetaData.append(2)

                                    if percent_of_neigborhood == 'all':
                                        sum_of_paired_cells = sum(paired_file_values['surrounding_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    if percent_of_neigborhood == 'immune':
                                        sum_of_paired_cells = sum(paired_file_values['surrounding_immune_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    if percent_of_neigborhood == 'tissue':
                                        sum_of_paired_cells = sum(paired_file_values['surrounding_tissue_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])].to_list())

                                    phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                                    surrounding_phenotype_df = str(surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))

                                    sum_of_selected_paired_cells = sum(paired_file_values[phenotype_df + ' neighbors: ' + surrounding_phenotype_df].to_list())

                                    if int(sum_of_selected_paired_cells) > 0:
                                        plot_list_2.append((sum_of_selected_paired_cells / sum_of_paired_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_selected_paired_cells / sum_of_paired_cells) * 100)
                                    if int(sum_of_selected_paired_cells) == 0 or int(sum_of_paired_cells) == 0:
                                        plot_list_2.append(0)
                                        df_final_values_phenotype_column.append(0)

                            else:
                                continue

                        if split_by == 'image':

                            UniqueImage = Cells_in_metadata.FileName_Image.unique()
                            DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                            for key in DataFrameDict_image.keys():
                                DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                            for image in UniqueImage:

                                Cells_in_image = DataFrameDict_image[image]

                                if percent_of_neigborhood == 'all':
                                    cells = Cells_in_image['surrounding_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                if percent_of_neigborhood == 'immune':
                                    cells = Cells_in_image['surrounding_immune_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                if percent_of_neigborhood == 'tissue':
                                    cells = Cells_in_image['surrounding_tissue_cells_per_image_counter ' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])]

                                phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                                surrounding_phenotype_df = str(surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))

                                selected_cells = Cells_in_image[phenotype_df + ' neighbors: ' + surrounding_phenotype_df]

                                if metadata == 1:

                                    if int(selected_cells) > 0:
                                        plot_list_1.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_1.append(0)

                                if metadata == 2:

                                    if int(selected_cells) > 0:
                                        plot_list_2.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_2.append(0)

                df_final_values_neigborhood[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neigboring: ' + surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2])] = df_final_values_phenotype_column
                df_final_values_neigborhood['FileName'] = FileName
                df_final_values_neigborhood['MetaData'] = MetaData
            
                plot_list = [plot_list_1, plot_list_2]

                try:

                    if statistics == 'mannwhitneyu':
                        U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')

                    if statistics == 'wilcoxon':
                        res = wilcoxon(plot_list[0], plot_list[1])
                        p = res.pvalue

                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    if statistics == 'wilcoxon' and split_by == 'file':

                        for x in range(0, len(plot_list[0])):
                            axs.plot([1, 2], [plot_list[0][x], plot_list[1][x]], color='gray', alpha=0.25)

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c=colors[0], s=45) #D01515
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c=colors[1], s=45)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.45, max_plot * 1.2, '*', fontsize=18, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.41, max_plot * 1.2, '**', fontsize=18, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.2, '***', fontsize=18, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.34, max_plot * 1.2, '****', fontsize=18, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '\n neigboring: ' + surrounding_phenotype[0] + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' neighboring cells')
                    axs.set_xticks([1, 2], ['Metadata 1','Metadata 2'], rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        final_output_folder = output_folder + '/' + 'Neighborhood analysis/' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' Percent of ' + percent_of + ' cells per ' + split_by + '/'

                        if os.path.isdir(final_output_folder) == True:
                            pass
                        else:
                            os.makedirs(final_output_folder)

                        export_file_path = final_output_folder + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neigboring ' + str(surrounding_phenotype[0]) + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]) + '.png'

                    if select_neighboring_cell_type_and_state[0] == True:

                        final_output_folder = output_folder + '/' + 'Selected neighborhood analysis/' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '/' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' Percent of ' + percent_of + ' cells per ' + split_by + '/'

                        if os.path.isdir(final_output_folder) == True:
                            pass
                        else:
                            os.makedirs(final_output_folder)

                        export_file_path = final_output_folder + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neigboring ' + str(surrounding_phenotype[0]) + ' pos ' + str(surrounding_phenotype[1]) + ' neg ' + str(surrounding_phenotype[2]) + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                except ValueError:
                    print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not have: ' + str(surrounding_phenotype) + ' as neighbor')
                    continue

            df_final_values_neigborhood.to_csv(final_output_folder + 'Neighborhood_plots.csv')

        plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
        plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
        plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
        plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
        plt.xlim(0, 0.5)
        plt.ylim(0, 0.5)

        if select_neighboring_cell_type_and_state[0] == False:

            plt.savefig(output_folder + '/' + 'Neighborhood analysis/' + 'p_value_legend.png', bbox_inches='tight', dpi=600)
            plt.close()

            analysis_df_neigborhood.to_csv(output_folder + '/' + 'Neighborhood analysis/' + 'Neighborhood.csv')

        if select_neighboring_cell_type_and_state[0] == True:

            plt.savefig(output_folder + '/' + 'Selected neighborhood analysis/' + 'p_value_legend.png', bbox_inches='tight', dpi=600)
            plt.close()

            analysis_df_neigborhood.to_csv(output_folder + '/' + 'Selected neighborhood analysis/' + select_neighboring_cell_type_and_state[1][0] + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + '/' + 'Selected neighborhood.csv')

#---------------------------------------------------------------------------------------------------------------------------------

    if analyse_states == True:

        for phenotype in phenotypes:

            columns_final_values = ['FileName', 'MetaData']
            df_final_values_states = pd.DataFrame(columns=columns_final_values)

            state_list = []

            if add_mixed_cells == True:

                for row in cell_type_state_and_neighbor_meta_data_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                            for state in type_and_state[1:]:
                                state_list.append(state)

                for row in cell_type_state_and_neighbor_meta_data_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    for type_and_state in cell_types_and_states:

                        if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                            for state in type_and_state[1:]:
                                state_list.append(state)

            if add_mixed_cells == False:

                for row in cell_type_state_and_neighbor_meta_data_1[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    if len(cell_types_and_states) == 2:

                        if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                            for state in cell_types_and_states[0][1:]:
                                state_list.append(state)

                for row in cell_type_state_and_neighbor_meta_data_2[df_column]:
                    cell_types_and_states = ast.literal_eval(row)

                    if len(cell_types_and_states) == 2:

                        if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                            for state in cell_types_and_states[0][1:]:
                                state_list.append(state)

            state_list = list(set(state_list))

            analysis_values_phenotype_states = ['FileName',
                                                'FileName_Image',
                                                'MetaData',
                                                'Phenotype_counter']

            for x in range(0, len(state_list)):
                analysis_values_phenotype_states.append(state_list[x])

            analysis_df_phenotype_states = pd.DataFrame(columns=analysis_values_phenotype_states)

            UniqueFileName = cell_type_state_and_neighbor_meta_data_1.ImageName.unique()
            DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

            for key in DataFrameDict_file.keys():
                DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_1[:][cell_type_state_and_neighbor_meta_data_1.ImageName == key].reset_index()

            for file in UniqueFileName:

                Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

                UniqueImage = Cells_in_file.ImageNumber.unique()
                DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                for key in DataFrameDict_image.keys():
                    DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index()

                for image in UniqueImage:

                    Cells_in_image = pd.DataFrame(DataFrameDict_image[image])
                    state_list_counter = [0 for x in range(0, len(state_list))]
                    Phenotype_counter = 0
                    analysis_df_row = []

                    for index, cell in Cells_in_image.iterrows():

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            if select_neighboring_cell_type_and_state[0] == False:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Phenotype_counter = Phenotype_counter + 1

                                        for state in type_and_state[1:]:
                                            state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                            if select_neighboring_cell_type_and_state[0] == True:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter > 0:

                                            Phenotype_counter = Phenotype_counter + 1

                                            for state in type_and_state[1:]:
                                                state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                        if add_mixed_cells == False:

                            if select_neighboring_cell_type_and_state[0] == False:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        Phenotype_counter = Phenotype_counter + 1

                                        for state in cell_types_and_states[0][1:]:
                                            state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                            if select_neighboring_cell_type_and_state[0] == True:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter > 0:

                                            Phenotype_counter = Phenotype_counter + 1

                                            for state in neighboring_cell_types_and_states[0][1:]:
                                                state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    analysis_df_row.append(file)
                    analysis_df_row.append(file + '_' + str(image))
                    analysis_df_row.append(1)
                    analysis_df_row.append(Phenotype_counter)

                    for state_count in state_list_counter:
                        analysis_df_row.append(state_count)

                    analysis_df_phenotype_states.loc[len(analysis_df_phenotype_states)] = analysis_df_row

            UniqueFileName = cell_type_state_and_neighbor_meta_data_2.ImageName.unique()
            DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

            for key in DataFrameDict_file.keys():
                DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_2[:][cell_type_state_and_neighbor_meta_data_2.ImageName == key].reset_index()

            for file in UniqueFileName:

                Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

                UniqueImage = Cells_in_file.ImageNumber.unique()
                DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                for key in DataFrameDict_image.keys():
                    DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index()

                for image in UniqueImage:

                    Cells_in_image = pd.DataFrame(DataFrameDict_image[image])
                    state_list_counter = [0 for x in range(0, len(state_list))]
                    Phenotype_counter = 0
                    analysis_df_row = []

                    for index, cell in Cells_in_image.iterrows():

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True:

                            if select_neighboring_cell_type_and_state[0] == False:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        Phenotype_counter = Phenotype_counter + 1

                                        for state in type_and_state[1:]:
                                            state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                            if select_neighboring_cell_type_and_state[0] == True:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]
                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            for neighboring_type_and_state in neighboring_cell_types_and_states:

                                                if neighboring_type_and_state[0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_type_and_state[1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_type_and_state[1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter > 0:

                                            Phenotype_counter = Phenotype_counter + 1

                                            for state in type_and_state[1:]:
                                                state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                        if add_mixed_cells == False:

                            if select_neighboring_cell_type_and_state[0] == False:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        Phenotype_counter = Phenotype_counter + 1

                                        for state in cell_types_and_states[0][1:]:
                                            state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                            if select_neighboring_cell_type_and_state[0] == True:

                                state_list_counter_cell = [0 for x in range(0, len(state_list))]

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:

                                        neighborhood_counter = 0
                                        for neighboring_cell in ast.literal_eval(cell['Neighborhood']):

                                            neighboring_cell_types_and_states = ast.literal_eval(Cells_in_image.iloc[neighboring_cell][df_column])

                                            if len(neighboring_cell_types_and_states) == 2:

                                                if neighboring_cell_types_and_states[0][0] == select_neighboring_cell_type_and_state[1][0] and set(select_neighboring_cell_type_and_state[1][1]).issubset(set(neighboring_cell_types_and_states[0][1:])) is True and set(select_neighboring_cell_type_and_state[1][2]).isdisjoint(set(neighboring_cell_types_and_states[0][1:])) is True:
                                                    neighborhood_counter = neighborhood_counter + 1

                                        if neighborhood_counter > 0:

                                            Phenotype_counter = Phenotype_counter + 1

                                            for state in neighboring_cell_types_and_states[0][1:]:
                                                state_list_counter_cell[state_list.index(state)] = 1

                                state_list_counter = [state_list_counter[i] + state_list_counter_cell[i] for i in range(len(state_list_counter_cell))]

                    analysis_df_row.append(file)
                    analysis_df_row.append(file + '_' + str(image))
                    analysis_df_row.append(2)
                    analysis_df_row.append(Phenotype_counter)

                    for state_count in state_list_counter:
                        analysis_df_row.append(state_count)

                    analysis_df_phenotype_states.loc[len(analysis_df_phenotype_states)] = analysis_df_row

            if select_area == True:

                analysis_df_phenotype_states = analysis_df_phenotype_states[~analysis_df_phenotype_states['FileName_Image'].isin(remove_row_list)]

            for state in state_list:

                df_final_values_phenotype_column = []
                FileName = []
                MetaData = []

                plot_list_1 = []
                plot_list_2 = []

                if statistics == 'mannwhitneyu':

                    UniqueMeta_data = analysis_df_phenotype_states.MetaData.unique()
                    DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_phenotype_states.MetaData.unique()}

                    for key in DataFrameDict_Metadata.keys():
                        DataFrameDict_Metadata[key] = analysis_df_phenotype_states[:][analysis_df_phenotype_states.MetaData == key].reset_index()

                    for metadata in UniqueMeta_data:

                        Cells_in_metadata = DataFrameDict_Metadata[metadata]

                        if split_by == 'file':

                            UniqueFile = Cells_in_metadata.FileName.unique()
                            DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                            for key in DataFrameDict_file.keys():
                                DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                            for file in UniqueFile:

                                Cells_in_file = DataFrameDict_file[file]

                                sum_of_cells = sum(Cells_in_file['Phenotype_counter'].to_list())

                                sum_of_states = sum(Cells_in_file[state].to_list())

                                if metadata == 1:

                                    MetaData.append(1)

                                    if int(sum_of_states) > 0:
                                        plot_list_1.append((sum_of_states / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_states / sum_of_cells) * 100)
                                    if int(sum_of_states) == 0 or int(sum_of_cells) == 0:
                                        plot_list_1.append(0)
                                        df_final_values_phenotype_column.append(0)

                                if metadata == 2:

                                    MetaData.append(2)

                                    if int(sum_of_states) > 0:
                                        plot_list_2.append((sum_of_states / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_states / sum_of_cells) * 100)
                                    if int(sum_of_states) == 0 or int(sum_of_cells) == 0:
                                        plot_list_2.append(0)
                                        df_final_values_phenotype_column.append(0)

                                FileName.append(file)

                        if split_by == 'image':

                            UniqueImage = Cells_in_metadata.FileName_Image.unique()
                            DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                            for key in DataFrameDict_image.keys():
                                DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                            for image in UniqueImage:

                                Cells_in_image = DataFrameDict_image[image]

                                cells = Cells_in_image['Phenotype_counter']

                                selected_cells = Cells_in_image[state]

                                if metadata == 1:

                                    if int(selected_cells) > 0:
                                        plot_list_1.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_1.append(0)

                                if metadata == 2:

                                    if int(selected_cells) > 0:
                                        plot_list_2.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_2.append(0)

                if statistics == 'wilcoxon':

                    UniqueMeta_data = analysis_df_phenotype_states.MetaData.unique()
                    DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_phenotype_states.MetaData.unique()}

                    for key in DataFrameDict_Metadata.keys():
                        DataFrameDict_Metadata[key] = analysis_df_phenotype_states[:][analysis_df_phenotype_states.MetaData == key].reset_index()

                    for metadata in UniqueMeta_data:

                        Cells_in_metadata = DataFrameDict_Metadata[metadata]

                        if split_by == 'file':

                            if metadata == 1:

                                UniqueFile = Cells_in_metadata.FileName.unique()
                                DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                                for key in DataFrameDict_file.keys():
                                    DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index()

                                for file in UniqueFile:

                                    FileName.append(file)
                                    MetaData.append(1)

                                    Cells_in_file = DataFrameDict_file[file]

                                    sum_of_cells = sum(Cells_in_file['Phenotype_counter'].to_list())

                                    sum_of_states = sum(Cells_in_file[state].to_list())

                                    if int(sum_of_states) > 0:
                                        plot_list_1.append((sum_of_states / sum_of_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_states / sum_of_cells) * 100)
                                    if int(sum_of_states) == 0 or int(sum_of_cells) == 0:
                                        plot_list_1.append(0)
                                        df_final_values_phenotype_column.append(0)

                                    paired_file = sample_parings[file]
                                    paired_file_values = analysis_df_phenotype_states.loc[analysis_df_phenotype_states['FileName'] == paired_file]

                                    FileName.append(paired_file)
                                    MetaData.append(2)

                                    sum_of_paired_cells = sum(paired_file_values['Phenotype_counter'].to_list())

                                    sum_of_paired_states = sum(paired_file_values[state].to_list())

                                    if int(sum_of_paired_states) > 0:
                                        plot_list_2.append((sum_of_paired_states / sum_of_paired_cells) * 100)
                                        df_final_values_phenotype_column.append((sum_of_paired_states / sum_of_paired_cells) * 100)
                                    if int(sum_of_paired_states) == 0 or int(sum_of_paired_cells) == 0:
                                        plot_list_2.append(0)
                                        df_final_values_phenotype_column.append(0)

                            else:
                                continue

                        if split_by == 'image':

                            UniqueImage = Cells_in_metadata.FileName_Image.unique()
                            DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

                            for key in DataFrameDict_image.keys():
                                DataFrameDict_image[key] = Cells_in_metadata[:][Cells_in_metadata.FileName_Image == key].reset_index()

                            for image in UniqueImage:

                                Cells_in_image = DataFrameDict_image[image]

                                cells = Cells_in_image['Phenotype_counter']

                                selected_cells = Cells_in_image[state]

                                if metadata == 1:

                                    if int(selected_cells) > 0:
                                        plot_list_1.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_1.append(0)

                                if metadata == 2:

                                    if int(selected_cells) > 0:
                                        plot_list_2.append(float((selected_cells / cells) * 100))
                                    if int(selected_cells) == 0 or int(cells) == 0:
                                        plot_list_2.append(0)

                df_final_values_states[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + str(state)] = df_final_values_phenotype_column
                df_final_values_states['FileName'] = FileName
                df_final_values_states['MetaData'] = MetaData

                plot_list = [plot_list_1, plot_list_2]

                try:

                    if statistics == 'mannwhitneyu':
                        U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')

                    if statistics == 'wilcoxon':
                        res = wilcoxon(plot_list[0], plot_list[1])
                        p = res.pvalue

                    max_plot = max([max(plot) for plot in plot_list])
                    min_plot = min([min(plot) for plot in plot_list])

                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                    if statistics == 'wilcoxon' and split_by == 'file':

                        for x in range(0, len(plot_list[0])):
                            axs.plot([1, 2], [plot_list[0][x], plot_list[1][x]], color='gray', alpha=0.25)

                    axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c=colors[0], s=45)
                    axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c=colors[1], s=45)
                    axs.boxplot(plot_list, widths=0.4)

                    if 0.005 < p < 0.05:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.45, max_plot * 1.2, '*', fontsize=18, color='black')
                    elif 0.0005 < p < 0.005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.41, max_plot * 1.2, '**', fontsize=18, color='black')
                    elif 0.00005 < p < 0.0005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.37, max_plot * 1.2, '***', fontsize=18, color='black')
                    elif p < 0.00005:
                        axs.plot([1.25, 1.75], [max_plot * 1.15, max_plot * 1.15], color='black')
                        axs.text(1.34, max_plot * 1.2, '****', fontsize=18, color='black')

                    axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot * 1.35)
                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '\n' + str(state))
                    axs.set_xticks([1, 2], ['Metadata 1','Metadata 2'], rotation=-90)

                    if select_neighboring_cell_type_and_state[0] == False:

                        final_output_folder = output_folder + '/' + 'Phenotype states/' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' per ' + split_by + '/'

                        if os.path.isdir(final_output_folder) == True:
                            pass
                        else:
                            os.makedirs(final_output_folder)

                        axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per ' + split_by)
                        export_file_path = final_output_folder + '/' + str(state) + '.png'

                    if select_neighboring_cell_type_and_state[0] == True:

                        final_output_folder = output_folder + '/' + 'Phenotype states/' + phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' neighboring ' + str(select_neighboring_cell_type_and_state[1][0]) + ' pos ' + str(select_neighboring_cell_type_and_state[1][1]) + ' neg ' + str(select_neighboring_cell_type_and_state[1][1]) + ' per ' + split_by + '/'

                        if os.path.isdir(final_output_folder) == True:
                            pass
                        else:
                            os.makedirs(final_output_folder)

                        axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per ' + split_by)
                        export_file_path = final_output_folder + '/' + str(state) + '.png'

                    fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                    plt.close()

                except ValueError:
                    print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not express: ' + str(state))
                    continue

            analysis_df_phenotype_states.to_csv(final_output_folder + 'States.csv')
            df_final_values_states.to_csv(final_output_folder + 'State_plots.csv')

        plt.text(0.1, 0.1, 'p < 0.05 = *', fontsize=25, color='black')
        plt.text(0.1, 0.2, 'p < 0.005 = **', fontsize=25, color='black')
        plt.text(0.1, 0.3, 'p < 0.0005 = ***', fontsize=25, color='black')
        plt.text(0.1, 0.4, 'p < 0.00005 = ****', fontsize=25, color='black')
        plt.xlim(0, 0.5)
        plt.ylim(0, 0.5)

        plt.savefig(output_folder + '/' + 'Phenotype states/' + 'p_value_legend.png', bbox_inches='tight', dpi=600)
        plt.close()

#too think about
'''
    if split_by == 'cell':

        for phenotype in phenotypes:

            full_phenotype_neighbours_list_1 = []  # empty list for all the for all the neighborhood phenotypes
            full_phenotype_neighbours_list_2 = []  # empty list for all the for all the neighborhood phenotypes

            for key, value in DataFrameDict_Cells.items():

                Cells_per_image = pd.DataFrame(DataFrameDict_Cells[key])

                for index, cell in Cells_per_image.iterrows():

                    cell_types = ast.literal_eval(cell[df_column])

                    if add_mixed_cells == True and select_neighbors[0] == False:

                        for cell_type in cell_types:

                            if cell_type[0] == phenotype[0]:

                                pass_counter = 0
                                cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                for index, cell_n in Cells_per_image.iterrows():

                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                        index_counter = -1
                                        for surrounding_phenotype in surrounding_phenotypes:
                                            index_counter = index_counter + 1

                                            for neighboring_cell_type in neighboring_cell_types:

                                                if neighboring_cell_type[0] == surrounding_phenotype[0]:
                                                    cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                phenotype_neighbours_list_1.append(cell_type_neighbours)

                    if add_mixed_cells == False and select_neighbors[0] == False:

                        if len(cell_types) == 2 and cell_types[0][0] == phenotype[0]:

                            pass_counter = 0
                            cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                            cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                            for index, cell_n in Cells_per_image.iterrows():

                                if int(cell_n['Cell_number']) in cell_neigboorhood:

                                    neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                    index_counter = -1
                                    for surrounding_phenotype in surrounding_phenotypes:
                                        index_counter = index_counter + 1

                                        if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == surrounding_phenotype[0]:
                                            cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                            phenotype_neighbours_list_1.append(cell_type_neighbours)

                    if add_mixed_cells == True and select_neighbors[0] == True:

                        for cell_type in cell_types:

                            if cell_type[0] == phenotype[0]:

                                pass_counter = 0
                                cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                for index, cell_n in Cells_per_image.iterrows():

                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                        index_counter = -1
                                        for surrounding_phenotype in surrounding_phenotypes:
                                            index_counter = index_counter + 1

                                            for neighboring_cell_type in neighboring_cell_types:

                                                if neighboring_cell_type[0] == surrounding_phenotype[0]:
                                                    cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                                if neighboring_cell_type[0] == select_neighbors[1]:
                                                    pass_counter = pass_counter + 1
                                if pass_counter > 0:
                                    phenotype_neighbours_list_1.append(cell_type_neighbours)

                                else:
                                    continue

                    if add_mixed_cells == False and select_neighbors[0] == True:

                        if len(cell_types) == 2 and cell_types[0][0] == phenotype[0]:

                            pass_counter = 0
                            cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                            cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                            for index, cell_n in Cells_per_image.iterrows():

                                if int(cell_n['Cell_number']) in cell_neigboorhood:

                                    neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                    index_counter = -1
                                    for surrounding_phenotype in surrounding_phenotypes:
                                        index_counter = index_counter + 1

                                        if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == surrounding_phenotype[0]:
                                            cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                        if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == select_neighbors[1]:
                                            pass_counter = pass_counter + 1

                            if pass_counter > 0:
                                phenotype_neighbours_list_1.append(cell_type_neighbours)

                            else:
                                continue

            if filename in meta_data_list_2:
                print('Running neighborhood statistics for ' + str(phenotype) + ' : ' + str(filedir))
                Cells = pd.read_csv(filedir)
                Cells['Meta_data_category'] = '2'

                for key in DataFrameDict_Cells.keys():
                    DataFrameDict_Cells[key] = Cells[:][Cells.ImageNumber == key].reset_index()

                for key, value in DataFrameDict_Cells.items():

                    Cells_per_image = pd.DataFrame(DataFrameDict_Cells[key])

                    for index, cell in Cells_per_image.iterrows():

                        cell_types = ast.literal_eval(cell[df_column])

                        if add_mixed_cells == True and select_neighbors[0] == False:

                            for cell_type in cell_types:

                                if cell_type[0] == phenotype[0]:

                                    pass_counter = 0
                                    cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                    cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                    for index, cell_n in Cells_per_image.iterrows():

                                        if int(cell_n['Cell_number']) in cell_neigboorhood:

                                            neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                            index_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:
                                                index_counter = index_counter + 1

                                                for neighboring_cell_type in neighboring_cell_types:

                                                    if neighboring_cell_type[0] == surrounding_phenotype[0]:
                                                        cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                    phenotype_neighbours_list_2.append(cell_type_neighbours)

                        if add_mixed_cells == False and select_neighbors[0] == False:

                            if len(cell_types) == 2 and cell_types[0][0] == phenotype[0]:

                                pass_counter = 0
                                cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                for index, cell_n in Cells_per_image.iterrows():

                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                        index_counter = -1
                                        for surrounding_phenotype in surrounding_phenotypes:
                                            index_counter = index_counter + 1

                                            if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == surrounding_phenotype[0]:
                                                cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                phenotype_neighbours_list_2.append(cell_type_neighbours)

                        if add_mixed_cells == True and select_neighbors[0] == True:

                            for cell_type in cell_types:

                                if cell_type[0] == phenotype[0]:

                                    pass_counter = 0
                                    cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                    cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                    for index, cell_n in Cells_per_image.iterrows():

                                        if int(cell_n['Cell_number']) in cell_neigboorhood:

                                            neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                            index_counter = -1
                                            for surrounding_phenotype in surrounding_phenotypes:
                                                index_counter = index_counter + 1

                                                for neighboring_cell_type in neighboring_cell_types:

                                                    if neighboring_cell_type[0] == surrounding_phenotype[0]:
                                                        cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                                    if neighboring_cell_type[0] == select_neighbors[1]:
                                                        pass_counter = pass_counter + 1
                                    if pass_counter > 0:
                                        phenotype_neighbours_list_2.append(cell_type_neighbours)

                                    else:
                                        continue

                        if add_mixed_cells == False and select_neighbors[0] == True:

                            if len(cell_types) == 2 and cell_types[0][0] == phenotype[0]:

                                pass_counter = 0
                                cell_type_neighbours = [0 for x in range(len(surrounding_phenotypes))]
                                cell_neigboorhood = ast.literal_eval(cell['Neighborhood'])  # imports the neighborhood    #Neighborhood

                                for index, cell_n in Cells_per_image.iterrows():

                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        neighboring_cell_types = ast.literal_eval(cell_n[df_column])

                                        index_counter = -1
                                        for surrounding_phenotype in surrounding_phenotypes:
                                            index_counter = index_counter + 1

                                            if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == surrounding_phenotype[0]:
                                                cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                            if len(neighboring_cell_types) == 2 and neighboring_cell_types[0][0] == select_neighbors[1]:
                                                pass_counter = pass_counter + 1

                                if pass_counter > 0:
                                    phenotype_neighbours_list_2.append(cell_type_neighbours)

                                else:
                                    continue

        if len(phenotype_neighbours_list_1) >= 1:
            for element in phenotype_neighbours_list_1:
                full_phenotype_neighbours_list_1.append(element)

        if len(phenotype_neighbours_list_2) >= 1:
            for element in phenotype_neighbours_list_2:
                full_phenotype_neighbours_list_2.append(element)

        full_phenotype_neighbours_list_1_percentages = []  # empty list for the percentages
        for neighbors in full_phenotype_neighbours_list_1:

            if sum(neighbors) == 0:  # checks if there are neighbors otherwise we devide by 0
                neighbors_percentages = neighbors
            else:
                neighbors_percentages = [(number / sum(neighbors)) * 100 for number in neighbors]  # calculates the percentage

            full_phenotype_neighbours_list_1_percentages.append(neighbors_percentages)  # appends ther percentages

        full_phenotype_neighbours_list_2_percentages = []  # empty list for the percentages
        for neighbors in full_phenotype_neighbours_list_2:

            if sum(neighbors) == 0:  # checks if there are neighbors otherwise we devide by 0
                neighbors_percentages = neighbors
            else:
                neighbors_percentages = [(number / sum(neighbors)) * 100 for number in neighbors]  # calculates the percentage

            full_phenotype_neighbours_list_2_percentages.append(neighbors_percentages)  # appends ther percentages

        x = np.arange(len(surrounding_phenotypes))  # amount of surrounding phenotypes
        width = 0.5  # sets the width if the individual bars
        multiplier = 0  # start multiplier

        fig, ax = plt.subplots(figsize=(len(surrounding_phenotypes) * 1.5, 5))  # initializes the figure

        for cell_type_number in range(len(full_phenotype_neighbours_list_1_percentages[0])):

            neighboring_cell_percentages_1 = [cell_lines[cell_type_number] for cell_lines in full_phenotype_neighbours_list_1_percentages]  # splits the full list by phenotypes
            neighboring_cell_percentages_2 = [cell_lines[cell_type_number] for cell_lines in full_phenotype_neighbours_list_2_percentages]  # splits the full list by phenotypes

            if statistics == 'mannwhitneyu':
                U1, p = mannwhitneyu(neighboring_cell_percentages_1, neighboring_cell_percentages_2, method='auto')

            if statistics == 'wilcoxon':
                res = wilcoxon(neighboring_cell_percentages_1, neighboring_cell_percentages_2)
                p = res.pvalue

                # creates the violin plots
            offset = width * multiplier
            violin_parts_1 = ax.violinplot(neighboring_cell_percentages_1, positions=[cell_type_number + offset], showmedians=True)
            violin_parts_2 = ax.violinplot(neighboring_cell_percentages_2, positions=[cell_type_number + offset + 1], showmedians=True)

            # sets the violin plot colors and style
            for pc in violin_parts_1['bodies']:
                pc.set_edgecolor('#000000')
                pc.set_facecolor('#cc0000')
                pc.set_linewidth(0.5)
                pc.set_alpha(1)
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
                vp = violin_parts_1[partname]
                vp.set_edgecolor('#000000')
                vp.set_linewidth(1)
            for pc in violin_parts_2['bodies']:
                pc.set_edgecolor('#000000')
                pc.set_linewidth(0.5)
                pc.set_facecolor('#1338BE')
                pc.set_alpha(1)
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
                vp = violin_parts_2[partname]
                vp.set_edgecolor('#000000')
                vp.set_linewidth(1)

            # adds stars based on the p value
            if 0.005 < p < 0.05:
                ax.text(cell_type_number + offset + 0.4, 90, '*', fontsize=12, color='black')
            elif 0.0005 < p < 0.005:
                ax.text(cell_type_number + offset + 0.35, 90, '**', fontsize=12, color='black')
            elif 0.00005 < p < 0.0005:
                ax.text(cell_type_number + offset + 0.25, 90, '***', fontsize=12, color='black')
            elif p < 0.00005:
                ax.text(cell_type_number + offset + 0.2, 90, '****', fontsize=12, color='black')

            multiplier += 2  # adds 2 to the multiplier

        ax.set_ylabel('Percent')  # adds the y label
        ax.set_title('Neighbouring cells of ' + str(phenotype[0]) + ' cells')  # title for the plot
        ax.set_xticks(x * 2 + width, sum(surrounding_phenotypes, []), rotation=-90)  # ticks for the x axis legend

        # adds a legend
        red_patch = mpatches.Patch(color='#cc0000', label='Before')
        blue_patch = mpatches.Patch(color='#1338BE', label='After')
        ax.legend(handles=[red_patch, blue_patch], loc="upper right")

        # adds the star ledgend
        textstr = '\n'.join((r'p < 0.05 = *', r'p < 0.005 = **', r'p < 0.0005 = ***', r'p < 0.00005 = ****'))
        ax.text(0.75, 0.86, textstr, transform=ax.transAxes, fontsize=8)

        for xtick, color in zip(ax.get_xticklabels(), colors):
            xtick.set_color(color)  # sets tick color

        ax.set_ylim(-5, 100 + 20)  # sets the y axis for the plot

        # cecks if the directory exists
        if os.path.isdir(output_folder) == True:
            pass
        else:
            os.makedirs(output_folder)

        if select_neighbors[0] == False:

            export_file_path = output_folder + '/' + 'Neighbouring cells of ' + str(phenotype[0]) + ' cells (median)'  # creates the export file path
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)  # saves the figure

        elif select_neighbors[0] == True:

            export_file_path = output_folder + '/' + 'Neighbouring cells of ' + str(phenotype[0]) + ' cells (median) neighbouring ' + str(select_neighbors[1]) + ' cells'  # creates the export file path
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)  # saves the figure

    if compare_means[0] == True:

        if not compare_means[1]:
            mean_cols = [col for col in cell_type_state_and_neighbor_1.columns if 'MeanIntensity' in col]

        else:
            mean_cols = ['Intensity_MeanIntensity_' + protein for protein in compare_means[1]]
'''
