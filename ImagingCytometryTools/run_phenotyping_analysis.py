import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import scandir as sd
import os
import pandas as pd
import ast
import seaborn as sns
import statistics as stat
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from scipy.stats import spearmanr

from assign_phenotypes_and_metadata import meta_data_list_1
from assign_phenotypes_and_metadata import meta_data_list_2

from assign_phenotypes_and_metadata import sample_merger

'''
Runs the statistics on the different phenotypes and cell states.

Here, one can select between calculating the overall abundance or cell category (immune/tissue cell) of a phenotype.
In addition to that, the overall cellular state can be accessed, with or without specific neighboring cells.
'''

def run_phenotyping_analysis(directory,
                             filestring,
                             df_column,
                             output_folder,
                             phenotypes,
                             percent_of='all',
                             add_mixed_cells=False,
                             statistics='mannwhitneyu',
                             select_area=False,
                             QC=False,
                             colors=['#D01515','#142EFB']):

    if statistics == 'wilcoxon':
        from assign_phenotypes_and_metadata import sample_parings

    counter_1 = 0
    counter_2 = 0
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
                    Cells = Cells.assign(Meta_data_category=1)

                    if counter_1 == 0:
                        cell_type_state_and_neighbor_meta_data_1 = Cells

                    else:
                        cell_type_state_and_neighbor_meta_data_1 = pd.concat([cell_type_state_and_neighbor_meta_data_1, Cells], ignore_index=True)

                    counter_1 = counter_1 + 1

                if filestring in filename_string and filename_string in meta_data_list_2:

                    print('Processing: ' + filedir)

                    Cells = pd.read_csv(filedir)
                    Cells = Cells[Cells.columns.drop(list(Cells.filter(regex='pixel_positive_area_')))]
                    Cells = Cells.assign(Meta_data_category=2)

                    if counter_2 == 0:
                        cell_type_state_and_neighbor_meta_data_2 = Cells

                    else:
                        cell_type_state_and_neighbor_meta_data_2 = pd.concat([cell_type_state_and_neighbor_meta_data_2, Cells], ignore_index=True)

                    counter_2 = counter_2 + 1

    # ---------------------------------------------------------------------------------------------------------------------------------

    # ------------------------------------- creates the dataframe for the final analysis -------------------------------------

    analysis_values_phenotyping = ['Percent_of_all_per_image',
                                   'Percent_of_immune_per_image',
                                   'Percent_of_tissue_per_image',
                                   'FileName',
                                   'FileName_Image',
                                   'MetaData']

    if select_area == True:
        analysis_values_phenotyping.append('Area')

    lineage_data = [lineage[0] for lineage in phenotypes]
    phenotype_data = [lineage[0] + ' pos ' + str(lineage[1]) + ' neg ' + str(lineage[2]) for lineage in phenotypes]

    for x in range(0, len(phenotypes)):
        analysis_values_phenotyping.append(lineage_data[x])
        analysis_values_phenotyping.append(phenotype_data[x])

    analysis_df_phenotypes = pd.DataFrame(columns=analysis_values_phenotyping)

    # ------------------------------------- Analysis of metadata 1 dataframe -------------------------------------

    UniqueFileName = cell_type_state_and_neighbor_meta_data_1.ImageName.unique()
    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

    for key in DataFrameDict_file.keys():
        DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_1[:][cell_type_state_and_neighbor_meta_data_1.ImageName == key].reset_index(drop=True)

    for file in UniqueFileName:

        Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

        UniqueImage = Cells_in_file.ImageNumber.unique()
        DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

        for key in DataFrameDict_image.keys():
            DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index(drop=True)

        for image in UniqueImage:

            if select_area == False:
                analysis_row = []

            if select_area == True:
                analysis_row_periportal = []
                analysis_row_liver = []
                analysis_row_sinusoidal = []

            phenotype_counter = -1
            for phenotype in phenotypes:
                phenotype_counter = phenotype_counter + 1

                if phenotype_counter == 0 and select_area == False:
                    Percent_of_all_per_image_counter = 0
                    Percent_of_immune_per_image_counter = 0
                    Percent_of_tissue_per_image_counter = 0

                if select_area == False:
                    percent_of_lineage_per_image_counter = 0
                    phenotype_per_image_counter_counter = 0

                if phenotype_counter == 0 and select_area == True:
                    Percent_of_all_per_image_counter_periportal = 0
                    Percent_of_immune_per_image_counter_periportal = 0
                    Percent_of_tissue_per_image_counter_periportal = 0

                    Percent_of_all_per_image_counter_liver = 0
                    Percent_of_immune_per_image_counter_liver = 0
                    Percent_of_tissue_per_image_counter_liver = 0

                    Percent_of_all_per_image_counter_sinusoidal = 0
                    Percent_of_immune_per_image_counter_sinusoidal = 0
                    Percent_of_tissue_per_image_counter_sinusoidal = 0

                if select_area == True:
                    percent_of_lineage_per_image_counter_periportal = 0
                    phenotype_per_image_counter_counter_periportal = 0

                    percent_of_lineage_per_image_counter_liver = 0
                    phenotype_per_image_counter_counter_liver = 0

                    percent_of_lineage_per_image_counter_sinusoidal = 0
                    phenotype_per_image_counter_counter_sinusoidal= 0

                Cells_in_image = pd.DataFrame(DataFrameDict_image[image])

                for index, cell in Cells_in_image.iterrows():

                    if QC == True:
                        cell_types_and_states = ast.literal_eval(cell[df_column])
                        QC_check = 0
                        for type_and_state in cell_types_and_states:
                            if len(type_and_state[1:]) > 0:
                                QC_check = 1

                        if QC_check == 0:
                            continue

                        else:

                            if select_area == True:

                                if cell['Selected Areas'] == 'periportal':

                                    if add_mixed_cells == True:

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                        for type_and_state in cell_types_and_states:

                                            if type_and_state[0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                            if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                                phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                                    if add_mixed_cells == False:

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + 1

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + 1

                                        if len(cell_types_and_states) == 2:

                                            if cell_types_and_states[0][0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                            if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                                phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                                if cell['Selected Areas'] == 'liver':

                                    if add_mixed_cells == True:

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                        for type_and_state in cell_types_and_states:

                                            if type_and_state[0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                            if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                                phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                                    if add_mixed_cells == False:

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + 1

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + 1

                                        if len(cell_types_and_states) == 2:

                                            if cell_types_and_states[0][0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                            if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                                phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                                if cell['Selected Areas'] == 'sinusoidal':

                                    if add_mixed_cells == True:

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                        if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                        for type_and_state in cell_types_and_states:

                                            if type_and_state[0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                            if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                                phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                                    if add_mixed_cells == False:

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                            Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + 1

                                        if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                            Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + 1

                                        if len(cell_types_and_states) == 2:

                                            if cell_types_and_states[0][0] == phenotype[0]:
                                                percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                            if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                                phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                            if select_area == False:

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

                    if QC == False:

                        cell_types_and_states = ast.literal_eval(cell[df_column])

                        if select_area == True:

                            if cell['Selected Areas'] == 'periportal':

                                if add_mixed_cells == True:

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                    for type_and_state in cell_types_and_states:

                                        if type_and_state[0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                        if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                            phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                                if add_mixed_cells == False:

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + 1

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + 1

                                    if len(cell_types_and_states) == 2:

                                        if cell_types_and_states[0][0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                        if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                            phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                            if cell['Selected Areas'] == 'liver':

                                if add_mixed_cells == True:

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                    for type_and_state in cell_types_and_states:

                                        if type_and_state[0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                        if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                            phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                                if add_mixed_cells == False:

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + 1

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + 1

                                    if len(cell_types_and_states) == 2:

                                        if cell_types_and_states[0][0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                        if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                            phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                            if cell['Selected Areas'] == 'sinusoidal':

                                if add_mixed_cells == True:

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                    if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                    for type_and_state in cell_types_and_states:

                                        if type_and_state[0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                        if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                            phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                                if add_mixed_cells == False:

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                        Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + 1

                                    if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                        Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + 1

                                    if len(cell_types_and_states) == 2:

                                        if cell_types_and_states[0][0] == phenotype[0]:
                                            percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                        if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                            phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                        if select_area == False:

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

                if select_area == True:

                    if phenotype_counter == 0:

                        Percent_of_all_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + Percent_of_tissue_per_image_counter_periportal
                        Percent_of_all_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + Percent_of_tissue_per_image_counter_liver
                        Percent_of_all_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + Percent_of_tissue_per_image_counter_sinusoidal

                        analysis_row_periportal.append(Percent_of_all_per_image_counter_periportal)
                        analysis_row_periportal.append(Percent_of_immune_per_image_counter_periportal)
                        analysis_row_periportal.append(Percent_of_tissue_per_image_counter_periportal)
                        analysis_row_periportal.append(file)
                        analysis_row_periportal.append(file + '_' + str(image))
                        analysis_row_periportal.append(1)
                        analysis_row_periportal.append('periportal')

                        analysis_row_liver.append(Percent_of_all_per_image_counter_liver)
                        analysis_row_liver.append(Percent_of_immune_per_image_counter_liver)
                        analysis_row_liver.append(Percent_of_tissue_per_image_counter_liver)
                        analysis_row_liver.append(file)
                        analysis_row_liver.append(file + '_' + str(image))
                        analysis_row_liver.append(1)
                        analysis_row_liver.append('liver')

                        analysis_row_sinusoidal.append(Percent_of_all_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(Percent_of_immune_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(Percent_of_tissue_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(file)
                        analysis_row_sinusoidal.append(file + '_' + str(image))
                        analysis_row_sinusoidal.append(1)
                        analysis_row_sinusoidal.append('sinusoidal')

                    analysis_row_periportal.append(percent_of_lineage_per_image_counter_periportal)
                    analysis_row_periportal.append(phenotype_per_image_counter_counter_periportal)

                    analysis_row_liver.append(percent_of_lineage_per_image_counter_liver)
                    analysis_row_liver.append(phenotype_per_image_counter_counter_liver)

                    analysis_row_sinusoidal.append(percent_of_lineage_per_image_counter_sinusoidal)
                    analysis_row_sinusoidal.append(phenotype_per_image_counter_counter_sinusoidal)

                if select_area == False:

                    if phenotype_counter == 0:

                        Percent_of_all_per_image_counter = Percent_of_immune_per_image_counter + Percent_of_tissue_per_image_counter

                        analysis_row.append(Percent_of_all_per_image_counter)
                        analysis_row.append(Percent_of_immune_per_image_counter)
                        analysis_row.append(Percent_of_tissue_per_image_counter)
                        analysis_row.append(file)
                        analysis_row.append(file + '_' + str(image))
                        analysis_row.append(1)

                    analysis_row.append(percent_of_lineage_per_image_counter)
                    analysis_row.append(phenotype_per_image_counter_counter)

            if select_area == True:
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_periportal
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_liver
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_sinusoidal

            if select_area == False:
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row

    # ------------------------------------- Analysis of metadata 2 dataframe -------------------------------------

    UniqueFileName = cell_type_state_and_neighbor_meta_data_2.ImageName.unique()
    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFileName}

    for key in DataFrameDict_file.keys():
        DataFrameDict_file[key] = cell_type_state_and_neighbor_meta_data_2[:][cell_type_state_and_neighbor_meta_data_2.ImageName == key].reset_index(drop=True)

    for file in UniqueFileName:

        Cells_in_file = pd.DataFrame(DataFrameDict_file[file])

        UniqueImage = Cells_in_file.ImageNumber.unique()
        DataFrameDict_image = {elem: pd.DataFrame() for elem in UniqueImage}

        for key in DataFrameDict_image.keys():
            DataFrameDict_image[key] = Cells_in_file[:][Cells_in_file.ImageNumber == key].reset_index(drop=True)

        for image in UniqueImage:

            if select_area == False:
                analysis_row = []

            if select_area == True:
                analysis_row_periportal = []
                analysis_row_liver = []
                analysis_row_sinusoidal = []

            phenotype_counter = -1
            for phenotype in phenotypes:
                phenotype_counter = phenotype_counter + 1

                if phenotype_counter == 0 and select_area == False:
                    Percent_of_all_per_image_counter = 0
                    Percent_of_immune_per_image_counter = 0
                    Percent_of_tissue_per_image_counter = 0

                if select_area == False:
                    percent_of_lineage_per_image_counter = 0
                    phenotype_per_image_counter_counter = 0

                if phenotype_counter == 0 and select_area == True:
                    Percent_of_all_per_image_counter_periportal = 0
                    Percent_of_immune_per_image_counter_periportal = 0
                    Percent_of_tissue_per_image_counter_periportal = 0

                    Percent_of_all_per_image_counter_liver = 0
                    Percent_of_immune_per_image_counter_liver = 0
                    Percent_of_tissue_per_image_counter_liver = 0

                    Percent_of_all_per_image_counter_sinusoidal = 0
                    Percent_of_immune_per_image_counter_sinusoidal = 0
                    Percent_of_tissue_per_image_counter_sinusoidal = 0

                if select_area == True:
                    percent_of_lineage_per_image_counter_periportal = 0
                    phenotype_per_image_counter_counter_periportal = 0

                    percent_of_lineage_per_image_counter_liver = 0
                    phenotype_per_image_counter_counter_liver = 0

                    percent_of_lineage_per_image_counter_sinusoidal = 0
                    phenotype_per_image_counter_counter_sinusoidal= 0

                Cells_in_image = pd.DataFrame(DataFrameDict_image[image])

                for index, cell in Cells_in_image.iterrows():

                    cell_types_and_states = ast.literal_eval(cell[df_column])

                    if select_area == True:

                        if cell['Selected Areas'] == 'periportal':

                            if add_mixed_cells == True:

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + (len(cell_types_and_states) - 1)

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                        phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + 1

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_periportal = Percent_of_tissue_per_image_counter_periportal + 1

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_periportal = percent_of_lineage_per_image_counter_periportal + 1

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                        phenotype_per_image_counter_counter_periportal = phenotype_per_image_counter_counter_periportal + 1

                        if cell['Selected Areas'] == 'liver':

                            if add_mixed_cells == True:

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + (len(cell_types_and_states) - 1)

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                        phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + 1

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_liver = Percent_of_tissue_per_image_counter_liver + 1

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_liver = percent_of_lineage_per_image_counter_liver + 1

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                        phenotype_per_image_counter_counter_liver = phenotype_per_image_counter_counter_liver + 1

                        if cell['Selected Areas'] == 'sinusoidal':

                            if add_mixed_cells == True:

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                if len(cell_types_and_states) >= 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + (len(cell_types_and_states) - 1)

                                for type_and_state in cell_types_and_states:

                                    if type_and_state[0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                    if type_and_state[0] == phenotype[0] and set(phenotype[1]).issubset(set(type_and_state[1:])) is True and set(phenotype[2]).isdisjoint(set(type_and_state[1:])) is True:
                                        phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                            if add_mixed_cells == False:

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45+':
                                    Percent_of_immune_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + 1

                                if len(cell_types_and_states) == 2 and cell_types_and_states[-1][0] == 'CD45-':
                                    Percent_of_tissue_per_image_counter_sinusoidal = Percent_of_tissue_per_image_counter_sinusoidal + 1

                                if len(cell_types_and_states) == 2:

                                    if cell_types_and_states[0][0] == phenotype[0]:
                                        percent_of_lineage_per_image_counter_sinusoidal = percent_of_lineage_per_image_counter_sinusoidal + 1

                                    if cell_types_and_states[0][0] == phenotype[0] and set(phenotype[1]).issubset(set(cell_types_and_states[0][1:])) is True and set(phenotype[2]).isdisjoint(set(cell_types_and_states[0][1:])) is True:
                                        phenotype_per_image_counter_counter_sinusoidal = phenotype_per_image_counter_counter_sinusoidal + 1

                    if select_area == False:

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

                if select_area == True:

                    if phenotype_counter == 0:

                        Percent_of_all_per_image_counter_periportal = Percent_of_immune_per_image_counter_periportal + Percent_of_tissue_per_image_counter_periportal
                        Percent_of_all_per_image_counter_liver = Percent_of_immune_per_image_counter_liver + Percent_of_tissue_per_image_counter_liver
                        Percent_of_all_per_image_counter_sinusoidal = Percent_of_immune_per_image_counter_sinusoidal + Percent_of_tissue_per_image_counter_sinusoidal

                        analysis_row_periportal.append(Percent_of_all_per_image_counter_periportal)
                        analysis_row_periportal.append(Percent_of_immune_per_image_counter_periportal)
                        analysis_row_periportal.append(Percent_of_tissue_per_image_counter_periportal)
                        analysis_row_periportal.append(file)
                        analysis_row_periportal.append(file + '_' + str(image))
                        analysis_row_periportal.append(2)
                        analysis_row_periportal.append('periportal')

                        analysis_row_liver.append(Percent_of_all_per_image_counter_liver)
                        analysis_row_liver.append(Percent_of_immune_per_image_counter_liver)
                        analysis_row_liver.append(Percent_of_tissue_per_image_counter_liver)
                        analysis_row_liver.append(file)
                        analysis_row_liver.append(file + '_' + str(image))
                        analysis_row_liver.append(2)
                        analysis_row_liver.append('liver')

                        analysis_row_sinusoidal.append(Percent_of_all_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(Percent_of_immune_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(Percent_of_tissue_per_image_counter_sinusoidal)
                        analysis_row_sinusoidal.append(file)
                        analysis_row_sinusoidal.append(file + '_' + str(image))
                        analysis_row_sinusoidal.append(2)
                        analysis_row_sinusoidal.append('sinusoidal')

                    analysis_row_periportal.append(percent_of_lineage_per_image_counter_periportal)
                    analysis_row_periportal.append(phenotype_per_image_counter_counter_periportal)

                    analysis_row_liver.append(percent_of_lineage_per_image_counter_liver)
                    analysis_row_liver.append(phenotype_per_image_counter_counter_liver)

                    analysis_row_sinusoidal.append(percent_of_lineage_per_image_counter_sinusoidal)
                    analysis_row_sinusoidal.append(phenotype_per_image_counter_counter_sinusoidal)

                if select_area == False:

                    if phenotype_counter == 0:

                        Percent_of_all_per_image_counter = Percent_of_immune_per_image_counter + Percent_of_tissue_per_image_counter

                        analysis_row.append(Percent_of_all_per_image_counter)
                        analysis_row.append(Percent_of_immune_per_image_counter)
                        analysis_row.append(Percent_of_tissue_per_image_counter)
                        analysis_row.append(file)
                        analysis_row.append(file + '_' + str(image))
                        analysis_row.append(2)

                    analysis_row.append(percent_of_lineage_per_image_counter)
                    analysis_row.append(phenotype_per_image_counter_counter)

            if select_area == True:
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_periportal
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_liver
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row_sinusoidal

            if select_area == False:
                analysis_df_phenotypes.loc[len(analysis_df_phenotypes)] = analysis_row

    # ------------------------------------- Final phenotypical analysis -------------------------------------

    analysis_df_phenotypes = analysis_df_phenotypes.loc[:, ~analysis_df_phenotypes.columns.duplicated()].copy()

    print(analysis_df_phenotypes.to_string())

    if select_area == True:

        UniqueArea = analysis_df_phenotypes.Area.unique()
        DataFrameDict_Area = {elem: pd.DataFrame() for elem in analysis_df_phenotypes.Area.unique()}

        for key in DataFrameDict_Area.keys():
            DataFrameDict_Area[key] = analysis_df_phenotypes[:][analysis_df_phenotypes.Area == key].reset_index()

        columns_final_values = ['FileName', 'MetaData']
        df_final_values_periportal = pd.DataFrame(columns=columns_final_values)
        df_final_values_liver = pd.DataFrame(columns=columns_final_values)
        df_final_values_sinusoidal = pd.DataFrame(columns=columns_final_values)

        for phenotype in phenotypes:

            df_final_values_phenotype_column_periportal = []
            FileName_periportal = []
            MetaData_periportal = []
            plot_list_1_periportal = []
            plot_list_2_periportal = []

            df_final_values_phenotype_column_liver = []
            FileName_liver = []
            MetaData_liver = []
            plot_list_1_liver = []
            plot_list_2_liver = []

            df_final_values_phenotype_column_sinusoidal = []
            FileName_sinusoidal = []
            MetaData_sinusoidal = []
            plot_list_1_sinusoidal = []
            plot_list_2_sinusoidal = []

            if statistics == 'mannwhitneyu':

                UniqueMeta_data = DataFrameDict_Area['periportal'].MetaData.unique()
                DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in DataFrameDict_Area['periportal'].MetaData.unique()}

                for key in DataFrameDict_Metadata.keys():
                    DataFrameDict_Metadata[key] = DataFrameDict_Area['periportal'][:][DataFrameDict_Area['periportal'].MetaData == key].reset_index()

                for metadata in UniqueMeta_data:

                    Cells_in_metadata = DataFrameDict_Metadata[metadata]

                    UniqueFile = Cells_in_metadata.FileName.unique()
                    DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                    for key in DataFrameDict_file.keys():
                        DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index(drop=True)

                    for file in UniqueFile:

                        Cells_in_file_periportal = DataFrameDict_file[file]
                        Cells_in_file_liver = DataFrameDict_Area['liver'].loc[DataFrameDict_Area['liver']['FileName'] == file]
                        Cells_in_file_sinusoidal = DataFrameDict_Area['sinusoidal'].loc[DataFrameDict_Area['sinusoidal']['FileName'] == file]

                        if percent_of == 'all':
                            sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_all_per_image'].to_list())
                            sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_all_per_image'].to_list())
                            sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_all_per_image'].to_list())

                        if percent_of == 'immune':
                            sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_immune_per_image'].to_list())
                            sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_immune_per_image'].to_list())
                            sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_immune_per_image'].to_list())

                        if percent_of == 'tissue':
                            sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_tissue_per_image'].to_list())
                            sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_tissue_per_image'].to_list())
                            sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_tissue_per_image'].to_list())

                        if percent_of == 'lineage':
                            sum_of_cells_periportal = sum(Cells_in_file_periportal[phenotype[0]].to_list())
                            sum_of_cells_liver = sum(Cells_in_file_liver[phenotype[0]].to_list())
                            sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal[phenotype[0]].to_list())

                        phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                        sum_of_selected_cells_periportal = sum(Cells_in_file_periportal[phenotype_df].to_list())
                        sum_of_selected_cells_liver = sum(Cells_in_file_liver[phenotype_df].to_list())
                        sum_of_selected_cells_sinusoidal = sum(Cells_in_file_sinusoidal[phenotype_df].to_list())

                        if metadata == 1:

                            MetaData_periportal.append(1)
                            MetaData_liver.append(1)
                            MetaData_sinusoidal.append(1)

                            if int(sum_of_selected_cells_periportal) > 0:
                                plot_list_1_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_liver) > 0:
                                plot_list_1_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_sinusoidal) > 0:
                                plot_list_1_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_periportal) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_1_periportal.append(0)
                                df_final_values_phenotype_column_periportal.append(0)

                            if int(sum_of_selected_cells_liver) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_1_liver.append(0)
                                df_final_values_phenotype_column_liver.append(0)

                            if int(sum_of_selected_cells_sinusoidal) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_1_sinusoidal.append(0)
                                df_final_values_phenotype_column_sinusoidal.append(0)

                        if metadata == 2:

                            MetaData_periportal.append(2)
                            MetaData_liver.append(2)
                            MetaData_sinusoidal.append(2)

                            if int(sum_of_selected_cells_periportal) > 0:
                                plot_list_2_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_liver) > 0:
                                plot_list_2_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_sinusoidal) > 0:
                                plot_list_2_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_periportal) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_2_periportal.append(0)
                                df_final_values_phenotype_column_periportal.append(0)

                            if int(sum_of_selected_cells_liver) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_2_liver.append(0)
                                df_final_values_phenotype_column_liver.append(0)

                            if int(sum_of_selected_cells_sinusoidal) == 0 or int((sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) == 0:
                                plot_list_2_sinusoidal.append(0)
                                df_final_values_phenotype_column_sinusoidal.append(0)

                        FileName_periportal.append(file)
                        FileName_liver.append(file)
                        FileName_sinusoidal.append(file)

            if statistics == 'wilcoxon':

                UniqueMeta_data = DataFrameDict_Area['periportal'].MetaData.unique()
                DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in DataFrameDict_Area['periportal'].MetaData.unique()}

                for key in DataFrameDict_Metadata.keys():
                    DataFrameDict_Metadata[key] = DataFrameDict_Area['periportal'][:][DataFrameDict_Area['periportal'].MetaData == key].reset_index()

                for metadata in UniqueMeta_data:

                    Cells_in_metadata = DataFrameDict_Metadata[metadata]

                    if metadata == 1:

                        UniqueFile = Cells_in_metadata.FileName.unique()
                        DataFrameDict_file = {elem: pd.DataFrame() for elem in UniqueFile}

                        for key in DataFrameDict_file.keys():
                            DataFrameDict_file[key] = Cells_in_metadata[:][Cells_in_metadata.FileName == key].reset_index(drop=True)

                        for file in UniqueFile:

                            FileName_periportal.append(file)
                            MetaData_periportal.append(1)
                            FileName_liver.append(file)
                            MetaData_liver.append(1)
                            FileName_sinusoidal.append(file)
                            MetaData_sinusoidal.append(1)

                            Cells_in_file_periportal = DataFrameDict_file[file]
                            Cells_in_file_liver = DataFrameDict_Area['liver'].loc[DataFrameDict_Area['liver']['FileName'] == file]
                            Cells_in_file_sinusoidal = DataFrameDict_Area['sinusoidal'].loc[DataFrameDict_Area['sinusoidal']['FileName'] == file]

                            if percent_of == 'all':
                                sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_all_per_image'].to_list())
                                sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_all_per_image'].to_list())
                                sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_all_per_image'].to_list())

                            if percent_of == 'immune':
                                sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_immune_per_image'].to_list())
                                sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_immune_per_image'].to_list())
                                sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_immune_per_image'].to_list())

                            if percent_of == 'tissue':
                                sum_of_cells_periportal = sum(Cells_in_file_periportal['Percent_of_tissue_per_image'].to_list())
                                sum_of_cells_liver = sum(Cells_in_file_liver['Percent_of_tissue_per_image'].to_list())
                                sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal['Percent_of_tissue_per_image'].to_list())

                            if percent_of == 'lineage':
                                sum_of_cells_periportal = sum(Cells_in_file_periportal[phenotype[0]].to_list())
                                sum_of_cells_liver = sum(Cells_in_file_liver[phenotype[0]].to_list())
                                sum_of_cells_sinusoidal = sum(Cells_in_file_sinusoidal[phenotype[0]].to_list())

                            phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                            sum_of_selected_cells_periportal = sum(Cells_in_file_periportal[phenotype_df].to_list())
                            sum_of_selected_cells_liver = sum(Cells_in_file_liver[phenotype_df].to_list())
                            sum_of_selected_cells_sinusoidal = sum(Cells_in_file_sinusoidal[phenotype_df].to_list())

                            if int(sum_of_selected_cells_periportal) == 0 or int(sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal) == 0:
                                plot_list_1_periportal.append(0)
                                df_final_values_phenotype_column_periportal.append(0)

                            if int(sum_of_selected_cells_liver) == 0 or int(sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal) == 0:
                                plot_list_1_liver.append(0)
                                df_final_values_phenotype_column_liver.append(0)

                            if int(sum_of_selected_cells_sinusoidal) == 0 or int(sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal) == 0:
                                plot_list_1_sinusoidal.append(0)
                                df_final_values_phenotype_column_sinusoidal.append(0)

                            if int(sum_of_selected_cells_periportal) > 0:
                                plot_list_1_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_periportal.append((sum_of_selected_cells_periportal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_liver) > 0:
                                plot_list_1_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_liver.append((sum_of_selected_cells_liver / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_cells_sinusoidal) > 0:
                                plot_list_1_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_sinusoidal.append((sum_of_selected_cells_sinusoidal / (sum_of_cells_periportal + sum_of_cells_liver + sum_of_cells_sinusoidal)) * 100)

                            paired_file = sample_parings[file]

                            paired_file_values_periportal = DataFrameDict_Area['periportal'].loc[DataFrameDict_Area['periportal']['FileName'] == paired_file]
                            paired_file_values_liver = DataFrameDict_Area['liver'].loc[DataFrameDict_Area['liver']['FileName'] == paired_file]
                            paired_file_values_sinusoidal = DataFrameDict_Area['sinusoidal'].loc[DataFrameDict_Area['sinusoidal']['FileName'] == paired_file]

                            FileName_periportal.append(paired_file)
                            MetaData_periportal.append(2)
                            FileName_liver.append(paired_file)
                            MetaData_liver.append(2)
                            FileName_sinusoidal.append(paired_file)
                            MetaData_sinusoidal.append(2)

                            if percent_of == 'all':
                                sum_of_paired_cells_periportal = sum(paired_file_values_periportal['Percent_of_all_per_image'].to_list())
                                sum_of_paired_cells_liver = sum(paired_file_values_liver['Percent_of_all_per_image'].to_list())
                                sum_of_paired_cells_sinusoidal = sum(paired_file_values_sinusoidal['Percent_of_all_per_image'].to_list())

                            if percent_of == 'immune':
                                sum_of_paired_cells_periportal = sum(paired_file_values_periportal['Percent_of_immune_per_image'].to_list())
                                sum_of_paired_cells_liver = sum(paired_file_values_liver['Percent_of_immune_per_image'].to_list())
                                sum_of_paired_cells_sinusoidal = sum(paired_file_values_sinusoidal['Percent_of_immune_per_image'].to_list())

                            if percent_of == 'tissue':
                                sum_of_paired_cells_periportal = sum(paired_file_values_periportal['Percent_of_tissue_per_image'].to_list())
                                sum_of_paired_cells_liver = sum(paired_file_values_liver['Percent_of_tissue_per_image'].to_list())
                                sum_of_paired_cells_sinusoidal = sum(paired_file_values_sinusoidal['Percent_of_tissue_per_image'].to_list())

                            if percent_of == 'lineage':
                                sum_of_paired_cells_periportal = sum(paired_file_values_periportal[phenotype[0]].to_list())
                                sum_of_paired_cells_liver = sum(paired_file_values_liver[phenotype[0]].to_list())
                                sum_of_paired_cells_sinusoidal = sum(paired_file_values_sinusoidal[phenotype[0]].to_list())

                            phenotype_df = str(phenotype[0] + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))

                            sum_of_selected_paired_cells_periportal = sum(paired_file_values_periportal[phenotype_df].to_list())
                            sum_of_selected_paired_cells_liver = sum(paired_file_values_liver[phenotype_df].to_list())
                            sum_of_selected_paired_cells_sinusoidal = sum(paired_file_values_sinusoidal[phenotype_df].to_list())

                            if int(sum_of_selected_paired_cells_periportal) == 0 or int(sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal) == 0:
                                plot_list_2_periportal.append(0)
                                df_final_values_phenotype_column_periportal.append(0)

                            if int(sum_of_selected_paired_cells_liver) == 0 or int(sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal) == 0:
                                plot_list_2_liver.append(0)
                                df_final_values_phenotype_column_liver.append(0)

                            if int(sum_of_selected_paired_cells_sinusoidal) == 0 or int(sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal) == 0:
                                plot_list_2_sinusoidal.append(0)
                                df_final_values_phenotype_column_sinusoidal.append(0)

                            if int(sum_of_selected_paired_cells_periportal) > 0:
                                plot_list_2_periportal.append((sum_of_selected_paired_cells_periportal / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_periportal.append((sum_of_selected_paired_cells_periportal / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_paired_cells_liver) > 0:
                                plot_list_2_liver.append((sum_of_selected_paired_cells_liver / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_liver.append((sum_of_selected_paired_cells_liver / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)

                            if int(sum_of_selected_paired_cells_sinusoidal) > 0:
                                plot_list_2_sinusoidal.append((sum_of_selected_paired_cells_sinusoidal / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)
                                df_final_values_phenotype_column_sinusoidal.append((sum_of_selected_paired_cells_sinusoidal / (sum_of_paired_cells_periportal + sum_of_paired_cells_liver + sum_of_paired_cells_sinusoidal)) * 100)

                    else:
                        continue

            df_final_values_periportal[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])] = df_final_values_phenotype_column_periportal
            df_final_values_periportal['FileName'] = FileName_periportal
            df_final_values_periportal['MetaData'] = MetaData_periportal
            plot_list_periportal = [plot_list_1_periportal, plot_list_2_periportal]

            df_final_values_liver[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])] = df_final_values_phenotype_column_liver
            df_final_values_liver['FileName'] = FileName_liver
            df_final_values_liver['MetaData'] = MetaData_liver
            plot_list_liver = [plot_list_1_liver, plot_list_2_liver]

            df_final_values_sinusoidal[str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2])] = df_final_values_phenotype_column_sinusoidal
            df_final_values_sinusoidal['FileName'] = FileName_sinusoidal
            df_final_values_sinusoidal['MetaData'] = MetaData_sinusoidal
            plot_list_sinusoidal = [plot_list_1_sinusoidal, plot_list_2_sinusoidal]

            try:

                if statistics == 'mannwhitneyu':
                    U1, p = mannwhitneyu(plot_list_periportal[0], plot_list_periportal[1], method='auto')

                if statistics == 'wilcoxon':
                    result = wilcoxon(plot_list_periportal[0], plot_list_periportal[1])
                    p = result.pvalue

                max_plot = max([max(plot) for plot in plot_list_periportal])
                min_plot = min([min(plot) for plot in plot_list_periportal])

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                if statistics == 'wilcoxon':

                    for x in range(0, len(plot_list_periportal[0])):
                        axs.plot([1, 2], [plot_list_periportal[0][x], plot_list_periportal[1][x]], color='gray', alpha=0.25)

                axs.scatter([1 for x in range(len(plot_list_periportal[0]))], plot_list_periportal[0], c=colors[0], s=45)
                axs.scatter([2 for x in range(len(plot_list_periportal[1]))], plot_list_periportal[1], c=colors[1], s=45)
                axs.boxplot(plot_list_periportal, widths=0.4)

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
                axs.set_xticks([1, 2], ['Metadata 1', 'Metadata 2'], rotation=-90)


                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_periportal.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_periportal.png'

                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                plt.close()

            except ValueError:
                print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not exist')
                continue

            # --------------------------------------------

            try:

                if statistics == 'mannwhitneyu':
                    U1, p = mannwhitneyu(plot_list_liver[0], plot_list_liver[1], method='auto')

                if statistics == 'wilcoxon':
                    result = wilcoxon(plot_list_liver[0], plot_list_liver[1])
                    p = result.pvalue

                max_plot = max([max(plot) for plot in plot_list_liver])
                min_plot = min([min(plot) for plot in plot_list_liver])

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                if statistics == 'wilcoxon':

                    for x in range(0, len(plot_list_liver[0])):
                        axs.plot([1, 2], [plot_list_liver[0][x], plot_list_liver[1][x]], color='gray', alpha=0.25)

                axs.scatter([1 for x in range(len(plot_list_liver[0]))], plot_list_liver[0], c=colors[0], s=45)
                axs.scatter([2 for x in range(len(plot_list_liver[1]))], plot_list_liver[1], c=colors[1], s=45)
                axs.boxplot(plot_list_liver, widths=0.4)

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
                axs.set_xticks([1, 2], ['Metadata 1', 'Metadata 2'], rotation=-90)

                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_liver.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_liver.png'

                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                plt.close()

            except ValueError:
                print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not exist')
                continue

            # --------------------------------------------

            try:

                if statistics == 'mannwhitneyu':
                    U1, p = mannwhitneyu(plot_list_sinusoidal[0], plot_list_sinusoidal[1], method='auto')

                if statistics == 'wilcoxon':
                    result = wilcoxon(plot_list_sinusoidal[0], plot_list_sinusoidal[1])
                    p = result.pvalue

                max_plot = max([max(plot) for plot in plot_list_sinusoidal])
                min_plot = min([min(plot) for plot in plot_list_sinusoidal])

                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

                if statistics == 'wilcoxon':

                    for x in range(0, len(plot_list_sinusoidal[0])):
                        axs.plot([1, 2], [plot_list_sinusoidal[0][x], plot_list_sinusoidal[1][x]], color='gray', alpha=0.25)

                axs.scatter([1 for x in range(len(plot_list_sinusoidal[0]))], plot_list_sinusoidal[0], c=colors[0], s=45)
                axs.scatter([2 for x in range(len(plot_list_sinusoidal[1]))], plot_list_sinusoidal[1], c=colors[1], s=45)
                axs.boxplot(plot_list_sinusoidal, widths=0.4)

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
                axs.set_xticks([1, 2], ['Metadata 1', 'Metadata 2'], rotation=-90)

                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_sinusoidal.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '_sinusoidal.png'

                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)
                plt.close()

            except ValueError:
                print('The Phenotype: ' + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + ' does not exist')
                continue
    
    if select_area == False:

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

            if statistics == 'wilcoxon':

                UniqueMeta_data = analysis_df_phenotypes.MetaData.unique()
                DataFrameDict_Metadata = {elem: pd.DataFrame() for elem in analysis_df_phenotypes.MetaData.unique()}

                for key in DataFrameDict_Metadata.keys():
                    DataFrameDict_Metadata[key] = analysis_df_phenotypes[:][analysis_df_phenotypes.MetaData == key].reset_index()

                for metadata in UniqueMeta_data:

                    Cells_in_metadata = DataFrameDict_Metadata[metadata]

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

                if statistics == 'wilcoxon':

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
                axs.set_xticks([1, 2], ['Metadata 1', 'Metadata 2'], rotation=-90)
                
                if percent_of == 'all' or percent_of == 'immune' or percent_of == 'tissue':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + percent_of + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + percent_of + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '.png'

                if percent_of == 'lineage':

                    final_output_folder = output_folder + '/' + 'Phenotyping/' + 'Percent of ' + str(phenotype[0]) + ' cells per file/'

                    if os.path.isdir(final_output_folder) == True:
                        pass
                    else:
                        os.makedirs(final_output_folder)

                    axs.set_title(str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]))
                    axs.set_ylabel('Percent of ' + phenotype[0] + ' cells per file')
                    export_file_path = final_output_folder + str(phenotype[0]) + ' pos ' + str(phenotype[1]) + ' neg ' + str(phenotype[2]) + '.png'


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

    analysis_df_phenotypes.to_csv(final_output_folder + 'Phenotypes.csv')

    if select_area == False:
        df_final_values.to_csv(final_output_folder + 'Phenotype_plots.csv')

    if select_area == True:
        df_final_values_periportal.to_csv(final_output_folder + 'Phenotype_plots_periportal.csv')
        df_final_values_liver.to_csv(final_output_folder + 'Phenotype_plots_liver.csv')
        df_final_values_sinusoidal.to_csv(final_output_folder + 'Phenotype_plots_sinusoidal.csv')

        print(df_final_values_periportal.head(30).to_string())
        print('-')
        print(df_final_values_liver.head(30).to_string())
        print('-')
        print(df_final_values_sinusoidal.head(30).to_string())
        print('-')