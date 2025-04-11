import scandir as sd
import os
import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu

from assign_stuff import meta_data_list_1
from assign_stuff import meta_data_list_2

'''
Runs the statistics on the different cell types and they neighborhood.
'''

def run_neighborhood_statistics_ISO(directory, filestring, df_column, output_folder, phenotypes, surrounding_phenotypes, colors, add_mixed_cells = False, select_neighbors = [False]):

    #goes throw all the selected phenotypes
    for phenotype in phenotypes:

        full_phenotype_neighbours_list_1 = [] #empty list for all the for all the neighborhood phenotypes
        full_phenotype_neighbours_list_2 = [] #empty list for all the for all the neighborhood phenotypes

        for paths, dirs, files in sd.walk(directory):  # goes throw all files and folders in given directory
    
            for file in os.listdir(paths):  # goes throw all files in a folder
                filedir = os.path.join(paths, file)  # returns full file directory
    
                if filedir.endswith(".csv"):  # returns all files that end with .csv
    
                    # gives you the file name
                    filename = os.path.basename(file)
                    filename_string = str(filename)
    
                    if filestring in filename_string:  # returns all files that contain the provided file string

                        phenotype_neighbours_list_1 = [] # sub list for an individual file
                        phenotype_neighbours_list_2 = [] # sub list for an individual file

                        if filename in meta_data_list_1:
                            print('Running neighborhood statistics for ' + str(phenotype) + ' : ' + str(filedir))

                            Cells = pd.read_csv(filedir) #imports the file
                            Cells['Meta_data_category'] = '1' #adds the metadata category into the data frame

                            UniqueNames_Cells = Cells.ImageNumber.unique()  # creates a list of all unique images
                            DataFrameDict_Cells = {elem: pd.DataFrame() for elem in UniqueNames_Cells}  # creates a dictionary of all unique images

                            for key in DataFrameDict_Cells.keys():
                                DataFrameDict_Cells[key] = Cells[:][Cells.ImageNumber == key].reset_index() #resets the index for the data frame dictionary

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

                                                            if len(neighboring_cell_types)== 2 and neighboring_cell_types[0][0] == surrounding_phenotype[0]:
                                                                cell_type_neighbours[index_counter] = cell_type_neighbours[index_counter] + 1

                                                            if len(neighboring_cell_types)== 2 and neighboring_cell_types[0][0] == select_neighbors[1]:
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

        full_phenotype_neighbours_list_1_percentages = [] #empty list for the percentages
        for neighbors in full_phenotype_neighbours_list_1:

            if sum(neighbors) == 0: #checks if there are neighbors otherwise we devide by 0
                neighbors_percentages = neighbors
            else:
                neighbors_percentages = [(number / sum(neighbors)) * 100 for number in neighbors] #calculates the percentage
                
            full_phenotype_neighbours_list_1_percentages.append(neighbors_percentages) #appends ther percentages

        full_phenotype_neighbours_list_2_percentages = [] #empty list for the percentages
        for neighbors in full_phenotype_neighbours_list_2:

            if sum(neighbors) == 0: #checks if there are neighbors otherwise we devide by 0
                neighbors_percentages = neighbors
            else:
                neighbors_percentages = [(number / sum(neighbors)) * 100 for number in neighbors] #calculates the percentage
                
            full_phenotype_neighbours_list_2_percentages.append(neighbors_percentages) #appends ther percentages

        x = np.arange(len(surrounding_phenotypes)) #amount of surrounding phenotypes
        width = 0.5 # sets the width if the individual bars
        multiplier = 0 # start multiplier

        fig, ax = plt.subplots(figsize=(len(surrounding_phenotypes) * 1.5, 5)) # initializes the figure

        for cell_type_number in range(len(full_phenotype_neighbours_list_1_percentages[0])):

            neighboring_cell_percentages_1 = [cell_lines[cell_type_number] for cell_lines in full_phenotype_neighbours_list_1_percentages] # splits the full list by phenotypes
            neighboring_cell_percentages_2 = [cell_lines[cell_type_number] for cell_lines in full_phenotype_neighbours_list_2_percentages] # splits the full list by phenotypes

            U1, p = mannwhitneyu(neighboring_cell_percentages_1, neighboring_cell_percentages_2, method='auto') #runs the mann whitney u test

            #creates the violin plots
            offset = width * multiplier
            violin_parts_1 = ax.violinplot(neighboring_cell_percentages_1, positions=[cell_type_number + offset], showmedians=True)
            violin_parts_2 = ax.violinplot(neighboring_cell_percentages_2, positions=[cell_type_number + offset + 1], showmedians=True)

            #sets the violin plot colors and style
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

            #adds stars based on the p value
            if 0.005 < p < 0.05:
                ax.text(cell_type_number + offset + 0.4, 90, '*', fontsize=12, color='black')
            elif 0.0005 < p < 0.005:
                ax.text(cell_type_number + offset + 0.35, 90, '**', fontsize=12, color='black')
            elif 0.00005 < p < 0.0005:
                ax.text(cell_type_number + offset + 0.25, 90, '***', fontsize=12, color='black')
            elif p < 0.00005:
                ax.text(cell_type_number + offset + 0.2, 90, '****', fontsize=12, color='black')

            multiplier += 2 #adds 2 to the multiplier

        ax.set_ylabel('Percent') #adds the y label
        ax.set_title('Neighbouring cells of ' + str(phenotype[0]) + ' cells') #title for the plot
        ax.set_xticks(x * 2 + width, sum(surrounding_phenotypes, []), rotation=-90) #ticks for the x axis legend

        #adds a legend
        red_patch = mpatches.Patch(color='#cc0000', label='Before')
        blue_patch = mpatches.Patch(color='#1338BE', label='After')
        ax.legend(handles=[red_patch, blue_patch], loc="upper right")

        #adds the star ledgend
        textstr = '\n'.join((r'p < 0.05 = *', r'p < 0.005 = **', r'p < 0.0005 = ***', r'p < 0.00005 = ****'))
        ax.text(0.75, 0.86, textstr, transform=ax.transAxes, fontsize=8)

        for xtick, color in zip(ax.get_xticklabels(), colors):
            xtick.set_color(color) #sets tick color

        ax.set_ylim(-5, 100 + 20) # sets the y axis for the plot

        #cecks if the directory exists
        if os.path.isdir(output_folder) == True:
            pass
        else:
            os.makedirs(output_folder)

        if select_neighbors[0] == False:

            export_file_path = output_folder + '/' + 'Neighbouring cells of ' + str(phenotype[0]) + ' cells (median)' #creates the export file path
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600) #saves the figure

        elif select_neighbors[0] == True:

            export_file_path = output_folder + '/' + 'Neighbouring cells of ' + str(phenotype[0]) + ' cells (median) neighbouring ' + str(select_neighbors[1]) + ' cells'  # creates the export file path
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)  # saves the figure