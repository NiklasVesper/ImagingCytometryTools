import scandir as sd
import os
import pandas as pd
import ast
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from assign_stuff import meta_data_list_1
from assign_stuff import meta_data_list_2

'''
Runs the statistics on the different phenotypes and cell states.

Here, one can select between calculating the overall abundance or cell category (immune/tissue cell) of a phenotype.
In addition to that, the overall cellular state can be accessed, with or without specific neighboring cells.
'''

def run_phenotype_statistics_ISO(directory, filestring, df_column, output_folder, cell_types, percent_of = ['All'], add_mixed_cells = False, select_neighbors = [False]):

    df_meta_data_1 = pd.DataFrame([1]) #empty data frame for the file
    df_meta_data_1_column_list = [] #empty list for the file names

    df_meta_data_2 = pd.DataFrame([1]) #empty data frame for the file
    df_meta_data_2_column_list = [] #empty list for the file names
    
    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths) : #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                if filestring in filename_string: #returns all files that contain the provided file string

                    #checks whether the file is in the metadata list
                    if filename in meta_data_list_1:
                        print('Running phenotype statistics for: ' + str(filedir))
                        Cells = pd.read_csv(filedir) #imports the file
                        Cells['Meta_data_category'] = '1' #adds the metadata category into the data frame

                        df_meta_data_1_column_list.append(filename_string) #appends the file name
                        
                        if select_neighbors[0] == True and percent_of[0] == 'lineage':

                            Selected_Cells = pd.DataFrame().reindex_like(Cells) #data frame of the cells selected for a specific neighbor
                            Selected_Cells.drop(Selected_Cells.index, inplace=True)

                            #goes throw each cell
                            for index, cell in Cells.iterrows():

                                cell_neigboorhood = ast.literal_eval(cell['Neigboorhood']) #imports the neighborhood    #Neighborhood
                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']] #selects the cells in the same image

                                #goes throw each cell in the image
                                for index, cell_n in Cells_n.iterrows():

                                    #selects the neighboring cells
                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        cell_n_types = ast.literal_eval(cell_n['cell_types_and_states']) #imports the neighboring cell type

                                        # adds the cells with a specific neighbor to the new data frame
                                        if add_mixed_cells == False:
                                            if len(cell_n_types) == 2 and select_neighbors[1] == cell_n_types[0][0]:
                                                Selected_Cells.loc[len(Selected_Cells)] = cell

                                        elif add_mixed_cells == True:
                                            for cell_n_type in cell_n_types:
                                                if select_neighbors[1] == cell_n_type[0]:
                                                    Selected_Cells.loc[len(Selected_Cells)] = cell

                            phenotype_df = Selected_Cells[df_column]  # selects the analysis column
                            df_meta_data_1 = pd.concat([df_meta_data_1, phenotype_df], axis=1, ignore_index=True)  # adds the file data frame into the analysis dataframe

                        else:
                            phenotype_df = Cells[df_column] #selects the analysis column
                            df_meta_data_1 = pd.concat([df_meta_data_1, phenotype_df], axis=1, ignore_index=True) #adds the file data frame into the analysis dataframe

                    # checks whether the file is in the metadata list
                    if filename in meta_data_list_2:
                        print('Running phenotype statistics for: ' + str(filedir))
                        Cells = pd.read_csv(filedir) #imports the file
                        Cells['Meta_data_category'] = '2' #adds the metadata category into the data frame

                        df_meta_data_2_column_list.append(filename_string) #appends the file name

                        if select_neighbors[0] == True and percent_of[0] == 'lineage':

                            Selected_Cells = pd.DataFrame().reindex_like(Cells)  # data frame of the cells selected for a specific neighbor
                            Selected_Cells.drop(Selected_Cells.index, inplace=True)

                            # goes throw each cell
                            for index, cell in Cells.iterrows():

                                cell_neigboorhood = ast.literal_eval(cell['Neigboorhood']) # imports the neighborhood
                                Cells_n = Cells[Cells['ImageNumber'] == cell['ImageNumber']] # selects the cells in the same image

                                # goes throw each cell in the image
                                for index, cell_n in Cells_n.iterrows():

                                    # selects the neighboring cells
                                    if int(cell_n['Cell_number']) in cell_neigboorhood:

                                        cell_n_types = ast.literal_eval(cell_n['cell_types_and_states']) # imports the neighboring cell type

                                        # adds the cells with a specific neighbor to the new data frame
                                        if add_mixed_cells == False:
                                            if len(cell_n_types) == 2 and select_neighbors[1] == cell_n_types[0][0]:
                                                Selected_Cells.loc[len(Selected_Cells)] = cell

                                        elif add_mixed_cells == True:
                                            for cell_n_type in cell_n_types:
                                                if select_neighbors[1] == cell_n_type[0]:
                                                    Selected_Cells.loc[len(Selected_Cells)] = cell

                            phenotype_df = Selected_Cells[df_column]  # selects the analysis column
                            df_meta_data_2 = pd.concat([df_meta_data_2, phenotype_df], axis=1, ignore_index=True)  # adds the file data frame into the analysis dataframe
                            
                        else:
                            phenotype_df = Cells[df_column] #selects the analysis column
                            df_meta_data_2 = pd.concat([df_meta_data_2, phenotype_df], axis=1, ignore_index=True) #adds the file data frame into the analysis dataframe

    df_meta_data_1.drop(df_meta_data_1.columns[[0]], axis=1, inplace=True) #drops the first column used to initialize the data frame
    df_meta_data_1.columns = df_meta_data_1_column_list #adds the file names as columns

    df_meta_data_2.drop(df_meta_data_2.columns[[0]], axis=1, inplace=True) #drops the first column used to initialize the data frame
    df_meta_data_2.columns = df_meta_data_2_column_list #adds the file names

    # takes all cells as a reference point and adds the counts of cells with just one cell type
    if percent_of[0] == 'All' and add_mixed_cells == False:

        #creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount','TotalIdentifiedTissueCellCount', 'TotalCellCount','FileName','MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns = cell_types) #creates a dataframe out the list

        #goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            #goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                #counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                #goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    #adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    #adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            #appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        #goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            #resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = [] # empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            #runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            #adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot*0.8, str(round(p,5)), fontsize=12, color='black')
            axs.set_ylim(min_plot-((max_plot-min_plot)*0.1), max_plot+((max_plot-min_plot)*0.1))
            axs.set_ylabel("Percent of all cells")
            axs.set_title(str(cell_type))

            #creates the folder
            if os.path.isdir(output_folder+ '/' + 'Percent of all cells') == True:
                pass
            else:
                os.makedirs(output_folder+ '/' + 'Percent of all cells')

            #saves the figure
            export_file_path = output_folder + '/' + 'Percent of all cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            #exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of all cells' + '/' + 'percent_of_all_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes all cells as a reference point and adds the counts of mixed cells
    if percent_of[0] == 'All' and add_mixed_cells == True:

        # creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount', 'TotalIdentifiedTissueCellCount', 'TotalCellCount', 'FileName', 'MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=cell_types) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = []# empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of all cells")
            axs.set_title(str(cell_type))

            # creates the folder
            if os.path.isdir(output_folder + '/' + 'Percent of all cells') == True:
                pass
            else:
                os.makedirs(output_folder + '/' + 'Percent of all cells')

            # saves the figure
            export_file_path = output_folder + '/' + 'Percent of all cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            # exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of all cells' + '/' + 'percent_of_all_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes all immune cells as a reference point and adds the counts of cells with just one immune cell type
    if percent_of[0] == 'CD45+' and add_mixed_cells == False:

        # creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount', 'TotalIdentifiedTissueCellCount', 'TotalCellCount', 'FileName', 'MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=cell_types) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = []  # empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalIdentifiedImmuneCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of immune cells")
            axs.set_title(str(cell_type))

            # creates the folder
            if os.path.isdir(output_folder + '/' + 'Percent of immune cells') == True:
                pass
            else:
                os.makedirs(output_folder + '/' + 'Percent of immune cells')

            # saves the figure
            export_file_path = output_folder + '/' + 'Percent of immune cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            # exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of immune cells' + '/' + 'percent_of_immune_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes all immune cells as a reference point and adds the counts of mixed immune cells
    if percent_of[0] == 'CD45+' and add_mixed_cells == True:

        # creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount', 'TotalIdentifiedTissueCellCount', 'TotalCellCount', 'FileName', 'MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=cell_types) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = [] # empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalIdentifiedImmuneCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of immune cells")
            axs.set_title(str(cell_type))

            # creates the folder
            if os.path.isdir(output_folder + '/' + 'Percent of immune cells') == True:
                pass
            else:
                os.makedirs(output_folder + '/' + 'Percent of immune cells')

            # saves the figure
            export_file_path = output_folder + '/' + 'Percent of immune cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            # exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of immune cells' + '/' + 'percent_of_immune_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes all tissue cells as a reference point and adds the counts of cells with just one tissue cell type
    if percent_of[0] == 'CD45-' and add_mixed_cells == False:

        # creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount', 'TotalIdentifiedTissueCellCount', 'TotalCellCount', 'FileName', 'MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=cell_types) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = [] # empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalIdentifiedTissueCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of tissue cells")
            axs.set_title(str(cell_type))

            # creates the folder
            if os.path.isdir(output_folder + '/' + 'Percent of tissue cells') == True:
                pass
            else:
                os.makedirs(output_folder + '/' + 'Percent of tissue cells')

            # saves the figure
            export_file_path = output_folder + '/' + 'Percent of tissue cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            # exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of tissue cells' + '/' + 'percent_of_tissue_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes all tissue cells as a reference point and adds the counts of mixed tissue cells
    if percent_of[0] == 'CD45-' and add_mixed_cells == True:

        # creates a datframe for the analysis
        analysis_values = ['TotalIdentifiedImmuneCellCount', 'TotalIdentifiedTissueCellCount', 'TotalCellCount', 'FileName', 'MetaData'] #adds the needed data
        cell_types.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=cell_types) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell type
            for cell_type in cell_types[:-len(analysis_values)]:

                # counters
                row_count = 0
                cell_type_count = 0
                total_identified_immune_cell_count = 0
                total_identified_tissue_cell_count = 0

                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    row_count = row_count + 1
                    cell_type_data = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the cell type count
                    if len(cell_type_data) == 2 and cell_type_data[0][0] == cell_type:
                        cell_type_count = cell_type_count + 1
                    if len(cell_type_data) > 2:
                        for mixed_cell in cell_type_data:
                            if mixed_cell[0] == cell_type:
                                cell_type_count = cell_type_count + 1

                    # adds a counter to the total immune cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45+':
                        total_identified_immune_cell_count = total_identified_immune_cell_count + 1

                    # adds a counter to the total tissue cell count
                    if len(cell_type_data) >= 2 and cell_type_data[-1][0] == 'CD45-':
                        total_identified_tissue_cell_count = total_identified_tissue_cell_count + 1

                analysis_row.append(cell_type_count) #append the cell type counter

            # appends the counters
            analysis_row.append(total_identified_immune_cell_count)
            analysis_row.append(total_identified_tissue_cell_count)
            analysis_row.append(row_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_type in cell_types[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = []# empty list for the plot
            for key, value in DataFrameDict_analysis.items():

                cell_type_list = DataFrameDict_analysis[key][cell_type].to_list() #creates a list of the celltypes
                total_cell_list = DataFrameDict_analysis[key]['TotalIdentifiedTissueCellCount'].to_list() #creates a list of the total cell count
                percent_cells = [int(ct) / int(c) for ct, c in zip(cell_type_list, total_cell_list)] #calculates the cell percentage
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of tissue cells")
            axs.set_title(str(cell_type))

            # creates the folder
            if os.path.isdir(output_folder + '/' + 'Percent of tissue cells') == True:
                pass
            else:
                os.makedirs(output_folder + '/' + 'Percent of tissue cells')

            # saves the figure
            export_file_path = output_folder + '/' + 'Percent of tissue cells' + '/' + str(cell_type)
            fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

            # exports the final analysis data frame
            export_file_path = output_folder + '/' + 'Percent of tissue cells' + '/' + 'percent_of_tissue_cells_phenotype_analysis.csv'
            analysis_df_final.to_csv(export_file_path)

    # takes lineage cells as a reference point and adds the counts of cell subsets
    if percent_of[0] == 'lineage' and add_mixed_cells == False:

        state_list = [] #creates a unique list for all the states

        # goes throw each column (file)
        for column in df_meta_data_1:

            # goes throw each row
            for row in df_meta_data_1[column].dropna():
                cell_types = ast.literal_eval(row) #imports the cell type

                if len(cell_types) == 2 and cell_types[0][0] == percent_of[1]: #selects the correct cell type

                    #gets all unique states and adds them to the list
                    for state in cell_types[0][1:]:
                        state_list.append(state)

        # goes throw each column (file)
        for column in df_meta_data_1:

            # goes throw each row
            for row in df_meta_data_1[column].dropna():
                cell_types = ast.literal_eval(row)

                if len(cell_types) == 2 and cell_types[0][0] == percent_of[1]: #selects the correct cell type

                    # gets all unique states and adds them to the list
                    for state in cell_types[0][1:]:
                        state_list.append(state)

        state_list = list(set(state_list)) #creates a set to remove al multiples
        analysis_values = ['CellTypeCount', 'FileName', 'MetaData'] #adds the needed data
        state_list.extend(analysis_values) #adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=state_list) #creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell state
            for cell_state in state_list[:-len(analysis_values)]:

                # counters
                cell_type_count = 0
                cell_state_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    cell_types = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the appropriate cell state
                    if len(cell_types) == 2 and cell_types[0][0] == percent_of[1]:
                        cell_type_count = cell_type_count + 1

                        for state in cell_types[0][1:]:

                            if state == cell_state:
                                cell_state_count = cell_state_count + 1

                analysis_row.append(cell_state_count) #append the cell type counter

            # appends the counters
            analysis_row.append(cell_type_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell state
            for cell_state in state_list[:-len(analysis_values)]:

                # counters
                cell_type_count = 0
                cell_state_count = 0
                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    cell_types = ast.literal_eval(row) #imports the cell type

                    # adds a counter to the appropriate cell state
                    if len(cell_types) == 2 and cell_types[0][0] == percent_of[1]:
                        cell_type_count = cell_type_count + 1

                        for state in cell_types[0][1:]:

                            if state == cell_state:
                                cell_state_count = cell_state_count + 1

                analysis_row.append(cell_state_count) #append the cell type counter

            # appends the counters
            analysis_row.append(cell_type_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row #appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_state in state_list[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = []# empty list for the plot
            for key, value in DataFrameDict_analysis.items():
                cell_type_list = DataFrameDict_analysis[key][cell_state].to_list() #creates a list of the cell staes
                total_cell_list = DataFrameDict_analysis[key]['CellTypeCount'].to_list() #creates a list of the cell count

                # calculates the cell state percentage
                percent_cells = []
                for ct, c in zip(cell_type_list, total_cell_list):
                    if int(c) == 0:
                        percent_cells.append(0)
                    else:
                        percent_cells.append(int(ct) / int(c))
                
                percent_cells_final = [x * 100 for x in percent_cells] #multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) #initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of " + str(percent_of[1]) + " cells")
            axs.set_title(str(cell_state))

            if select_neighbors[0] == False:

                # creates the folder
                if os.path.isdir(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells') == True:
                    pass
                else:
                    os.makedirs(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells')

                # saves the figure
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells' + '/' + 'Percent of ' + str(percent_of[1]) + ' cells ' + str(cell_state)
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

                # exports the final analysis data frame
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells' + '/' + 'percent_of_' + str(percent_of[1]) + '_cells_phenotype_analysis.csv'
                analysis_df_final.to_csv(export_file_path)

            if select_neighbors[0] == True:

                # creates the folder
                if os.path.isdir(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells') == True:
                    pass
                else:
                    os.makedirs(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells')

                # saves the figure
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells' + '/' + 'Percent of ' + str(percent_of[1]) + ' cells ' + str(cell_state) + ' neighbouring ' + str(select_neighbors[1]) + ' cells'
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

                # exports the final analysis data frame
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells' + '/' + 'percent_of_' + str(percent_of[1]) + '_cells_neighbouring_' + str(select_neighbors[1]) + '_cells_phenotype_analysis.csv'
                analysis_df_final.to_csv(export_file_path)

    # takes mixed lineage cells as a reference point and adds the counts of cell subsets
    if percent_of[0] == 'lineage' and add_mixed_cells == True:

        state_list = [] # creates a unique list for all the states

        # goes throw each column (file)
        for column in df_meta_data_1:

            # goes throw each row
            for row in df_meta_data_1[column].dropna():
                cell_types = ast.literal_eval(row)

                for cell_type in cell_types: # imports the cell type

                    if cell_type[0] == percent_of[1]:  # selects the correct cell type

                        # gets all unique states and adds them to the list
                        for state in cell_types[0][1:]:
                            state_list.append(state)

        # goes throw each column (file)
        for column in df_meta_data_1:

            # goes throw each row
            for row in df_meta_data_1[column].dropna():
                cell_types = ast.literal_eval(row)

                for cell_type in cell_types:

                    if cell_type[0] == percent_of[1]: # selects the correct cell type

                        # gets all unique states and adds them to the list
                        for state in cell_types[0][1:]:
                            state_list.append(state)

        state_list = list(set(state_list)) # creates a set to remove al multiples
        analysis_values = ['CellTypeCount', 'FileName', 'MetaData'] # adds the needed data
        state_list.extend(analysis_values) # adds the analysis to the wanted cell types
        analysis_df = pd.DataFrame(columns=state_list) # creates a dataframe out the list

        # goes throw each column (file)
        for column in df_meta_data_1:

            analysis_row = []
            # goes throw each cell state
            for cell_state in state_list[:-len(analysis_values)]:

                # counters
                cell_type_count = 0
                cell_state_count = 0

                # goes throw each row
                for row in df_meta_data_1[column].dropna():

                    cell_types = ast.literal_eval(row) # imports the cell type
                    
                    for cell_type in cell_types:

                        # adds a counter to the appropriate cell state
                        if cell_type[0] == percent_of[1]:
                            cell_type_count = cell_type_count + 1
    
                            for state in cell_type[1:]:
    
                                if state == cell_state:
                                    cell_state_count = cell_state_count + 1

                analysis_row.append(cell_state_count) # append the cell type counter

            # appends the counters
            analysis_row.append(cell_type_count)
            analysis_row.append(column)
            analysis_row.append(1)

            analysis_df.loc[len(analysis_df)] = analysis_row # appends the final row into the analysis df

        # goes throw each column (file)
        for column in df_meta_data_2:

            analysis_row = []
            # goes throw each cell state
            for cell_state in state_list[:-len(analysis_values)]:

                # counters
                cell_type_count = 0
                cell_state_count = 0
                # goes throw each row
                for row in df_meta_data_2[column].dropna():

                    cell_types = ast.literal_eval(row) # imports the cell type

                    # adds a counter to the appropriate cell state
                    for cell_type in cell_types:

                        if cell_type[0] == percent_of[1]:
                            cell_type_count = cell_type_count + 1

                            for state in cell_type[1:]:

                                if state == cell_state:
                                    cell_state_count = cell_state_count + 1

                analysis_row.append(cell_state_count) # append the cell type counter

            # appends the counters
            analysis_row.append(cell_type_count)
            analysis_row.append(column)
            analysis_row.append(2)

            analysis_df.loc[len(analysis_df)] = analysis_row # appends the final row into the analysis df

        analysis_df_final = analysis_df

        # goes throw each cell type
        for cell_state in state_list[:-len(analysis_values)]:

            DataFrameDict_analysis = {elem: pd.DataFrame() for elem in analysis_df_final.MetaData.unique()} # splits the analysis data frame by metadata

            # resets the index for the data frame dictionary
            for key in DataFrameDict_analysis.keys():
                DataFrameDict_analysis[key] = analysis_df_final[:][analysis_df_final.MetaData == key].reset_index()

            plot_list = [] # empty list for the plot
            for key, value in DataFrameDict_analysis.items():
                cell_type_list = DataFrameDict_analysis[key][cell_state].to_list() # creates a list of the cell staes
                total_cell_list = DataFrameDict_analysis[key]['CellTypeCount'].to_list() # creates a list of the cell count

                # calculates the cell state percentage
                percent_cells = []
                for ct, c in zip(cell_type_list, total_cell_list):
                    if int(c) == 0:
                        percent_cells.append(0)
                    else:
                        percent_cells.append(int(ct) / int(c))

                percent_cells_final = [x * 100 for x in percent_cells] # multiplies each value by 100 for the final percentage
                plot_list.append(percent_cells_final)

            # runs a mannwhitneyu test on the different metadata catergories
            U1, p = mannwhitneyu(plot_list[0], plot_list[1], method='auto')
            max_plot = max([max(plot) for plot in plot_list])
            min_plot = min([min(plot) for plot in plot_list])

            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 4)) # initilizes the figure

            # adds the data into the image
            axs.scatter([1 for x in range(len(plot_list[0]))], plot_list[0], c='#D01515', s=15)
            axs.scatter([2 for x in range(len(plot_list[1]))], plot_list[1], c='#142EFB', s=15)
            axs.boxplot(plot_list)
            axs.text(1.3, max_plot * 0.8, str(round(p, 5)), fontsize=12, color='black')
            axs.set_ylim(min_plot - ((max_plot - min_plot) * 0.1), max_plot + ((max_plot - min_plot) * 0.1))
            axs.set_ylabel("Percent of " + str(percent_of[1]) + " cells")
            axs.set_title(str(cell_state))

            if select_neighbors[0] == False:

                # creates the folder
                if os.path.isdir(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells') == True:
                    pass
                else:
                    os.makedirs(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells')

                # saves the figure
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells' + '/' + 'Percent of ' + str(percent_of[1]) + ' cells ' + str(cell_state)
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

                # exports the final analysis data frame
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells' + '/' + 'percent_of_' + str(percent_of[1]) + '_cells_phenotype_analysis.csv'
                analysis_df_final.to_csv(export_file_path)

            if select_neighbors[0] == True:

                # creates the folder
                if os.path.isdir(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells') == True:
                    pass
                else:
                    os.makedirs(output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells')

                # saves the figure
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells' + '/' + 'Percent of ' + str(percent_of[1]) + ' cells ' + str(cell_state) + ' neighbouring ' + str(select_neighbors[1]) + ' cells'
                fig.savefig(export_file_path, bbox_inches='tight', dpi=600)

                # exports the final analysis data frame
                export_file_path = output_folder + '/' + 'Percent of ' + str(percent_of[1]) + ' cells neighbouring ' + str(select_neighbors[1]) + ' cells' + '/' + 'percent_of_' + str(percent_of[1]) + '_cells_neighbouring_' + str(select_neighbors[1]) + '_cells_phenotype_analysis.csv'
                analysis_df_final.to_csv(export_file_path)