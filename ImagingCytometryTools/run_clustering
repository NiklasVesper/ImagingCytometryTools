import pandas as pd
import numpy as np
import sklearn.cluster
import phenograph
import scandir as sd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

from ImagingCytometryTools.get_data_from_files import get_markers_from_csv

def run_clustering(directory, filestring, new_folder_name, output_folder_analysis, clustering_mode, clustering_channels = 'all', heat_map_channels='all', show_mean = True, log_scale = True, split_by_cluster = False, save_data = True):

    file_list = []
    file_names = []
    for paths, dirs, files in sd.walk(directory): #goes throw all files and folders in given directory

        for file in os.listdir(paths): #goes throw all files in a folder
            filedir = os.path.join(paths, file) #returns full file directory

            if filedir.endswith(".csv"): #returns all files that end with .csv

                #gives you the file name
                filename = os.path.basename(file)
                filename_string = str(filename)

                #gives you the directory name
                dir_name = os.path.dirname(filedir)
                dir_name_string = str(dir_name)

                #gives you the folder name
                folder_name = os.path.basename(dir_name)
                folder_name_string = str(folder_name)

                if filestring in filename_string: #returns all files that contain the provided file string
                    file_list.append(filedir)
                    file_names.append(filename_string)

    #creates the folder for the merged files and the heatmap
    if save_data == True:
        if os.path.isdir(output_folder_analysis) == True:
            pass
        else:
            os.makedirs(output_folder_analysis)

    #merges all the found files into one data frame
    files_merged = pd.read_csv(file_list[0])
    for file in file_list[1:]:
        files_merged = pd.concat([files_merged, pd.read_csv(file)])

    #selects the channels for the clustering
    if clustering_channels == 'all':
        files_merged_clustering = files_merged.loc[:, get_markers_from_csv(files_merged)]
    else:
        files_merged_clustering = files_merged.loc[:,clustering_channels] #removes all the channels that are not wanted for the clustering
    if heat_map_channels == 'all':
        heat_map_channels_vis = files_merged.loc[:, get_markers_from_csv(files_merged)]
    else:
        heat_map_channels_vis = files_merged.loc[:,heat_map_channels]  #removes all the channels that are not wanted for the clustering
     
    data_matrix = files_merged_clustering.values #creates a data matrix for the clustering

    if clustering_mode[0] == 'KMeans':

        print('Running KMeans clustering for: ' + str(file_names))

        K_means_clustering = sklearn.cluster.KMeans(n_clusters = clustering_mode[1]).fit(data_matrix) #runs the clustering
        lables = K_means_clustering.labels_ #extracts the cluster labels
        files_merged['KMeans_clustering'] = lables #appends all the found clusters into the data frame

        #saves the clustering per file location
        if save_data == True:

            #splits the merged file by individual file directories
            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index() #creates the data frame to be saved

                output_folder = str(key)[:-3] + new_folder_name #creates the output folder

                #checks if the output folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                #creates the file name directory
                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') +'.csv'

                #checks if the file name directory allready exists
                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path) #safes the data frame

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv' #creates the file name directory for the merged file

            #checks if the file name directory allready exists
            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        #finds all the unique clusters and creates a dictionary out of it
        Unique_clusters = files_merged.KMeans_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        #gets the means for each individual cluster
        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.KMeans_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame() #creates an empty dataframe for the cluster means
        for chanel in heat_map_channels_vis:
            columb = [] #empty list for the mean intensities
            key_list = [] #list for cluster numbers

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean()) #gets the mean by cluster
                key_list.append(key)

            heat_map_df[chanel] = columb #appends the mean

        heat_map_df['cluster'] = key_list #appends the cluster number

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv' #creates a export file directory for the heatmap data frame

            # checks if the new file directory already exists
            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1] #sorts the data frame by cluster number

        #uses a log scale foir the heatmap
        if log_scale == True:
            sns.heatmap(heat_map_df_image, annot = show_mean, norm = LogNorm(), cmap="coolwarm") #shows means

        #uses a linear scale foir the heatmap
        if log_scale == False:
            sns.heatmap(heat_map_df_image, annot=show_mean, cmap="coolwarm")  # shows means

        #visualizes the heatmap
        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        elif save_data == False:
            plt.show()

    elif clustering_mode[0] == 'AgglomerativeClustering':

        print('Running Agglomerative clustering for:' + str(file_names))

        Agglomerative_Clustering = sklearn.cluster.AgglomerativeClustering(n_clusters = clustering_mode[1]).fit(data_matrix)
        lables = Agglomerative_Clustering.labels_

        files_merged['Agglomerative_clustering'] = lables

        # saves the clustering per file location
        if save_data == True:

            # splits the merged file by individual file directories
            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()  # creates the data frame to be saved

                output_folder = str(key)[:-3] + new_folder_name  # creates the output folder

                # checks if the output folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                # creates the file name directory
                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                # checks if the file name directory allready exists
                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)  # safes the data frame

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'  # creates the file name directory for the merged file

            # checks if the file name directory allready exists
            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        # finds all the unique clusters and creates a dictionary out of it
        Unique_clusters = files_merged.Agglomerative_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        # gets the means for each individual cluster
        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.Agglomerative_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()  # creates an empty dataframe for the cluster means
        for chanel in heat_map_channels_vis:
            columb = []  # empty list for the mean intensities
            key_list = []  # list for cluster numbers

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())  # gets the mean by cluster
                key_list.append(key)

            heat_map_df[chanel] = columb  # appends the mean

        heat_map_df['cluster'] = key_list  # appends the cluster number

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'  # creates a export file directory for the heatmap data frame

            # checks if the new file directory already exists
            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]  # sorts the data frame by cluster number

        # uses a log scale foir the heatmap
        if log_scale == True:
            sns.heatmap(heat_map_df_image, annot=show_mean, norm=LogNorm(), cmap="coolwarm")  # shows means

        # uses a linear scale foir the heatmap
        if log_scale == False:
            sns.heatmap(heat_map_df_image, annot=show_mean, cmap="coolwarm")  # shows means

        # visualizes the heatmap
        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        elif save_data == False:
            plt.show()

    elif clustering_mode[0] == 'OPTICS':

        print('Running OPTICS clustering for:' + str(file_names))

        OPTICS_clustering = sklearn.cluster.OPTICS(min_samples = clustering_mode[1]).fit(data_matrix)
        lables = OPTICS_clustering.labels_

        files_merged['OPTICS_clustering'] = lables
        # saves the clustering per file location
        if save_data == True:

            # splits the merged file by individual file directories
            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()  # creates the data frame to be saved

                output_folder = str(key)[:-3] + new_folder_name  # creates the output folder

                # checks if the output folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                # creates the file name directory
                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                # checks if the file name directory allready exists
                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)  # safes the data frame

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'  # creates the file name directory for the merged file

            # checks if the file name directory allready exists
            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        # finds all the unique clusters and creates a dictionary out of it
        Unique_clusters = files_merged.OPTICS_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        # gets the means for each individual cluster
        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.OPTICS_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()  # creates an empty dataframe for the cluster means
        for chanel in heat_map_channels_vis:
            columb = []  # empty list for the mean intensities
            key_list = []  # list for cluster numbers

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())  # gets the mean by cluster
                key_list.append(key)

            heat_map_df[chanel] = columb  # appends the mean

        heat_map_df['cluster'] = key_list  # appends the cluster number

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'  # creates a export file directory for the heatmap data frame

            # checks if the new file directory already exists
            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]  # sorts the data frame by cluster number

        # uses a log scale foir the heatmap
        if log_scale == True:
            sns.heatmap(heat_map_df_image, annot=show_mean, norm=LogNorm(), cmap="coolwarm")  # shows means

        # uses a linear scale foir the heatmap
        if log_scale == False:
            sns.heatmap(heat_map_df_image, annot=show_mean, cmap="coolwarm")  # shows means

        # visualizes the heatmap
        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        elif save_data == False:
            plt.show()

    elif clustering_mode[0] == 'BisectingKMeans':

        print('Running BisectingKMeans clustering for:' + str(file_names))

        Bisecting_K_means_clustering = sklearn.cluster.BisectingKMeans(n_clusters = clustering_mode[1]).fit(data_matrix)
        lables = Bisecting_K_means_clustering.labels_

        files_merged['BisectingKMeans_clustering'] = lables
        # saves the clustering per file location
        if save_data == True:

            # splits the merged file by individual file directories
            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()  # creates the data frame to be saved

                output_folder = str(key)[:-3] + new_folder_name  # creates the output folder

                # checks if the output folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                # creates the file name directory
                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                # checks if the file name directory allready exists
                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)  # safes the data frame

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'  # creates the file name directory for the merged file

            # checks if the file name directory allready exists
            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        # finds all the unique clusters and creates a dictionary out of it
        Unique_clusters = files_merged.BisectingKMeans_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        # gets the means for each individual cluster
        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.BisectingKMeans_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()  # creates an empty dataframe for the cluster means
        for chanel in heat_map_channels_vis:
            columb = []  # empty list for the mean intensities
            key_list = []  # list for cluster numbers

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())  # gets the mean by cluster
                key_list.append(key)

            heat_map_df[chanel] = columb  # appends the mean

        heat_map_df['cluster'] = key_list  # appends the cluster number

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'  # creates a export file directory for the heatmap data frame

            # checks if the new file directory already exists
            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]  # sorts the data frame by cluster number

        # uses a log scale foir the heatmap
        if log_scale == True:
            sns.heatmap(heat_map_df_image, annot=show_mean, norm=LogNorm(), cmap="coolwarm")  # shows means

        # uses a linear scale foir the heatmap
        if log_scale == False:
            sns.heatmap(heat_map_df_image, annot=show_mean, cmap="coolwarm")  # shows means

        # visualizes the heatmap
        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        elif save_data == False:
            plt.show()

    elif clustering_mode[0] == 'PhenoGraph':

        print('PhenoGraph clustering for:' + str(file_names))

        communities = phenograph.cluster(data_matrix, k = clustering_mode[1])
        lables = communities

        files_merged['PhenoGraph_clustering'] = lables[0]
        # saves the clustering per file location
        if save_data == True:

            # splits the merged file by individual file directories
            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()  # creates the data frame to be saved

                output_folder = str(key)[:-3] + new_folder_name  # creates the output folder

                # checks if the output folder exists
                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                # creates the file name directory
                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                # checks if the file name directory allready exists
                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)  # safes the data frame

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'  # creates the file name directory for the merged file

            # checks if the file name directory allready exists
            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        # finds all the unique clusters and creates a dictionary out of it
        Unique_clusters = files_merged.PhenoGraph_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        # gets the means for each individual cluster
        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.PhenoGraph_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()  # creates an empty dataframe for the cluster means
        for chanel in heat_map_channels_vis:
            columb = []  # empty list for the mean intensities
            key_list = []  # list for cluster numbers

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())  # gets the mean by cluster
                key_list.append(key)

            heat_map_df[chanel] = columb  # appends the mean

        heat_map_df['cluster'] = key_list  # appends the cluster number

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'  # creates a export file directory for the heatmap data frame

            # checks if the new file directory already exists
            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]  # sorts the data frame by cluster number

        # uses a log scale foir the heatmap
        if log_scale == True:
            sns.heatmap(heat_map_df_image, annot=show_mean, norm=LogNorm(), cmap="coolwarm")  # shows means

        # uses a linear scale foir the heatmap
        if log_scale == False:
            sns.heatmap(heat_map_df_image, annot=show_mean, cmap="coolwarm")  # shows means

        # visualizes the heatmap
        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        elif save_data == False:
            plt.show()

    if cuustering_mode == 'FlowSOM':
        print('to do')
        #https://github.com/saeyslab/FlowSOM_Python
        

    else:
        print('This clustering algorithm is not available. You can use: KMeans, AgglomerativeClustering, OPTICS, BisectingKMeans and PhenoGraph')

    if save_data == True:
        print('I phenotyped all the cells. -------------------------------------------------------------------------------')
