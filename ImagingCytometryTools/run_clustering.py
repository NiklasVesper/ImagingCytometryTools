import pandas as pd
import numpy as np
import sklearn.cluster
import phenograph
import scandir as sd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm, Normalize

from ImagingCytometryTools.get_data_from_files import get_markers_from_csv

'''
https://scikit-learn.org/stable/modules/clustering.html
https://scikit-learn.org/stable/api/sklearn.cluster.html#module-sklearn.cluster
https://github.com/dpeerlab/phenograph
'''

def run_clustering(directory, filestring, new_folder_name, output_folder_analysis, clustering_mode, clustering_channels='all', heat_map_channels='all', show_mean=True, log_scale=True, split_by_cluster=False, save_data=True):

    file_list = []
    file_names = []

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

                if filestring in filename_string:

                    file_list.append(filedir)
                    file_names.append(filename_string)

    if save_data == True:
        if os.path.isdir(output_folder_analysis) == True:
            pass
        else:
            os.makedirs(output_folder_analysis)

    files_merged = pd.read_csv(file_list[0])
    for file in file_list[1:]:
        files_merged = pd.concat([files_merged, pd.read_csv(file)])

        print('Running clustering for: ' + str(file))

    if clustering_channels == 'all':
        files_merged_clustering = files_merged.loc[:, get_markers_from_csv(files_merged)]
    else:
        files_merged_clustering = files_merged.loc[:, clustering_channels]

    if heat_map_channels == 'all':
        heat_map_channels_vis = files_merged.loc[:, get_markers_from_csv(files_merged)]
    else:
        heat_map_channels_vis = files_merged.loc[:, heat_map_channels]

    data_matrix = files_merged_clustering.values

    if clustering_mode[0] == 'KMeans':

        K_means_clustering = sklearn.cluster.KMeans(n_clusters=clustering_mode[1]).fit(data_matrix)
        lables = K_means_clustering.labels_
        files_merged['KMeans_clustering'] = lables

        if save_data == True:

            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()

                output_folder = str(key)[:-3] + new_folder_name

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'

            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        Unique_clusters = files_merged.KMeans_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.KMeans_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()
        for chanel in heat_map_channels_vis:
            columb = []
            key_list = []

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())
                key_list.append(key)

            heat_map_df[chanel] = columb

        heat_map_df['cluster'] = key_list

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'

            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]
        heat_map_df_image_norm = (heat_map_df_image - heat_map_df_image.min()) / (heat_map_df_image.max() - heat_map_df_image.min())

        if log_scale == True:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, norm=LogNorm(), cmap="coolwarm")

        if log_scale == False:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, cmap="coolwarm")

        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        if save_data == False:
            plt.show()
            plt.close()

    if clustering_mode[0] == 'OPTICS':

        OPTICS_clustering = sklearn.cluster.OPTICS(min_samples=clustering_mode[1]).fit(data_matrix)
        lables = OPTICS_clustering.labels_

        files_merged['OPTICS_clustering'] = lables

        if save_data == True:

            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()

                output_folder = str(key)[:-3] + new_folder_name

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'

            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        Unique_clusters = files_merged.OPTICS_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.OPTICS_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()
        for chanel in heat_map_channels_vis:
            columb = []
            key_list = []

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())
                key_list.append(key)

            heat_map_df[chanel] = columb

        heat_map_df['cluster'] = key_list

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'

            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]
        heat_map_df_image_norm = (heat_map_df_image - heat_map_df_image.min()) / (heat_map_df_image.max() - heat_map_df_image.min())

        if log_scale == True:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, norm=LogNorm(), cmap="coolwarm")

        if log_scale == False:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, cmap="coolwarm")

        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        if save_data == False:
            plt.show()
            plt.close()

    if clustering_mode[0] == 'PhenoGraph':

        communities = phenograph.cluster(data_matrix, k=clustering_mode[1])
        lables = communities

        files_merged['PhenoGraph_clustering'] = lables[0]

        if save_data == True:

            Unique_Directory = files_merged.FileDirectory.unique()
            DataFrame_Unique_Directory = {elem: pd.DataFrame() for elem in Unique_Directory}

            for key in DataFrame_Unique_Directory.keys():
                DataFrame_Unique_Directory[key] = files_merged[:][files_merged.FileDirectory == key].reset_index()

                output_folder = str(key)[:-3] + new_folder_name

                if os.path.isdir(output_folder) == True:
                    pass
                else:
                    os.makedirs(output_folder)

                first_folder = key[:-len(os.path.basename(key)) - 1]
                second_folder = first_folder[:-len(os.path.basename(first_folder)) - 1]
                export_file_path = str(key)[:-3] + new_folder_name + '/' + str(os.path.basename(second_folder)) + '_' + new_folder_name.replace(' ', '_') + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_Directory[key].to_csv(export_file_path)

            export_file_path_all = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.csv'

            if os.path.isfile(export_file_path_all) == True:
                print('The file: ' + export_file_path_all + ' already exists!')
            else:
                files_merged.to_csv(export_file_path_all)

        Unique_clusters = files_merged.PhenoGraph_clustering.unique()
        DataFrame_Unique_clusters = {elem: pd.DataFrame() for elem in Unique_clusters}

        for key in DataFrame_Unique_clusters.keys():
            DataFrame_Unique_clusters[key] = files_merged[:][files_merged.PhenoGraph_clustering == key].reset_index()

            if save_data == True and split_by_cluster == True:
                export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_' + str(key) + '.csv'

                if os.path.isfile(export_file_path) == True:
                    print('The file: ' + export_file_path + ' already exists!')
                else:
                    DataFrame_Unique_clusters[key].to_csv(export_file_path)

        heat_map_df = pd.DataFrame()
        for chanel in heat_map_channels_vis:
            columb = []
            key_list = []

            for key in DataFrame_Unique_clusters.keys():
                columb.append(DataFrame_Unique_clusters[key][chanel].mean())
                key_list.append(key)

            heat_map_df[chanel] = columb

        heat_map_df['cluster'] = key_list

        if save_data == True:
            export_file_path = output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '_heat_map' + '.csv'

            if os.path.isfile(export_file_path) == True:
                print('The file: ' + export_file_path + ' already exists!')
            else:
                heat_map_df.to_csv(export_file_path)

        heat_map_df_image = heat_map_df.sort_values(by='cluster', ascending=True).reset_index().iloc[:, 1:-1]
        heat_map_df_image_norm = (heat_map_df_image - heat_map_df_image.min()) / (heat_map_df_image.max() - heat_map_df_image.min())

        if log_scale == True:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, norm=LogNorm(), cmap="coolwarm")

        if log_scale == False:
            sns.heatmap(heat_map_df_image_norm, annot=heat_map_df_image, cmap="coolwarm")

        if save_data == True:
            plt.savefig(output_folder_analysis + '/' + new_folder_name.replace(' ', '_') + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        if save_data == False:
            plt.show()
            plt.close()

    if clustering_mode[0] == 'FlowSOM':
        'https://github.com/saeyslab/FlowSOM_Python'

    if save_data == True:
        print('I phenotyped all the cells. -------------------------------------------------------------------------------')