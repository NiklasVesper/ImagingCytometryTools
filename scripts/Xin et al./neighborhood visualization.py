import pandas as pd
import re
import statistics
import matplotlib.pyplot as plt
import scandir as sd
import os


folder_dir = r'F:\ATF6 FAK\All\High\CD11'

#for multiple files

for paths, dirs, files in sd.walk(folder_dir): #goes throw all files and folders in given directory

    neighborhood_Tissue_percent = []
    neighborhood_Immune_percent = []

    neighborhood_CD4_percent = []
    neighborhood_CD8_percent = []
    neighborhood_CD20_percent = []
    neighborhood_CD15_percent = []
    neighborhood_CD68_percent = []
    neighborhood_DC_percent = []

    for file in os.listdir(paths): #goes throw all files in a folder
        filedir = os.path.join(paths, file) #returns fuull file directory

        if filedir.endswith("CD11c_neighborhood.csv"):
            print(filedir)

            filename = os.path.basename(file)
            filename_string = str(filename)
            filedir_string = str(filedir)[:-len(filename)]

            neighborhood_of_cells = pd.read_csv(filedir)


            for index, cell in neighborhood_of_cells.iterrows():


                if cell['Cell_types'] == "['Immune cell', 'DC cell']":


                    Neigboorhood = re.split(r'], ', cell['CD11 neigboors'])

                    Cell_list = []
                    for x in Neigboorhood:
                        Cell_list.append(x.replace('[', '').replace(']', ''))

                    number_of_cells = len(Cell_list)

                    Tissue_count = []

                    Immune_count = []
                    CD4_count = []
                    CD8_count = []
                    B_count = []
                    Granulocyte_count = []
                    Myeloid_count = []
                    DC_count = []

                    for y in Cell_list:

                        if y == "'Tissue cell'":
                            Tissue_count.append('1 count')

                        # ----------

                        if y == "'Immune cell'":
                            Immune_count.append('1 count')

                        if y == "'Immune cell', 'CD4 T cell'":
                            Immune_count.append('1 count')
                            CD4_count.append('1 count')

                        if y == "'Immune cell', 'CD8 T cell'":
                            Immune_count.append('1 count')
                            CD8_count.append('1 count')

                        # ----------

                        if y == "'Immune cell', 'B cell'":
                            Immune_count.append('1 count')
                            B_count.append('1 count')

                        if y == "'Immune cell', 'Granulocyte'":
                            Immune_count.append('1 count')
                            Granulocyte_count.append('1 count')

                        if y == "'Immune cell', 'Myeloid cell'":
                            Immune_count.append('1 count')
                            Myeloid_count.append('1 count')

                        if y == "'Immune cell', 'DC cell'":
                            Immune_count.append('1 count')
                            DC_count.append('1 count')


                    Tissue_percent = (len(Tissue_count) / len(Cell_list)) * 100
                    Immune_percent = (len(Immune_count) / len(Cell_list)) * 100

                    neighborhood_Tissue_percent.append(Tissue_percent)
                    neighborhood_Immune_percent.append(Immune_percent)

                    CD4_percent = (len(CD4_count) / len(Cell_list)) * 100
                    CD8_percent = (len(CD8_count) / len(Cell_list)) * 100
                    B_percent = (len(B_count) / len(Cell_list)) * 100
                    Granulocyte_percent = (len(Granulocyte_count) / len(Cell_list)) * 100
                    Myeloid_percent = (len(Myeloid_count) / len(Cell_list)) * 100
                    DC_percent = (len(DC_count) / len(Cell_list)) * 100

                    neighborhood_CD4_percent.append(CD4_percent)
                    neighborhood_CD8_percent.append(CD8_percent)
                    neighborhood_CD20_percent.append(B_percent)
                    neighborhood_CD15_percent.append(Granulocyte_percent)
                    neighborhood_CD68_percent.append(Myeloid_percent)
                    neighborhood_DC_percent.append(DC_percent)

            #--------

    try:

        categories = ['Tissue', 'Immune']

        values = [statistics.mean(neighborhood_Tissue_percent),statistics.mean(neighborhood_Immune_percent)]

        errors = [statistics.stdev(neighborhood_Tissue_percent),statistics.stdev(neighborhood_Immune_percent)]

        colors = ['gray', '#44933a']

        # Create the bar plot
        fig, ax = plt.subplots()

        # Plot bars with error bars
        ax.bar(categories, values, yerr=errors, capsize=5, color=colors, edgecolor='black')

        # Add labels and title
        ax.set_ylabel('Percent of Neighboring cells')
        ax.set_title('')

        # ax.set_xticklabels(categories, rotation=30)
        ax.set_ylim(0, 110)

        # Show the plot
        plt.show()
        fig.savefig(folder_dir + "/CD11c_neighborhood_vis",dpi=300)
        plt.close()

        categories = ['CD4', 'CD8', 'CD20', 'CD15', 'CD68','CD11c']

        values = [statistics.mean(neighborhood_CD4_percent),
                  statistics.mean(neighborhood_CD8_percent),
                  statistics.mean(neighborhood_CD20_percent),
                  statistics.mean(neighborhood_CD15_percent),
                  statistics.mean(neighborhood_CD68_percent),
                  statistics.mean(neighborhood_DC_percent)
                  ]

        errors = [statistics.stdev(neighborhood_CD4_percent),
                  statistics.stdev(neighborhood_CD8_percent),
                  statistics.stdev(neighborhood_CD20_percent),
                  statistics.stdev(neighborhood_CD15_percent),
                  statistics.stdev(neighborhood_CD68_percent),
                  statistics.stdev(neighborhood_DC_percent),
                  ]

        colors = ['#E40303', '#FF8C00', '#FFED00', '#008026', '#24408E','#732982']

        # Create the bar plot
        fig, ax = plt.subplots()

        # Plot bars with error bars
        ax.bar(categories, values, yerr=errors, capsize=5, color=colors, edgecolor='black')

        # Add labels and title

        ax.set_ylabel('Percent of Neighboring cells')
        ax.set_title('')
        ax.set_ylim(0, 60)

        plt.show()
        fig.savefig(folder_dir + "/CD11c_immune_neighborhood_vis", dpi=300)
        plt.close()

        neigboorhood_df = pd.DataFrame()

        neigboorhood_df['CD4'] = neighborhood_CD4_percent
        neigboorhood_df['CD8'] = neighborhood_CD8_percent
        neigboorhood_df['CD20'] = neighborhood_CD20_percent
        neigboorhood_df['CD15'] = neighborhood_CD15_percent
        neigboorhood_df['CD68'] = neighborhood_CD68_percent
        neigboorhood_df['CD11c'] = neighborhood_DC_percent

        neigboorhood_df.to_csv(folder_dir + "/" + "CD11c_immune_neighborhood.csv")


    except statistics.StatisticsError:
        print('------')
        print(filedir)
        print('This file hase nothing to calculate')
        print('------')
        continue
