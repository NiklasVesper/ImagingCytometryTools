import pandas as pd
import re
import statistics
import matplotlib.pyplot as plt
import scandir as sd
import os

'''
The "neighborhood visualization" script generates plots to visualize the phenotypes of surrounding cells in relation to one specific cellular phenotype
The script takes files that were previously generated with the "phenotyping" script.

Continue now with the The "neighborhood testing" script.
'''

folder_dir = r'F:\ATF6 FAK\All\High\CD8' #folder directory

#for multiple files

for paths, dirs, files in sd.walk(folder_dir): #goes through all files and folders in given directory

    #lists for the phenotypes
    neighborhood_Tissue_percent = []
    neighborhood_Immune_percent = []

    neighborhood_CD4_percent = []
    neighborhood_CD8_percent = []
    neighborhood_CD20_percent = []
    neighborhood_CD15_percent = []
    neighborhood_CD68_percent = []
    neighborhood_DC_percent = []

    for file in os.listdir(paths): #goes through all files in a folder
        filedir = os.path.join(paths, file) #returns full file directory

        if filedir.endswith("_CD8_neighborhood.csv"): #checks if the file has the proper condition in its name lists for the phenotypes
            '''
            Other conditions are: _CD20_neighborhood.csv / _CD11c_neighborhood.csv 
            '''
            print(filedir)

            filename = os.path.basename(file) #gives you the file name
            filename_string = str(filename) #turns the filename into a string
            filedir_string = str(filedir)[:-len(filename)] #gives you the file directory as a string

            neighborhood_of_cells = pd.read_csv(filedir) #opens the found directory

            for index, cell in neighborhood_of_cells.iterrows(): #loops through all cells

                if cell['Cell_types'] == "['Immune cell', 'CD8 T cell']": #checks for the cell conditions one is interested in

                    '''
                    Other conditions are: "['Immune cell', 'B cell']" / "['Immune cell', 'DC cell']" 
                    '''

                    Neigboorhood = re.split(r'], ', cell['CD8 neighbors'])

                    '''
                    Other conditions are: 'CD20 neighbors' / 'CD11c neighbors'
                    '''
                    #replaces the brackets in the strings
                    Cell_list = []
                    for x in Neigboorhood:
                        Cell_list.append(x.replace('[', '').replace(']', ''))

                    #lists for cell counts
                    number_of_cells = len(Cell_list)

                    #list to hold the counts of neighboring cells
                    Tissue_count = []
                    Immune_count = []
                    CD4_count = []
                    CD8_count = []
                    B_count = []
                    Granulocyte_count = []
                    Myeloid_count = []
                    DC_count = []

                    #appends counters to the phenotype lists
                    for c in Cell_list:

                        if c == "'Tissue cell'":
                            Tissue_count.append('1 count')
                            
                        if c == "'Immune cell'":
                            Immune_count.append('1 count')

                        if c == "'Immune cell', 'CD4 T cell'":
                            Immune_count.append('1 count')
                            CD4_count.append('1 count')

                        if c == "'Immune cell', 'CD8 T cell'":
                            Immune_count.append('1 count')
                            CD8_count.append('1 count')

                        if c == "'Immune cell', 'B cell'":
                            Immune_count.append('1 count')
                            B_count.append('1 count')

                        if c == "'Immune cell', 'Granulocyte'":
                            Immune_count.append('1 count')
                            Granulocyte_count.append('1 count')

                        if c == "'Immune cell', 'Myeloid cell'":
                            Immune_count.append('1 count')
                            Myeloid_count.append('1 count')

                        if c == "'Immune cell', 'DC cell'":
                            Immune_count.append('1 count')
                            DC_count.append('1 count')

                    #calculates the percentages in regard to the neighboring cells
                    Tissue_percent = (len(Tissue_count) / len(Cell_list)) * 100
                    Immune_percent = (len(Immune_count) / len(Cell_list)) * 100
                    CD4_percent = (len(CD4_count) / len(Cell_list)) * 100
                    CD8_percent = (len(CD8_count) / len(Cell_list)) * 100
                    B_percent = (len(B_count) / len(Cell_list)) * 100
                    Granulocyte_percent = (len(Granulocyte_count) / len(Cell_list)) * 100
                    Myeloid_percent = (len(Myeloid_count) / len(Cell_list)) * 100
                    DC_percent = (len(DC_count) / len(Cell_list)) * 100

                    #appends the percentages of neighboring cells into a list for later analysis
                    neighborhood_Tissue_percent.append(Tissue_percent)
                    neighborhood_Immune_percent.append(Immune_percent)
                    neighborhood_CD4_percent.append(CD4_percent)
                    neighborhood_CD8_percent.append(CD8_percent)
                    neighborhood_CD20_percent.append(B_percent)
                    neighborhood_CD15_percent.append(Granulocyte_percent)
                    neighborhood_CD68_percent.append(Myeloid_percent)
                    neighborhood_DC_percent.append(DC_percent)

    #--------

    try:
        #basic visualization if a cell is more surrounded by other immune cells of tissue cells

        #settings for the neighborhood plot
        categories = ['Tissue', 'Immune']
        values = [statistics.mean(neighborhood_Tissue_percent),statistics.mean(neighborhood_Immune_percent)]
        errors = [statistics.stdev(neighborhood_Tissue_percent),statistics.stdev(neighborhood_Immune_percent)]
        colors = ['gray', '#44933a']

        #generates the plot
        fig, ax = plt.subplots()
        ax.bar(categories, values, yerr=errors, capsize=5, color=colors, edgecolor='black')
        ax.set_ylabel('Percent of Neighboring cells')
        ax.set_title('')
        ax.set_ylim(0, 110)

        #shows the plot and saves it
        plt.show()
        fig.savefig(folder_dir + "/CD8_neighborhood_vis",dpi=300)
        plt.close()

        #visualization of the surrounding phenotypes

        #settings for the neighborhood plot
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

        #generates the plot
        fig, ax = plt.subplots()
        ax.bar(categories, values, yerr=errors, capsize=5, color=colors, edgecolor='black')
        ax.set_ylabel('Percent of Neighboring cells')
        ax.set_title('')
        ax.set_ylim(0, 60)

        #shows the plot and saves it
        plt.show()
        fig.savefig(folder_dir + "/CD8_immune_neighborhood_vis", dpi=300)
        plt.close()

        #pandas data frame to export the neighborhood for further analysis
        neigboorhood_df = pd.DataFrame()
        neigboorhood_df['CD4'] = neighborhood_CD4_percent
        neigboorhood_df['CD8'] = neighborhood_CD8_percent
        neigboorhood_df['CD20'] = neighborhood_CD20_percent
        neigboorhood_df['CD15'] = neighborhood_CD15_percent
        neigboorhood_df['CD68'] = neighborhood_CD68_percent
        neigboorhood_df['CD11c'] = neighborhood_DC_percent

        neigboorhood_df.to_csv(folder_dir + "/" + "CD8_immune_neighborhood.csv") #exports the pandas data frame

    #handels errors and continues in case there is nothing to calculate
    except statistics.StatisticsError:
        print('------')
        print(filedir)
        print('This file hase nothing to calculate')
        print('------')
        continue
