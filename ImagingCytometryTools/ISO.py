import pandas as pd
import ast
import shapely.geometry
import numpy as np

'''
Function that tries to correct cells containing marker combinations that are not possible.

This is done by cecking neighboring cells whether they are assigned a cell type that could explain such a marker combination. 
Overlap is then verified by checking whether both cells share the same positive pixel area. 
'''

def ISO(Cells, allowed_distance = 10, allowed_iterations = 3):

    #adds new columns to store the information
    Cells['cell_types_and_states_ISO'] = Cells['cell_types_and_states']
    Cells['ISO_visualisation_and_information'] = np.nan

    fixed_cell_count = 42 #42? Yes, yes I thought it over quite thoroughly its 42.
    final_count = [] #list for the final cell count
    iteration = 1 #iteration counter

    while fixed_cell_count >=1 and iteration <= allowed_iterations: #runs the loop till no new cells are found or the maximum of iterations is reached

        if iteration == 1:
            cell_count = 0 #counter for the total cellcount that can be fixed

        fixed_cell_count = 0

        for index, cell in Cells.iterrows(): # goes throw the cells

            cell_types_and_state = ast.literal_eval(cell['cell_types_and_states_ISO']) #importsd cells sates

            if len(cell_types_and_state) >= 3: #finds all the cells

                #counts the total amount of cells only during the first iteration
                if iteration == 1:
                    cell_count = cell_count + 1

                #print('start cell: ' + str(cell_types_and_state))

                outline = ast.literal_eval(cell['Cell_outline']) #gets the cell outline
                inner_cell_area = ast.literal_eval(cell['Cell_pixel_area']) #gets the inner cell area
                cell_area = outline + inner_cell_area #creates the total cell area by combining them

                #generates the area to search for the allowed overlap
                cell_outline = shapely.geometry.Polygon(np.array(outline))
                cell_allowed_overlap = cell_outline.buffer(allowed_distance)

                cells_with_allowed_overlap = [] #empty list for the cells with the allowed overlap
                for index_other_cell, other_cell in Cells.iterrows(): #goes throw the other cells

                    if cell['ImageNumber'] == other_cell['ImageNumber']: #takes the cells that have the same image number
                        outline_other_cell = ast.literal_eval(other_cell['Cell_outline']) #gets the cell outline
                        other_cell_area = shapely.geometry.Polygon(np.array(outline_other_cell)) #creates a polygon out of the cell outline

                        if cell_allowed_overlap.intersects(other_cell_area) == True: #checks if the cell has the allowed overlap
                            other_cell_types_and_state = ast.literal_eval(other_cell['cell_types_and_states_ISO']) #gets the cell type of the cell

                            #if the cell type is determined it can be uses to check whether there s overlap
                            if len(other_cell_types_and_state) == 2:
                                cells_with_allowed_overlap.append(index_other_cell)

                #empty lists for the final evaluations
                cells_who_explain_mixture = []
                cells_who_explain_mixture_information = []

                cell_types_and_state_ckeck_list = [cell_lineage[0] for cell_lineage in cell_types_and_state[:-1]] #creates a list for all cell lineage markers within one cell

                #checks each cell if the can explain the marker combination
                for cell_with_allowed_overlap in cells_with_allowed_overlap:

                    cell_lineage_to_check = ast.literal_eval(Cells['cell_types_and_states_ISO'].loc[cell_with_allowed_overlap])[0][0] #gets the cell lineage marker form the cell that needs to be checked

                    if cell_lineage_to_check in cell_types_and_state_ckeck_list: #Verifies that a cell in close has the correct marker

                        #gets the cellular area of the cell
                        outline_other_cell_to_ckeck = ast.literal_eval(Cells.loc[cell_with_allowed_overlap,'Cell_outline'])
                        inner_cell_area_to_ckeck = ast.literal_eval(Cells.loc[cell_with_allowed_overlap,'Cell_pixel_area'])
                        cell_area_to_ckeck = outline_other_cell_to_ckeck + inner_cell_area_to_ckeck

                        pixel_positive_areas_to_check = ast.literal_eval(Cells.loc[cell_with_allowed_overlap,'pixel_positive_area_' + cell_lineage_to_check[:-1]]) #gets the correct pixel positive area

                        #goes throw all the pixel positive areas and checks if one can be found that explains the marker mixture
                        for pixel_positive_area_to_check in pixel_positive_areas_to_check:

                            pixel_positive_area_to_check_switch = [(pixel_position[1], pixel_position[0]) for pixel_position in pixel_positive_area_to_check]

                            #checks whether a pixel positive area is found that is shared by both cells
                            if bool(set(pixel_positive_area_to_check_switch) & set(cell_area_to_ckeck)) == True and bool(set(pixel_positive_area_to_check_switch) & set(cell_area)) == True:

                                cells_who_explain_mixture.append(cell_lineage_to_check[:-1])
                                cells_who_explain_mixture_information.append([cell_lineage_to_check[:-1],cell_with_allowed_overlap,pixel_positive_area_to_check_switch])

                #if all cells needed tho explain the mixture are found the new cell type is assigned
                if len(set(cells_who_explain_mixture)) == len(cell_types_and_state_ckeck_list)-1:
                    fixed_cell_count = fixed_cell_count + 1
                    real_cell_type = []
                    list_mixture = [cell_who_explain_mixture[0] + '+' for cell_who_explain_mixture in cells_who_explain_mixture_information]
                    list_og = [cell_type_and_state[0] for cell_type_and_state in cell_types_and_state[:-1]]
                    list_real = list(set(list_mixture).symmetric_difference(set(list_og)))

                    for cell_type_and_state in cell_types_and_state:

                        if list_real[0] == cell_type_and_state[0]:
                            real_cell_type.append(cell_type_and_state)
                            real_cell_type.append(['CD45+'])

                    #append the corrected cell types into a new column
                    Cells.loc[index,'cell_types_and_states_ISO'] = str(real_cell_type)
                    Cells.loc[index, 'ISO_visualisation_and_information'] = str(cells_who_explain_mixture_information)

        iteration = iteration + 1
        final_count.append(fixed_cell_count)

    print('The total number of cells that needed to be fixed were: ' + str(cell_count) + ' and the number of cells that were fixed are: ' + str(final_count))

    return(Cells)