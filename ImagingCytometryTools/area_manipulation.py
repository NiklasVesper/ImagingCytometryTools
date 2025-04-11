import pandas as pd
import ast
import shapely.geometry
import numpy as np

'''
Selects cells intersecting with specially chosen areas.
'''

def select_areas(Cells,selected_areas):

    Selected_Cells = pd.DataFrame().reindex_like(Cells)  # data frame of the cells selected for a specific neighbor
    Selected_Cells.drop(Selected_Cells.index, inplace=True)

    for selected_area in selected_areas:

        for index, cell in Cells.iterrows(): # goes throw the cells

            outline = ast.literal_eval(cell['Cell_outline'])  # gets the cell outline
            inner_cell_area = ast.literal_eval(cell['Cell_pixel_area'])  # gets the inner cell area
            cell_area = outline + inner_cell_area  # creates the total cell area by combining them
            pixel_positive_areas_to_check = ast.literal_eval(cell['pixel_positive_area_'+selected_area])  # gets the cell outline

            for pixel_positive_area_to_check in pixel_positive_areas_to_check:

                pixel_positive_area_to_check_switch = [(pixel_position[1], pixel_position[0]) for pixel_position in pixel_positive_area_to_check]

                if bool(set(pixel_positive_area_to_check_switch) & set(cell_area)) == True:

                    Selected_Cells.loc[len(Selected_Cells)] = cell

                    break

    return(Selected_Cells)