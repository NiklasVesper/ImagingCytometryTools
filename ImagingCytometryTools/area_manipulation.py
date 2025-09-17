import pandas as pd
import ast
import shapely.geometry
import numpy as np

'''
Selects cells intersecting with specially chosen areas.
'''

#selects areas bases on pixel positive areas
def select_areas_by_pixel_positive_area(Cells, selected_areas):

    selected_cells = []
    for index, cell in Cells.iterrows():

        outline = ast.literal_eval(cell['Cell_outline'])
        inner_cell_area = ast.literal_eval(cell['Cell_pixel_area'])
        cell_area = outline + inner_cell_area

        for selected_area in selected_areas:

            pixel_positive_areas_to_check = ast.literal_eval(cell['pixel_positive_area_' + selected_area])

            for pixel_positive_area_to_check in pixel_positive_areas_to_check:

                pixel_positive_area_to_check_switch = [(pixel_position[1], pixel_position[0]) for pixel_position in pixel_positive_area_to_check]

                if bool(set(pixel_positive_area_to_check_switch) & set(cell_area)) == True:
                    selected_cells.append('in the area')
                    break

                else:
                    selected_cells.append('outside the area')

    Cells['Selected Areas'] = selected_cells

    return(Cells)

def select_areas_from_dictionary(Cells, selected_areas):

    selected_cells = []
    for index, cell in Cells.iterrows():

        cell_outline = ast.literal_eval(cell['Cell_outline'])
        area_to_check_key = str(cell['ImageName']) + '_' + str(cell['ImageNumber'])
        areas_to_check = selected_areas[area_to_check_key]
        
        if not areas_to_check:
            selected_cells.append('outside the area')

        else:

            cell_outline_polygon = shapely.geometry.Polygon(cell_outline)

            for area_to_check in areas_to_check:

                area_to_check_polygon = shapely.geometry.Polygon(area_to_check)

                if cell_outline_polygon.within(area_to_check_polygon) == True:

                    selected_cells.append('in the area')

                else:
                    selected_cells.append('outside the area')

    Cells['Selected Areas'] = selected_cells

    return(Cells)

def select_areas_from_both(Cells, selected_areas_dic, selected_areas_ppa):

    selected_cells = []
    for index, cell in Cells.iterrows():
        
        cell_outline = ast.literal_eval(cell['Cell_outline'])
        inner_cell_area = ast.literal_eval(cell['Cell_pixel_area'])
        cell_area = cell_outline + inner_cell_area
        
        area_to_check_key = str(cell['ImageName']) + '_' + str(cell['ImageNumber'])
        areas_to_check = selected_areas_dic[area_to_check_key]
        
        cell_outline_polygon = shapely.geometry.Polygon(cell_outline)

        validation_counter = 0
        
        for area_to_check in areas_to_check:

            area_to_check_polygon = shapely.geometry.Polygon(area_to_check)

            if cell_outline_polygon.within(area_to_check_polygon) == True:

                for selected_area_ppa in selected_areas_ppa:

                    pixel_positive_areas_to_check = ast.literal_eval(cell['pixel_positive_area_' + selected_area_ppa])

                    for pixel_positive_area_to_check in pixel_positive_areas_to_check:

                        pixel_positive_area_to_check_switch = [(pixel_position[1], pixel_position[0]) for pixel_position in pixel_positive_area_to_check]

                        if bool(set(pixel_positive_area_to_check_switch) & set(cell_area)) == True:
                            validation_counter = validation_counter + 1
                            
        if validation_counter == 0:
            selected_cells.append('outside the area')
            
        else:
            selected_cells.append('in the area')
            
    Cells['Selected Areas'] = selected_cells

    return(Cells)


    
    