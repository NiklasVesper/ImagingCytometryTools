from ImagingCytometryToolsIMC.run_compensation_IMC import run_compensation_IMC
from ImagingCytometryToolsIMC.generate_images_from_IMC_data import generate_images_from_IMC_data

from ImagingCytometryToolsIMC.run_neighborhood_analysis import run_neighborhood_analysis
from ImagingCytometryToolsIMC.run_subcellular_location_analysis import run_subcellular_location_analysis

from ImagingCytometryToolsIMC.run_assign_cell_types import run_assign_phenotypes
from ImagingCytometryToolsIMC.generate_image_galaries import generate_image_galaries

from ImagingCytometryToolsIMC.export_pipeline_settings import export_pipeline_settings

'''

'''

file_dir_comp_mat = r'C:\Users\vesper\Desktop\spillover_ira_3.txt' # imports spillover matrix
folder_dir_data = r'F:\NiklasILoveYou' # folder directory

#To do (add row comp/improve file handling?)
run_compensation_IMC(folder_dir_data, file_dir_comp_mat, use_GPU = True)
#run_compensation_MIBI
#run_compensation_FI

generate_images_from_IMC_data(folder_dir_data, generate_ome_tiff = False) #(Done)

#-------------------------------------------Segmentation----------------------------------------------------------------
'''
https://github.com/NiklasVesper/ImagingCytometryTools/tree/main/Cellprofiler%20pipelines

To do (add example pipelines)
'''

#To do (add outlines / remove unnececary chnnels)
file_neighborhood = r'RunCellpose_C'
output_folder_neighborhood = r'neighborhood'
run_neighborhood_analysis(folder_dir_data,file_neighborhood,output_folder_neighborhood, use_own_neighborhood_radius = [True, 50])

#To do (advanced/two cells)
file_subcellular_location = [r'RunCellpose_C','RunCellpose_N','Cytoplasm']
output_folder_subcellular_location = r'subcellular localisation'
run_subcellular_location_analysis(folder_dir_data,file_subcellular_location,output_folder_subcellular_location, mode = 'basic')

file_neighborhood_subcell = r'subcellular_localisation'
output_folder_neighborhood_subcell = r'subcellular and neigboorhood'
run_neighborhood_analysis(folder_dir_data,file_neighborhood_subcell,output_folder_neighborhood_subcell, add_outline = False) #(Done)

file_phenotypes = r'subcellular_and_neigboorhood'
output_folder_phenotypes = r'subcellular and neigboorhood phenotypes'
#run_assign_phenotypes(folder_dir_data,file_phenotypes,output_folder_phenotypes) #(Done)

#To do (RGB + outline + neiboorhood)
#single_file = 'No'
#generate_image_galaries(folder_dir_data,'phenotypes')

#set_pixel_positive_area_threshholds(folder_dir_data,Marker,Area,Intensety,Pixel_distance)
#run_ISO()

#single_file = 'Yes'
#generate_image_galaries(folder_dir_data)

#run_assign_neiboohood_phenotypes()

#run_phenotype_statistics()
#run_neiboohood_statistics()
#run_area_statistics()

#from assign_cell_types import cell_type_thresholds
#export_pipeline_settings(cell_type_thresholds)