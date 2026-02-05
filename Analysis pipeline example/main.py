if __name__ == "__main__":
  
# Pre-processing-----------------------------------------------------------------------------------------------------------
    folder_directory_with_spill_mat = r'path\Spillovermatrix_IMC.txt'
    folder_directory_with_data = r'path' #Careful if the spillover matrix and data is in the same directory, you will get an error message!

    # Compensation---------------------------------------------------------------------------------------------------------
    from ImagingCytometryTools.run_compensation_IMC import run_compensation_IMC


    run_compensation_IMC(folder_directory_with_data,
                         folder_directory_with_spill_mat,
                         use_GPU=False)


    # Image and folder generation------------------------------------------------------------------------------------------
    from ImagingCytometryTools.generate_images_from_IMC_data import generate_images_from_IMC_data
    from ImagingCytometryTools.generate_folders_for_CellProfiler import generate_folders_for_CellProfiler
    from ImagingCytometryTools.split_images_for_CellProfiler import split_images_for_CellProfiler

    generate_images_from_IMC_data(folder_directory_with_data, 
                                  generate_ome_tiff=False,
                                  change_name_for_MCD_viewer=False)
    
    generate_folders_for_CellProfiler(folder_directory_with_data,
                                      generate_folders_per_mcd_file=True)

    split_images_for_CellProfiler(folder_directory_with_data,
                                  image_size_threshold=600000)

# Segmentation and general image processing--------------------------------------------------------------------------------

    '''
    For the Segmentation and any additional image processing a CellProfiler pipeline is used.
    In the future there also will be an end to end solution avalible.
    An example pipeline can be found under: https://github.com/NiklasVesper/ImagingCytometryTools/tree/main/Cellprofiler%20pipelines

    For this CellProfiler pipeline Cellpose 2 is required as a plugin. Instructions for the plugin installation can be found under:
    https://github.com/CellProfiler/CellProfiler-plugins/tree/master

    CellProfiler: https://doi.org/10.1186/s12859-021-04344-9
    Cellpose 2: https://doi.org/10.1038/s41592-022-01663-4
    '''

# Post-processing----------------------------------------------------------------------------------------------------------

    # Generating neighborhood----------------------------------------------------------------------------------------------
    from ImagingCytometryTools.run_generate_neighborhood import run_generate_neighborhood
    from ImagingCytometryTools.run_test_outline_identification_and_matching import run_test_outline_identification_and_matching

    distance = 5
    file_neighborhood = r'RunCellpose_C'
    output_folder_neighborhood = r'full cell neighborhood ' + str(distance)
    run_generate_neighborhood(folder_directory_with_data,
                              file_neighborhood,
                              output_folder_neighborhood,
                              add_cellular_information_and_outline=True,
                              use_own_neighborhood_radius=[True, distance])
    
    file_test = r'full_cell_neighborhood_5'
    output_folder_neighborhood_analysis_test = r'neighborhood test outline'
    run_test_outline_identification_and_matching(folder_directory_with_data,
                                                 file_test,
                                                 output_folder_neighborhood_analysis_test,
                                                 show_individual_cells=False,
                                                 show_nucleus=False)

    file_test = r'full_cell_neighborhood_5'
    output_folder_neighborhood_analysis_test = r'neighborhood test outline individual cells'
    run_test_outline_identification_and_matching(folder_directory_with_data,
                                                 file_test,
                                                 output_folder_neighborhood_analysis_test,
                                                 show_individual_cells=True,
                                                 show_nucleus=False)
  
  # Generating subcellular localization and neighborhood-----------------------------------------------------------------
    from ImagingCytometryTools.run_generate_subcellular_localization import run_generate_subcellular_localization
    from ImagingCytometryTools.run_test_outline_identification_and_matching import run_test_outline_identification_and_matching

    file_subcellular_location = [r'RunCellpose_C', 'RunCellpose_N', 'Cytoplasm']
    output_folder_subcellular = r'subcellular information'
    run_generate_subcellular_localization(folder_directory_with_data,
                                          file_subcellular_location,
                                          output_folder_subcellular,
                                          mode='basic',
                                          nucleus_count=1)
    
    distance = 5
    file_neighborhood_subcell = r'subcellular_information'
    output_folder_neighborhood_subcell = r'subcellular and neighborhood ' + str(distance)
    run_generate_neighborhood(folder_directory_with_data,
                              file_neighborhood_subcell,
                              output_folder_neighborhood_subcell,
                              add_cellular_information_and_outline=False,
                              use_own_neighborhood_radius=[True, distance])
    
    file_test = r'subcellular_and_neigboorhood_5'
    output_folder_subcellular_localisation_analysis_test = r'subcellular test outline'
    run_test_outline_identification_and_matching(folder_directory_with_data,
                                                 file_test,
                                                 output_folder_subcellular_localisation_analysis_test,
                                                 show_individual_cells=False,
                                                 show_nucleus=False)

    file_test = r'subcellular_and_neigboorhood_5'
    output_folder_subcellular_localisation_analysis_test_nucleus = r'subcellular test outline individual cells'
    run_test_outline_identification_and_matching(folder_directory_with_data,
                                                 file_test,
                                                 output_folder_subcellular_localisation_analysis_test_nucleus,
                                                 show_individual_cells=True,
                                                 show_nucleus=True)

    # Phenotyping----------------------------------------------------------------------------------------------------------
    from ImagingCytometryTools.run_assign_cell_types_and_states import run_assign_cell_types_and_states
    from ImagingCytometryTools.run_generate_image_galleries import run_generate_image_galleries

    
    file_phenotypes = r'full_cell_neighborhood_5'
    output_folder_phenotypes = r'full cell phenotypes and neighborhood 5'
    run_assign_cell_types_and_states(folder_directory_with_data,
                                     file_phenotypes,
                                     output_folder_phenotypes,
                                     overwrite_existing_files=True)

    file_phenotypes = r'full_cell_phenotypes_and_neighborhood_5'
    analysis = r'cell_types_and_states'
    output_folder_images = r'image gallery CD8+ ATF6 full image'
    max_pixel_images = ['Eu151Di','Pt194Di'] #'Dy162Di','Pt194Di'
    proteins = ['ATF6','Nucleus_HH3'] #'CD8','Nucleus_HH3'
    contrast_multipliers = [1750,3000] #2000,3000
    run_generate_image_galleries(folder_directory_with_data,
                                 file_phenotypes,
                                 analysis,
                                 output_folder_images,
                                 max_pixel_images,
                                 proteins,
                                 normalisation_value='mean',
                                 contrast_multiplier=contrast_multipliers,
                                 select_cell_type_and_state=[True,['CD8+',['ATF6+'],[]]],
                                 add_mixed_cells=True,
                                 generate_crops=False,
                                 crop_size=40,
                                 show_neighborhood_radius=True,
                                 select_neighboring_cell_type_and_state=[False,['CD20+',[],[]]],
                                 show_neighboring_cells=False)
