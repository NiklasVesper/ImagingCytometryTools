if __name__ == "__main__":
  
# Pre-processing-----------------------------------------------------------------------------------------------------------
    folder_directory_with_spill_mat = r'C:\Users\vesper\Desktop\Spillovermatrix_IMC_2.txt'
    folder_directory_with_data = r'D:\ATF6\Round 1' #Careful if the spillover matrix and data is in the same directory, you will get an error message!

    # Compensation---------------------------------------------------------------------------------------------------------
    from ImagingCytometryTools.run_compensation_IMC import run_compensation_IMC


    run_compensation_IMC(folder_directory_with_data,
                         folder_directory_with_spill_mat,
                         use_GPU=False)


    # Image and folder generation------------------------------------------------------------------------------------------
    from ImagingCytometryTools.generate_images_from_IMC_data import generate_images_from_IMC_data
    from ImagingCytometryTools.generate_folders_for_CellProfiler import generate_folders_for_CellProfiler


    generate_images_from_IMC_data(folder_directory_with_data, 
                                  generate_ome_tiff=False,
                                  change_name_for_MCD_viewer=False)
    
    generate_folders_for_CellProfiler(folder_directory_with_data,
                                      generate_folders_per_mcd_file=True)
