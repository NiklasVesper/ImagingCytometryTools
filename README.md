<!-- ImagingCytometryTools -->
## ImagingCytometryTools (beta-1.0.2)

A key problem of high throughput multiplex tissue imaging data analysis is the segmentation of cells. Generalist segmentation algorithms like Cellpose and Deepcell significantly 
improved on that, but these algorithms only allow cellular and/or nuclear segmentation that are independent of each other (1, 2). 

ImagingCytometryTools is a collection of scripts that allows the user to generate subcellular information out of a Cellpose segmentation in conjunction with a CellProfiler pipeline via a novel algorithm (3, 4). The primary advantage of this approach is an improved characterization of transcription factor translocation or the localization of oncogenes. The pipeline has a focus on imaging mass cytometry, but it can easily be used for any other imaging method.

The foundation of the analysis is a <a href="https://github.com/CellProfiler/CellProfiler-plugins/tree/master"><strong>CellProfiler pipeline with Cellpose as a 
plugin</strong></a>.

**If you use this pipeline or its associated package please cite:** <a href="https"><strong>Li, Lebeaupin, et al. (Nature 2025, in press)</strong></a>

An ImagingCytometryTools python package with more functionalities and an improved version of the associated CellProfiler pipeline are under active development 
and soon to be released.

To use a preliminary package install <a href="https://git-scm.com/downloads"><strong>git</strong></a> and run:

    pip install git+https://github.com/NiklasVesper/ImagingCytometryTools.git

If you are interested in using this pipeline, have any questions or want to help or contribute to its development you can contact me under:
niklas.vesper@uniklinik-freiburg.de

## Visual explanation of the algorithm to generate subcellular information 

![](https://github.com/NiklasVesper/ImagingCytometryTools/blob/main/Images/algorithm.png)

## References 

(1) Pachitariu, M., & Stringer, C. (2022). Cellpose 2.0: how to train your own model. Nature methods, 19(12), 1634-1641.

(2) Greenwald, N. F., Miller, G., Moen, E., Kong, A., Kagel, A., Dougherty, T., ... & Van Valen, D. (2022). Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning. Nature biotechnology, 40(4), 555-565.

(3) Stirling, D. R., Swain-Bowden, M. J., Lucas, A. M., Carpenter, A. E., Cimini, B. A., & Goodman, A. (2021). CellProfiler 4: improvements in speed, utility and usability. BMC bioinformatics, 22(1), 433.

(4) Li, Lebeaupin, et al. (Nature 2025, in press)

## Release history

1.0.1: Added image compensation, image generation and helper functions to set up the CellProfiler pipeline.<br />
1.0.2: Updated generation of neighborhood and subcellular data as wella test functions  
