<!-- ImagingCytometryTools -->
## ImagingCytometryTools

A key problem of high throughput multiplex tissue imaging data analysis is the segmentation of cells. Generalist segmentation algorithms like Cellpose and Deepcell significantly 
improved on that, but these algorithms only allow cellular and/or nuclear segmentation that are independent of each other (1,2). 

ImagingCytometryTools is a collection of scripts that allows the user to generate subcellular information out of a Cellpose segmentation in conjunction with a CellProfiler pipeline via a novel algorithm (3). The primary advantage of this approach is an improved characterization of transcription factor translocation or the localization of oncogenes. The pipeline has a focus on imaging mass cytometry, but it can be adapted for any other imaging method.

The foundation of the analysis is a <a href="https://github.com/CellProfiler/CellProfiler-plugins/tree/master"><strong>CellProfiler pipeline with Cellpose as a 
plugin</strong></a>.

**If you use this pipeline or its associated package please cite:** Li, Lebeaupin, et al.

An ImagingCytometryTools python package with more functionalities and an improved version of the associated CellProfiler pipeline are under active development 
and soon to be released.

To use a preliminary package install <a href="https://git-scm.com/downloads"><strong>git</strong></a> and run:

    pip install git+https://github.com/NiklasVesper/ImagingCytometryTools.git

If you are interested in using this pipeline or want to help or contribute to its development you can contact me under:<br />
niklas.vesper@uniklinik-freiburg.de

## Visual explanation of the algorithm 

![](https://github.com/NiklasVesper/ImagingCytometryTools/blob/main/Images/algorithm.png)

## References 

(1) Pachitariu, M., & Stringer, C. (2022). Cellpose 2.0: how to train your own model. Nature methods, 19(12), 1634-1641.

(2) Greenwald, N. F., Miller, G., Moen, E., Kong, A., Kagel, A., Dougherty, T., ... & Van Valen, D. (2022). Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning. Nature biotechnology, 40(4), 555-565.

(3) Stirling, D. R., Swain-Bowden, M. J., Lucas, A. M., Carpenter, A. E., Cimini, B. A., & Goodman, A. (2021). CellProfiler 4: improvements in speed, utility and usability. BMC bioinformatics, 22(1), 433.
