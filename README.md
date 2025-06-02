ImagingCytometryTools

The main problem of tissue data analysis is the segmentation of cells.
Generalist segmentation algorithms like Cellpose and Deepcell significantly improved on that,
but these algorithms only allow cellular and/or nuclear segmentation that are independent from each other.

ImagingCytometryTools is a collection of scripts for the analysis of multiplex tissue imaging data,
with a focus on Imaging Mass Cytometry, but it can be adapted for any other imaging method.

The foundation the analysis is a CellProfiler pipeline with Cellpose as a plugin: 
(https://github.com/CellProfiler/CellProfiler-plugins/tree/master)



The main advantage of this analysis is that it allowes the user to generate subcellular information out of the cellular segmentation made by Cellpose.


An ImagingCytometryTools Python package and an improved version of the associated analysis Cellprofiler pipeline is still under development.

If you are interested in using this pipeline or want to help or contribute to its development you can contact me under:
niklas.vesper@uniklinik-freiburg.de
