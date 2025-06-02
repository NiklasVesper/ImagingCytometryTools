<h3 align="left">ImagingCytometryTools</h3>

The key problem of multiplex tissue imaging data analysis is the segmentation of cells.
Generalist segmentation algorithms like Cellpose and Deepcell significantly improved on that,
but these algorithms only allow cellular and/or nuclear segmentation that are independent from each other.

ImagingCytometryTools is a collection of scripts that allows the user to generate subcellular information out of a cellular segmentation. The pipeline has a focus on Imaging Mass Cytometry, but it can be adapted for any other imaging method.

The foundation the analysis is a <a href="https://github.com/CellProfiler/CellProfiler-plugins/tree/master"><strong>CellProfiler pipeline with Cellpose as a plugin.</strong></a>

An ImagingCytometryTools Python package with more functionalities and an improved version of the associated analysis Cellprofiler pipeline are under active development.

If you are interested in using this pipeline or want to help or contribute to its development you can contact me under:
niklas.vesper@uniklinik-freiburg.de
