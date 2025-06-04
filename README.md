<!-- ImagingCytometryTools -->
## ImagingCytometryTools

A key problem of multiplex tissue imaging data analysis is the segmentation of cells. Generalist segmentation algorithms like Cellpose and Deepcell significantly improved on that, but these algorithms only allow cellular and/or nuclear segmentation that are independent of each other. ImagingCytometryTools is a collection of scripts that allows the user to generate subcellular information out of a Cellpose segmentation via a novel algorithm. The primary advantage of this approach is an improved characterization of transcription factor translocation. The pipeline has a focus on Imaging Mass Cytometry, but it can be adapted for any other imaging methods.

The foundation the analysis is a CellProfiler pipeline with <a href="https://github.com/CellProfiler/CellProfiler-plugins/tree/master"><strong>Cellpose as a plugin.</strong></a>

If you use this pipeline or its associated package please cite:
Li, Lebeaupin, et al.

To use the preliminary package install git and run:

    pip install git+https://github.com/NiklasVesper/ImagingCytometryTools.git

An ImagingCytometryTools python package with more functionalities and an improved version of the associated analysis CellProfiler pipeline are under active development and soon to be released.

If you are interested in using this pipeline or want to help or contribute to its development you can contact me under:<br />
niklas.vesper@uniklinik-freiburg.de
