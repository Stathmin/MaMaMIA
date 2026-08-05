# MaMaMIA

`MaMaMIA` (Max Mapping for Meticulous Introgression Analysis) detects interspecific introgressions from Genotyping-by-Sequencing (GBS) coverage data.

It builds a Read Count Array (`RCA`) from windowed donor and recipient coverage and GC-content tables, corrects GC bias with a zero-inflated negative binomial model, and segments pairwise coverage differences with circular binary segmentation (CBS) into an Introgression Segment Array (`ISA`). Putative introgressions can be plotted, filtered, and exported as BEDPE.

`MaMaMIA` ships with an example wheat dataset (`triticum`: *T. timopheevii* donor vs *T. aestivum* Chinese Spring T2T recipient). See `vignette("Overview", package = "MaMaMIA")` for a full walkthrough.

## Installation

You can install `MaMaMIA` from [GitHub](https://github.com/alermol/mamamia) with:

```r
# install.packages("devtools")
devtools::install_github("alermol/mamamia")
```

## Disclaimer

This package is still under active development, the content is therefore subject to change. 

## Contact

Suggestions and bug reports: please [open an issue](https://github.com/alermol/mamamia/issues).
