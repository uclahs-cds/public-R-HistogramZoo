# HistogramZoo

HistogramZoo is a generalized framework for histogram segmentation and statistical characterization written in R.

## Installation

Using devtools in R:
```R
library(devtools)
install_github('https://github.com/uclahs-cds/public-R-HistogramZoo')
```

From source:
```shell script
git clone https://github.com/uclahs-cds/public-R-HistogramZoo.git
R CMD INSTALL public-R-HistogramZoo
```

## Quick Start
```R

library(HistogramZoo)

x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))

results = segment_and_fit(x, eps = 0.005)

results_table = summarize_results(results)

```

## Additional features
HistogramZoo supports common genomic data types including bigWig files and BED files. HistogramZoo supports both genomic and transcriptomic coordinate systems. For more information, refer to the HistogramZoo [vignette](https://github.com/uclahs-cds/public-R-HistogramZoo/tree/main/vignettes/HistogramZoo.html).
