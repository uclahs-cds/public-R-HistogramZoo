# ConsensusPeaks

ConsensusPeaks is a generalized framework for histogram segmentation and statistical characterization written in R.

## Installation

Using devtools in R:
```R
library(devtools)
install_github('https://github.com/uclahs-cds/public-R-ConsensusPeaks')
```

From source:
```shell script
git clone https://github.com/uclahs-cds/public-R-ConsensusPeaks.git
R CMD INSTALL public-R-ConsensusPeaks
```

## Quick Start
```R

library(ConsensusPeaks)

x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))

results = segment_and_fit(x, eps = 0.005)

results_table = summarize_results(results)

```

## Additional features
ConsensusPeaks supports common genomic data types including bigWig files and BED files. ConsensusPeaks supports both genomic and transcriptomic coordinate systems. For more information, refer to the ConsensusPeaks [vignette](https://github.com/uclahs-cds/public-R-ConsensusPeaks/tree/main/vignettes/ConsensusPeaks.html).

## Publication

Coming soon!
