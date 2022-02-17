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

x = c(0, 0, 1, 2, 3, 2, 1, 0, 0)

histograms = segment.fit.agnostic(x)

results = summarize.results.agnostic(histograms)

print(results)

```

## Vignette
For more information, refer to the ConsensusPeaks [vignette](https://github.com/uclahs-cds/public-R-ConsensusPeaks/tree/main/vignettes/ConsensusPeaks.html).

## Publication

Coming soon!
