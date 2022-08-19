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

set.seed(271828)

# Generating Data
my_data <- rnorm(10000, mean = 50, sd = 20)
histogram_data <- observations_to_histogram(my_data, histogram_bin_width=5)

# Segmentation and fitting distributions
results <- segment_and_fit(histogram_data, eps = 1)

# Summarizing results
results_table <- summarize_results(results)

# Plotting
create_coverageplot(results)

```

## Additional features
HistogramZoo supports common genomic data types including bigWig files and BED files. HistogramZoo supports both genomic and transcriptomic coordinate systems.
