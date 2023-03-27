# HistogramZoo

1. [Description](#description)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Resources](#resources)
5. [Getting help](#getting-help)
6. [Citation information](#citation-information)
7. [License](#license)

## Description
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

From CRAN (coming soon!)

## Quick start

HistogramZoo supports common genomic data types including bigWig files and BED files. HistogramZoo supports both genomic and transcriptomic coordinate systems.

```R

library(HistogramZoo);

set.seed(271828);

# Generating Data
gaussian_data <- rnorm(10000, mean = 50, sd = 15);
histogram_data <- observations_to_histogram(gaussian_data, histogram_bin_width=5);

# Segmentation and fitting distributions
results <- segment_and_fit(histogram_data, eps = 1);

# Summarizing results
results_table <- summarize_results(results);

# Plotting
create_coverageplot(results);

```
|   region_id  |   segment_id  |   start  |   end  |   interval_count  |   interval_sizes  |   interval_starts  |   histogram_start  |   histogram_end  |   value  |   metric  |   dist               |   dist_param1  |   dist_param2  |   dist_param1_name  |   dist_param2_name  |         |       |
|--------------|---------------|----------|--------|-------------------|-------------------|--------------------|--------------------|------------------|----------|-----------|----------------------|----------------|----------------|---------------------|---------------------|---------|-------|
|   -21-102    |   1           |   19     |   79   |   1               |   61              |                    |   1                |                  |   9      |   20      |   0.955492468663381  |   consensus    |   norm         |   33.032410958755   |   14.881052093018   |   mean  |   sd  |

![coverage plot](readme_imports/hz_output.pdf)

## Resources

The function reference and a detailed vignette are available [here](https://uclahs-cds.github.io/public-R-HistogramZoo/). A CRAN page is currently in progress.

## Getting help

Looking for guidance or support with HistogramZoo? Look no further.

* Check out our [Discussions](https://github.com/uclahs-cds/public-R-HistogramZoo/discussions) page!
* Submit bugs :bug:, suggest new features :cherry_blossom: or see current work :mechanical_arm: at our [Issues](https://github.com/uclahs-cds/public-R-HistogramZoo/issues) page.

## Citation information

You have stumbled upon an unpublished software :shushing_face: :shushing_face: :shushing_face:. We are currently preparing the manuscript for HistogramZoo. Please befriend us to learn more or check back later for updated citation information.

## Licence

Authors: Helen Zhu, Stefan Eng, Paul C. Boutros (PBoutros@mednet.ucla.edu)

HistogramZoo is licensed under the GNU General Public License version 3.0. See the file LICENSE.md for the terms of the GNU GPL license.

HistogramZoo is a generalized framework for histogram segmentation and statistical characterization written in R.

Copyright (C) University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
