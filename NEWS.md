# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

# HistogramZoo v1.4.0

## [1.4.2] - 2023-07-03

### Added
- Add `+.Histogram` support to add bin width 1 histograms together.

## [1.4.1] - 2023-05-02

### Added
- ModelFit class to allow pretty-printing of `fit_distributions` results

### Fixed
- Bugs associated with fitting truncated distributions

## [1.4.0] - 2023-04-19

### Added
- S3 methods to allow numeric/Histogram/GenomicHistogram input in key segment_and_fit functions [Pull #121](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/121)
  - consecutive_start and consecutive_end for GenomicHistogram to keep track of disjoint intervals
  - shifted gamma function to allow for non-zero starts

## Changed
- Added summary statistics estimation for GenomicHistogram objects [Pull #124](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/124)
   - Added `bin_width` and `midpoint` estimation for Histogram and GenomicHistogram objects
- Enforcing continuity of bins in Histogram and GenomicHistogram objects when subsetting [Pull #125](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/125)

## Fixed
- Fix vignette/README/image file paths for pkgdown [Pull #115](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/115)

# HistogramZoo v1.3.0

## Added
- Added summary statistics (mean, sd, var, skew) estimation from histogram objects [Pull #113](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/113)

## Changed
- Enforcing bin width and addition of intron-traversing bins in GenomicHistogram [Pull #108](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/108)

# HistogramZoo v1.2.0

## Added
- Added a classic histogram plotting function [Pull #111](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/111)

# HistogramZoo v1.1.0

## Added
- Maxiumum likelihood estimation for distribution parameters [Pull #97](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/97)
- A vignette! [Pull #81](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/81)
- Added `gamma_flip` as an available distribution [Pull #62](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/62)

## Changed
- New rules for `Histogram` and `GenomicHistogram` [Pull #76](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/76)
  - `Histogram` uses a strictly continuous Base 0 system while `GenomicHistogram` uses a potentially disjoint Base 1 system
  - Both are restricted to having a fixed `bin_width` in (n-1) bins
- Add pkgdown build on new release [Pull #56](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/56)
- Updated `find_consensus_model` to include user-specified weights and metric prioritization & addition of new RRA method [Pull #53](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/53)

## Fixed
- Bug in `segment_and_fit` - finer segmentation is required to initialize `meaningful_gaps_local` [Pull #73](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/73)
- Refactored `create_trackplot` to use `BoutrosLab.plotting.general::create.scatterplot` as a base function thus increasing efficiency [Pull #71](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/71)
- Refactored `find_uniform_segment` to reduce code redundancy [Pull #59](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/59)

# HistogramZoo v1.0.0

## Added
- First release!
