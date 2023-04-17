
# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added

### Changed
- S3 methods to allow numeric/Histogram/GenomicHistogram input in key segment_and_fit functions

### Fixed
- Minor fixes to vignette file path locations for successful pkgdown build

## [1.3.0] - 2023-04-14

### Added
- Added summary statistics (mean, sd, var, skew) estimation from histogram objects [Pull #113](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/113)

### Changed
- Enforcing bin width and addition of intron-traversing bins in GenomicHistogram [Pull #108](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/108)

## [1.2.0] - 2023-04-12

### Added
- Added a classic histogram plotting function [Pull #111](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/111)

## [1.1.0] - 2023-03-09

### Added
- Maxiumum likelihood estimation for distribution parameters [Pull #97](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/97)
- A vignette! [Pull #81](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/81)
- Added `gamma_flip` as an available distribution [Pull #62](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/62)

### Changed
- New rules for `Histogram` and `GenomicHistogram` [Pull #76](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/76)
  - `Histogram` uses a strictly continuous Base 0 system while `GenomicHistogram` uses a potentially disjoint Base 1 system
  - Both are restricted to having a fixed `bin_width` in (n-1) bins
- Add pkgdown build on new release [Pull #56](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/56)
- Updated `find_consensus_model` to include user-specified weights and metric prioritization & addition of new RRA method [Pull #53](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/53)

### Fixed
- Bug in `segment_and_fit` - finer segmentation is required to initialize `meaningful_gaps_local` [Pull #73](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/73)
- Refactored `create_trackplot` to use `BoutrosLab.plotting.general::create.scatterplot` as a base function thus increasing efficiency [Pull #71](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/71)
- Refactored `find_uniform_segment` to reduce code redundancy [Pull #59](https://github.com/uclahs-cds/public-R-HistogramZoo/pull/59)

## [1.0.0] - 2022-08-22

### Added
- First release!
