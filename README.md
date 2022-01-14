# exNGT

### version 0.1

---

## Overview

**exNGT** is an R software project to analyze the new NGT-2012 stable oxygen
isotope stack time series covering the last millennium up to the year 2011 CE,
which is built from an extended North Greenland Traverse (NGT) firn core
dataset. This project includes the compilation and the processing of all relevant data
as well as all analyses related to the new NGT-2012 stack, including the
comparison with other climate datasets for the Arctic region.

This software (v0.1) and the included analyses are the basis for the manuscript
"Exceptional temperatures in central-north Greenland ice cores" by M. Hörhold et
al., which is under review for publication in _Nature_.

## Installation

- Download or clone the repository to your machine into a directory of your choice.
- Have a look at the `dependencies.R` within the project's main folder to
   install any of the required R packages not yet present on your system.

## Running analyses

For working with **exNGT** you need to source the initialization script
`init.R`, located within the project's main folder, by calling
`source("init.R")` in R. This loads all relevant R packages and the **exNGT**
function library in `lib/`, thereby preparing everything to start any of the
analyses provided under `analysis/`.
