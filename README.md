# Global, Regional, and National Economic Value of Reducing Amenable Tuberculosis Mortality

This repository contains the source code for the modelling study: **"Global, regional, and national economic value of reducing amenable tuberculosis mortality"**, submitted for publication.

The analysis quantifies the amenable tuberculosis (TB) mortality burden and its associated Value of Lost Welfare (VLW) using data from the Global Burden of Disease (GBD) Study 2021. The model employs a hybrid Monte Carlo simulation framework to account for uncertainty in mortality rates, case fatality rates, and economic parameters.

## 1. Methodology Overview

The code implements a stochastic model with the following key steps:

* **Uncertainty Simulation**: Uses a hybrid Monte Carlo approach (Gamma distribution for raw rates, Triangular distribution for derived parameters) to propagate uncertainty.
* **Amenable Mortality**: Estimated using a regional "Benchmark Case Fatality Rate (CFR)" method.
* **Economic Valuation**: Calculates VSLY (Value of Statistical Life Year) adjusted for income elasticity based on GDP per capita.
* **Aggregation**: Synthesizes global estimates using the Delta Method.

## 2. Repository Structure

* **main_analysis.R**: The core script for simulation, estimation, and economic valuation.
* **visualization.R**: Generates the maps, bar charts, and scatter plots presented in the manuscript.

## 3. Data Availability

The primary data used in this study were obtained from the Global Burden of Disease (GBD) Study 2021 and the World Bank. To reproduce the analysis, please download the required datasets and organize them as follows:

1. Create a folder named **"data"** in your project directory.
2. Place the following files into the **"data"** folder:
    * **death1.csv**: TB Mortality Estimates (Source: GBD 2021 Results Tool)
    * **risk_deathall.csv**: Mortality Attributable to Risk Factors (Source: GBD 2021 Results Tool)
    * **incidence1.csv**: TB Incidence Estimates (Source: GBD 2021 Results Tool)
    * **population1.csv**: Population Estimates (Source: GBD 2021 Results Tool)
    * **le.xlsx**: Life Expectancy Tables (Source: GBD 2021 Reference)
    * **gdp_data.csv**: GDP per capita (Source: World Bank Open Data)
    * **shapefiles/**: Subfolder containing world map shapefiles for visualization

## 4. Requirements

* **R Version**: 4.4.2
* **Dependencies**: tidyverse, readxl, triangle, yll, sf, and ggplot2.

## 5. Usage

1. **Setup Data**: Ensure all datasets are placed in the "data" folder.
2. **Run Analysis**: Execute `main_analysis.R` to generate results in the "results" folder.
3. **Visualize**: Execute `visualization.R` to generate figures in the "figures" folder.

## 6. License

This project is licensed under the MIT License.




