# Seasonal Loops between Local Outgoing Longwave Radiation and Surface Temperature

![OLR Loops](https://i.ibb.co/SdQTbjG/Screenshot-2021-05-02-at-12-54-53.png)

This repository provides code to reproduce the figures and values found in the project of the same name. 

### The interactive OLR loop viewer can be found [here](https://bdgrichards.github.io/OLR-Loop-Viewer/)

The project included both simple data analysis using Python Jupyter Notebooks and extensive offline calculations, also in Python, primarily using the [CLIMLAB](https://pypi.org/project/climlab/) package on a SLURM high performance computing cluster. The core Python scripts used on the cluster have been included in the Cluster Code directory.

The code uses ERA5 and CERES datasets, which are not included here due to their size, but which can be easily obtained from [the CERES data products website](https://ceres.larc.nasa.gov/data/) and [the ERA5 website](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5). The processed datasets, where the raw datasets have been re-gridded, averaged and had loop height calculations performed, are included for reference in the Processed Data Outputs directory . 

## Cluster Code 

Cluster code is provided for recreation of the offline calculations used to investigate the effects of holding atmospheric lapse rates and relative humidity profiles constant. The datasets created with offline calculations are referred to as 'calculated' or 'calc'. These scripts require a parent script to be run fully, as the given Python files must be run for every latitude / longitude grid point on a 1x1 degree grid, for every month between Jan 2001 and Dec 2020. The scripts primarily use ERA5 monthly averaged surface temperature, atmospheric relative humidity profiles and atmospheric temperature profiles which are not included here due to their size but can be easily accessed from the [the ERA5 website](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5).

As the cluster code is run for each individual spacial and temporal location independently, the resulting ~ 65,000 output files (one for each 1x1 grid location) are combined using the Xarray `open_mfdataset` method. This is the purpose of the `data_comb.py` script. 

No loop height calculations are performed on the cluster, these are performed in the data processing notebook.

## Jupyter Notebooks

The `Data processing.ipynb` notebook takes the output from the cluster scripts, as well as raw ERA5 and CERES data, and produces the reduced datasets required for plotting in the `Paper figures.ipynb` notebook. This includes re-gridding datasets that are at different spatial resolutions, removing excess temporal data and renaming variables and coordinates for consistency. 

Loop calculations are also performed in the data processing notebook, as well as linear regression, for all locations on the processed datasets. All of the resulting data are then saved to NetCDF files so that they can be passed to the plotting notebook. 

The plotting notebook is organised by figure number, with each `savefig` corresponding to a given figure. Additional values such as means, biases and mean absolute errors are also calculated.   

***

Feel free to contact me at <bdgrichards@gmail.com>
