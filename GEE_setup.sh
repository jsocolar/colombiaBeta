# bash script run via terminal
# Running python 3.8.2 via conda 4.8.3, shipped via Anaconda3-2020.02-MacOSX-x86_64

conda create --name gee_interface      # Create a conda environment

conda activate gee_interface           # Activate the environment

conda install -c conda-forge earthengine-api # Install the Earth Engine Python API

earthengine authenticate          # Authenticate your access with the command line tool

conda install pandas
conda install numpy