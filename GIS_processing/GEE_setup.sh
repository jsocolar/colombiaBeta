# bash script run via terminal
# Running python 3.8.2 via conda 4.8.3, shipped via Anaconda3-2020.02-MacOSX-x86_64
# to run, paste filepath to terminal. For example: /Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/GEE_setup.sh

conda env remove --name gee_interface

y | conda create --name gee_interface      # Create a conda environment

eval "$(conda shell.bash hook)"

conda activate gee_interface           # Activate the environment

conda install -c conda-forge earthengine-api # Install the Earth Engine Python API

earthengine authenticate          # Authenticate your access with the command line tool

conda install numpy
conda install pandas