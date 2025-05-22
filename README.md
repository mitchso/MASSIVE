Made in Python 3.13.0
# Description
Scripts for visualizing and quantifying nucleic acid MALDI experiments

# Installation
1. Clone the repo
2. Create and activate a new virtual environment, ideally with python 3.13.0
```
python3 -m venv [name]
source [name]/bin/activate
```
3. Install dependencies
```
python -m pip install -r requirements.txt
```

# Data requirements
Individual samples can be analyzed by providing a mass spectrum 

# Usage
Analysis is mainly performed using three objects, Experiment, Sample and Oligo.
Supplemental plotting functions can be accessed via the plotting module.
See notebooks contained in /demos for usage examples.

