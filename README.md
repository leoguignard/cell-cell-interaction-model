# Cell Cell Interaction Model

This repository contains the notebook and the source code that was used to run cell-cell interaction model described in [Contact-dependent cell communications drive morphological invariance during ascidian embryogenesis](https://www.biorxiv.org/content/early/2017/12/24/238741)

**It is important to notice that the notebook as it is proposed here will not work without the datasets that should be manually included in the folder Data**

The datasets are not included in this repository for space reasons. They can be accessed and downloaded, together with the static version of this repository, at http://www.crbm.cnrs.fr/Astec/Datasets/Model-cell-cell-interaction.zip.

## I - CONTENTS OF THE  REPOSITORY
  - Model-run.ipynb : the notebook which shows how to run the model 
  - Data: folder place-holder where the necessary datas to run the model should be placed
  - LICENCE: the licence terms you accept by using the model source code and the notebook
  - README.md: this file
  - loading_data.py: Python script that allows to load the output datas of the segmentation 
  - useful_functions.py: Set of useful functions that are used in the making of the model
  - ClassModel.py: Python class for the cell-cell interaction model

## II - INSTALLATION AND SOFTWARE REQUIREMENTS
In order to be able to run the notebook and the model you need to install several packages:
  * python 2.7 or higher  
    - Installation : Should be installed, use the command line `python -V` to get the installed version or visit https://www.python.org/
the following libraries are necessary:
  * pip, an installer for python (https://pypi.python.org/pypi/pip)
    - Installation : run command line `sudo apt-get install python-pip python-dev build-essential` 
  * numpy,scipy,matplotlib, different scientific packages for python  (http://www.numpy.org, http://www.scipy.org, http://matplotlib.org)
    - Installation : run command line `sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose`
  * jupyter notebook, a open-source framework to share scientific code (http://jupyter.org)
    - Installation : run command line `sudo pip install jupyter`
  * Routines for plotting area-weighted two- and three-circle venn diagrams  
    - Installation : run command line `sudo pip install matplotlib-venn`

## II.iii - TROUBELSHOOTING
  - The installation of jupyter may failed according of the version of your ipython : to correct this issue run `sudo pip install IPython==5.0 `
  - The notebook execution may failed with this error : AttributeError: 'AxesSubplot' object has no attribute 'violinplot'
you just have to update matplotlib running the command `sudo pip install matplotlib --upgrade`
  - The notebook execution may failed due to pandas version you just have to run the command `sudo pip install pandas --upgrade` and `sudo pip install natsort --upgrade`



## III - RUNNING THE NOTEBOOK
The code is written in python and can easly be run with a web interface using a Notebook (http://jupyter.org) 
You can find a user interface tutorial here (http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb)
To start the notebook, open a terminal and  run the command `jupyter notebook path/to/Package/Model-run.ipynb`


