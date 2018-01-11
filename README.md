********************************************************************************
**************************************** Model Cell Cell Interaction ***********
********************************************************************************

This file contains instructions to compile and run the model presented in the paper.

I - CONTENTS OF THE  REPOSITORY
-------------------------------
Once uncompressed, the folder contains the following elements:
  - Model-run.ipynb : the notebook which shows how to run the model 
  - Data: folder containing the necessary datas to build the model 
  - license.txt: the licence terms you accept by using the workflow 
  - README.md: this file 
  - loading_data.py: Python script that allows to load the output datas of the segmentation 
  - useful_functions.py: Set of useful functions that are used in the making of the model
  - ClassModel.py: Python class for the cell-cell model

II - INSTALLATION AND SOFTWARE REQUIREMENTS
-------------------------------------------
The Model Cell Cell Interaction can run only on Linux system and was tested on Unbuntu 14.04 64 bits.

In order to be able to compile the different source code you need to install several packages, mostly from the terminal application
to run the different codes:
	* python 2.7 or higher  
		-> Installation : Should be installed, use the command line `python -V` to get the installed version or visit  http://php.net 
the following libraries are necessary:
	* pip, an installer for python (https://pypi.python.org/pypi/pip)
		-> Installation : run command line `sudo apt-get install python-pip python-dev build-essential` 
	* numpy,scipy,matplotlib, different scientific packages for python  (http://www.numpy.org, http://www.scipy.org, http://matplotlib.org)
		->  Installation : run command line `sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose`
	* jupyter notebook, a open-source framework to share scientific code (http://jupyter.org)
		->  Installation : run command line `sudo pip install jupyter`
	* Routines for plotting area-weighted two- and three-circle venn diagrams  
		->  Installation : run command line `sudo pip install matplotlib-venn`


II.iii - TROUBELSHOOTING
-------------------------
   - The installation of jupyter may failed according of the version of your ipython : to correct this issue run `sudo pip install IPython==5.0 `
   - The notebook execution may failed with this error : AttributeError: 'AxesSubplot' object has no attribute 'violinplot'
	you just have to update matplotlib running the command `sudo pip install matplotlib --upgrade`
   - The notebook execution may failed due to pandas version you just have to run the command `sudo pip install pandas --upgrade` and `sudo pip install natsort --upgrade`



 III - RUNNING ASTEC- FIGURES
 ----------------------------
The code is written in python and can easly be run with a web interface using a Notebook (http://jupyter.org)  
You can find a user interface tutorial here (http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb)
To start the notebook, open a terminal and  run the command `jupyter notebook path/to/Package/ASTEC-JupyterNotebook.ipynb`


