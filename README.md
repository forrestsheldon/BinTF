# BinTF
 Normalisation and labeling utilities for transfected cells

 After cloning the repository you will need to install the
 appropriate Python and Julia packages.

 For Python, there is a conda environment file in the main
 directory called environment.yml. With conda installed:

 ```
 conda env create -f environment.yml
 conda activate BinTF
 ```

 should create and activate the environment.

 For Julia, with Julia installed and added to your path type

 ```
 julia
 ```

 to open the REPL, after which

 ```
 using Pkg
 Pkg.activate(".")
 Pkg.instantiate()
 ```

 should install all the required packages. Once that is done,
 you can view the Pluto notebooks in DevelopmentNotebooks/
 with

 ```
 julia

 using Pluto
 Pluto.run()
 ```
 and then selecting the desired notebook. To run the Streamlit
 app, in the main directory run

 ```
 streamlit run BarcodeContApp.py --server.fileWatcherType none
 ```

 Because of the large number of count files, the fileWatcher has
 to be turned off.
