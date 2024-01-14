# Sequential Change Detection in Sentinel-1 Imagery on the Google Earth Engine

Source files for the Docker image mort/eesar_devdocker,
with Jupyter widget interfaces for 
sequential SAR change detection and forest mapping on Google Earth Engine.

Author: mortcanty, August, 2022

Pull and/or run the container with 

    docker run -d -p 8888:8888 --name=eesar mort/eesar_devdocker  

Point your browser to http\://localhost\:8888 to see the Jupyter notebook home page. 
 
Open the Change Detection Notebook 

    interface.ipynb 

and run the first cell for authentication. This only has to be done once in a while.    

Stop the container with

    docker stop eesar 
     
Re-start with

    docker start eesar    
