# Dockerfile with Jupyter widget interface for 
# SAR Change detection
 
FROM     debian:buster

MAINTAINER Mort Canty "mort.canty@gmail.com"

ENV     REFRESHED_AT 2023-08-18
 
RUN apt-get update && apt-get upgrade

RUN apt-get install -y \
    python3 \
    build-essential \
    pandoc \
    python3-dev \
    python3-pygments \
    libssl-dev \
    libffi-dev \
    libpq-dev \
    gdal-bin \
    libgdal-dev \
    python3-gdal \
    python3-pip


# jupyter notebook with leaflet
RUN     pip3 install --upgrade pip
RUN     pip3 install numpy
RUN     pip3 install scipy
RUN     pip3 install jupyter
RUN     pip3 install ipyleaflet

RUN     pip3 install -U cython
RUN     pip3 install matplotlib==3.0.3

EXPOSE 8888

# setup for earthengine
RUN     pip3 install earthengine-api --upgrade

ENV     REFRESHED_AT 2022-08-28
# install python modules 
copy    setup.py  /home/setup.py
ADD     eesar /home/eesar
WORKDIR /home
RUN     python3 setup.py install  

#COPY    forest_cover.ipynb /home/forest_cover.ipynb
COPY    interface.ipynb /home/interface.ipynb

RUN     pip3 install geopy

# ipython notebook startup script
COPY    notebook.sh /
RUN     chmod u+x /notebook.sh

WORKDIR /home  
CMD     ["/notebook.sh"]
