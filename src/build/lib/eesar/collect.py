# collect.py 
# Utilities to assemble time series and 
# run the sequential change detectio algorithm

import time
import math
from collections import Counter
from eesar.sarseqalgorithm import change_maps 

import ee

def getS1collection(poly, orbitpass, startdate, enddate):
    return ee.ImageCollection('COPERNICUS/S1_GRD') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(startdate), ee.Date(enddate)) \
                      .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) \
                      .filter(ee.Filter.eq('resolution_meters', 10)) \
                      .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                      .filter(ee.Filter.eq('orbitProperties_pass', orbitpass)) 
                      
def get_vvvh(image):   
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values '''
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()                          

def minimum(a, b):      
    if a <= b:
        return a
    else:
        return b

def get_timestamp_list(collection):
    ''' make timestamps from image collection in YYYYMMDD format '''   
    acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()
    tsl = []
    for timestamp in acquisition_times:
        tmp = time.gmtime(int(timestamp)/1000)
        tsl.append(time.strftime('%x', tmp))        
    tsl= [x.replace('/','') for x in tsl]  
    tsl = ['T20'+x[4:]+x[0:4] for x in tsl]         
    return tsl


def clipList(current,prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    poly = ee.Dictionary(prev).get('poly') 
    enl = ee.Number(ee.Dictionary(prev).get('enl')) 
    imlist = imlist.add(ee.Image(current).multiply(enl).clip(poly))
    return ee.Dictionary({'imlist':imlist,'poly':poly,'enl':enl})
        
        
def make_mosaics(current, prev):
    ''' return equitemporal mosaicked images in plist '''
    mLen = ee.List(current)
    prev = ee.Dictionary(prev)
    pList = ee.List(prev.get('plist'))
    cList = ee.List(prev.get('clist'))   
    pList = pList.add( ee.ImageCollection(cList.slice(0,mLen)).mosaic() )    
    return ee.Dictionary({'plist':pList, 'clist':cList.slice(mLen)})

def assemble_and_run(poly, median = True, significance = 0.01, startdate='20180101', 
                     enddate='20190101', platform='A', orbitpass='DESCENDING', relativeorbitnumber=0):
    
    ''' Collect a time series from the archive involving more than one orbit path'''    
    
    def trim_list(current): 
        current = ee.Image(current)
        bns = current.bandNames().slice(0,count-1)
        return current.select(ee.List.sequence(0,count-2),bns)
    
    try:
        
        collection = getS1collection(poly, orbitpass, startdate, enddate)          
        if platform != 'Both':
            collection = collection.filter(ee.Filter.eq('platform_number', platform))
        if relativeorbitnumber > 0:
            collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(relativeorbitnumber)))    
        count = collection.size().getInfo()     
        if count==0:
            raise ValueError('No images found') 
        print('Images found: %i, platform: %s'%(count,platform))
        
        collection = collection.sort('system:time_start')                                                                      
        archive_crs = ee.Image(collection.first()).select(0).projection().crs().getInfo()
        
        rons = map( int, ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo() )
        rons = list(set(rons))
        rons.sort()
        print('Relative orbit numbers: %s'%str(rons))
        
        cmap_list = ee.List([])
        bmap_list = ee.List([])
        
        count = 500            
        for ron in rons:
           
            collection_ron = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', ron))
            
            timestamplist = get_timestamp_list(collection_ron)    

            ctr = Counter(timestamplist)    
            uniquetimestamps = list(set(timestamplist))
            uniquetimestamps.sort() 
            
            orbit_lengths = [ctr[timestamp] for timestamp in uniquetimestamps]
            print('Orbit: %i, lengths %s'%(ron, str(orbit_lengths)))
        
            cList = collection_ron.map(get_vvvh).toList(500)   
 
            # make list of combined (mosaicked) images along orbit path 
            mLen = ee.List(orbit_lengths) 
            first = ee.Dictionary({'plist': ee.List([]), 'clist': cList})  
            pList = ee.List(ee.Dictionary(mLen.iterate(make_mosaics, first)).get('plist'))

            first = ee.Dictionary({'imlist':ee.List([]),'enl':ee.Number(4.4),'poly':poly})        
            imList = ee.List(ee.Dictionary(pList.iterate(clipList, first)).get('imlist'))   
           
            count = minimum(imList.size().getInfo(),count) 
          
            #Run the algorithm ************************************************
            result = change_maps(imList, median, significance)
            #******************************************************************
            smap = ee.Image(result.get('smap')).byte()
            cmap = ee.Image(result.get('cmap')).byte()
            fmap = ee.Image(result.get('fmap')).byte() 
            bmap = ee.Image(result.get('bmap')).byte()              
            cmap_list = cmap_list.add( ee.Image.cat(cmap,smap,fmap).rename(['cmap','smap','fmap']) ) 
            bmap_list = bmap_list.add( bmap ) 
                
        # mosaic smap,cmap,fmap images   
        cmaps = ee.ImageCollection(cmap_list).mosaic()
        # truncate bitemporal maps to length of shortest series
        bmap_list = bmap_list.map(trim_list)
        # mosaic bitemporal maps over orbit paths
        bmaps = ee.ImageCollection(bmap_list).mosaic().rename(uniquetimestamps[1:count])
         
        return (cmaps, bmaps, count, rons, collection, archive_crs)   

                   
    except Exception as e:
        print('Error: %s'%e) 
           
if __name__ == '__main__':
    pass
