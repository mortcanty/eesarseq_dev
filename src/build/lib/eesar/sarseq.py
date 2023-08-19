# sarseq.py
# command line interface to sequential SAR change detection
# export bitemporal maps (bmap) to ImageCollection
# mort canty June, 2022

import sys
import getopt
import time
import math
from ast import literal_eval
import ee
from eesarseq import change_maps

#   --------------------------
#    Auxiliary functions
#   -------------------------

def multipoly2polylist(multipoly):
    ''' Convert a multipolygon to a list of polygons'''
    def fetchpoly(listelem):
        return ee.Geometry.Polygon(multipoly.coordinates().get(listelem))    
    size = multipoly.coordinates().size()
    polylist = ee.List.sequence(0, size.add(-1), 1)
    return polylist.map(fetchpoly)

def getS1collection(t1,  t2, aoi, platform = 'A'):
    ''' Return the longest sequence of S1 images within interval [t1,t2] which completely
        overlap aoi and have a unique relative orbit number and node ''' 
    try :
        collection =  ee.ImageCollection('COPERNICUS/S1_GRD') \
                          .filterBounds(aoi) \
                          .filterDate(ee.Date(t1), ee.Date(t2)) \
                          .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) \
                          .filter(ee.Filter.eq('resolution_meters', 10)) \
                          .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                          .filter(ee.Filter.contains(rightValue=aoi,leftField='.geo'))  
        if platform != 'BOTH':
            collection = collection.filter(ee.Filter.eq('platform_number', platform))                    
        count = collection.size().getInfo()  
        if count < 2: 
            raise ValueError('Less than 2 images found')
        collectionD = collection.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        countD = collectionD.size().getInfo()
        collectionA = collection.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))         
        countA = collectionA.size().getInfo() 
        if countA > countD:
            collection = collectionA 
            node = 'ASCENDING' 
            count = countA
        else:
            collection = collectionD 
            node = 'DESCENDING'      
            count = countD               
        relativeorbitnumbers = map(int,ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo())
        rons = list(set(relativeorbitnumbers))
        ron = rons[0]        
        if len(rons) == 2:    
            collection0 = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(rons[0])))      
            count0 = collection0.size().getInfo() 
            collection1 = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(rons[1])))      
            count1 = collection1.size().getInfo() 
            if count1 > count0:
                ron = rons[1]
                count = count1
                collection = collection1
            else:
                count = count0
                collection = collection0
        return (collection.sort('system:time_start'), count, node, ron)   
    except Exception as e:
        print('Error: %s'%e)                                                               

def get_vvvh(image):
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values '''
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()

def convert_timestamp_list(tsl):
    ''' Make timestamps in YYYYMMDD format '''            
    tsl= [x.replace('/','') for x in tsl]
    tsl = ['T20'+x[4:]+x[0:4] for x in tsl]         
    return tsl

def clipList(current,prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    aoi = ee.Dictionary(prev).get('aoi') 
    ctr = ee.Number(ee.Dictionary(prev).get('ctr'))   
    imlist =  imlist.add(ee.Image(current).multiply(4.4).clip(aoi))
    return ee.Dictionary({'imlist':imlist,'aoi':aoi,'ctr':ctr.add(1)})  

def assemble_and_run(t1, t2, aoi, platform, alpha):
    ''' Collect a time sequence and return the change maps '''
    try:    
        # Gather the time series
        collection, count, node, ron = getS1collection(t1, t2, aoi, platform)          
        print('Images found: %i \nNode: %s \nRelative orbit: %i \nPlatform: %s'%(count, node, ron, platform))
        acquisition_times = ee.List(collection.aggregate_array('system:time_start'))
        pList = collection.map(get_vvvh).toList(500)      
        first = ee.Dictionary({'imlist':ee.List([]),'aoi':aoi,'ctr':ee.Number(0)})          
        imList = ee.List(ee.Dictionary(pList.iterate(clipList,first)).get('imlist'))
        #Run the algorithm 
        return (change_maps(imList, median = True, alpha = alpha), acquisition_times)
    except Exception as e:
        print('Error: %s'%e) 
         
def export_as_image(result, assetpath, acquisition_times, aoi):
    ''' Export change maps as a single multi-band image '''
    smap = ee.Image(result.get('smap')).byte()
    cmap = ee.Image(result.get('cmap')).byte()
    fmap = ee.Image(result.get('fmap')).byte() 
    bmap = ee.Image(result.get('bmap')).byte() 
    timestamplist = []
    for timestamp in acquisition_times.getInfo():
        tmp = time.gmtime(int(timestamp)/1000)
        timestamplist.append(time.strftime('%x', tmp))         
    timestamplist = convert_timestamp_list(timestamplist) 
#   in case of duplicates                  
    timestamplist1 = [timestamplist[i] + '_' + str(i+1) for i in range(len(timestamplist))]             
    cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist1[1:])  
    assexport = ee.batch.Export.image.toAsset(cmaps.byte().clip(aoi),
                                description='assetExportTask', 
                                pyramidingPolicy={".default": 'sample'},
                                assetId=assetpath,scale=10,maxPixels=1e11)      
    assexport.start()
    print('Exporting change maps to %s\ntask id: %s'%(assetpath,str(assexport.id)))         
    
def export_as_image_collection(result, assetpath, acquisition_times, aoi):
    ''' Export the bitemporal change maps as images with properties 'system:time_start'
        and 'system:time_end' to an existing ImageCollection '''
    def tag_images(i):
        i = ee.Number(i)
        image = bmap.select(i)
        image = image.set('system:time_start', acquisition_times.get(i)) \
                     .set('system:time_end', acquisition_times.get(i.add(1))) \
                     .set('system:image_id', ee.String('bmap')) 
        return image       
    try:
        bmap = ee.Image(result.get('bmap')).byte()  
        size = bmap.bandNames().size()
        images = ee.List.sequence(0, size.add(-1))
        tagged_images = ee.List(images.map(tag_images))
        tagged_images = tagged_images
        for i in range(size.getInfo()):
            assetId = assetpath+'/bmap'+str(i)
            assexport = ee.batch.Export.image.toAsset(ee.Image(tagged_images.get(i)).clip(aoi),
                                 description='assetExportTask', 
                                 pyramidingPolicy={".default": 'sample'},
                                 assetId=assetId,scale=10,maxPixels=1e11)   
            assexport.start()
            print('Exporting change maps to %s task id: %s'%(assetId,str(assexport.id)))        
    except Exception as e:
        print('Error: %s'%e)  

#   -------------------
#   main routine
#   -------------------
def main():
    
    ee.Initialize()
     
    usage = '''Usage:
    ---------------------------------------------
    Perform sequential SAR change detection 
    
    python %s  [OPTIONS] assetname
    
    Options:
      -h                  this help
      -c <list>           polygon coordinates of area of interest (default: Houston AOI)
      -t <list>           time interval (default: ['2020-02-01','2022-12-31'])
      -a <float>          false positive rate (default: 0.01)
      -p <int>            platform ('A', 'B' or 'BOTH' default: 'A')
      -d <boolean>        If set: export as ee.ImageCollection, else: export as ee.Image (default: False)
    
    If -d is not set:
        For a sequence of k observations, change maps are exported to GEE assets as k+3-band image:
        band 1: cmap  interval of last recorded change (0 ... k-1) byte
        band 2: smap  interval of first recorded change (0 ... k-1) byte
        band 3: fmap  number of changes in all (0, ... k-1) byte
        bands 4 through k+3:  bitemporal change maps (labeled as end of interval time)  
                              consisting of loewner order of changes in each interval (0, 1, 2, 3) byte  
    Else:
        Export the bitemporal change maps as images with properties 'system:time_start' and 'system:time_end'
        to an existing ImageCollection      
    '''%sys.argv[0]

    options, args = getopt.getopt(sys.argv[1:], 'hc:t:a:p:d')

   
#  Houston AOI
    coords = [[[-96.06878489327627, 29.823701939611126],
               [-96.06878489327627, 29.1423007492751],
               [-95.00860422921377, 29.1423007492751],
               [-95.00860422921377, 29.823701939611126]]]
# #  Juelich AOI    
#     coords = [[[6.348381042480468, 50.88527552329112],
#                [6.479530334472656, 50.88527552329112],
#                [6.479530334472656, 50.94696387166774],
#                [6.348381042480468, 50.94696387166774],
#                [6.348381042480468, 50.88527552329112]]]
    t1, t2 =  ['2020-02-01','2022-12-31']
    alpha = 0.001
    platform = 'A'
    as_collection = False
    
    for option, value in options: 
        if option == '-h':
            print(usage)
            sys.exit(0)
        elif option == '-c':
            coords = literal_eval(value)
        elif option == '-t':
            t1, t2 = literal_eval(value)          
        elif option == '-a':
            alpha = literal_eval(value)  
        elif option == '-p':
            platform = value
        elif option == '-d':
            as_collection = True
                   
#   assetpath = 'projects/sentinel-change-detection/assets/test/' + args[0] 
    assetpath = args[0]
    
    aoi = ee.Geometry.Polygon(coords)
       
    print('-------------------------------')
    print('Sequential SAR Change Detection')
    print('-------------------------------')
    
    result, acquisition_times = assemble_and_run(t1, t2, aoi, platform, alpha)
    
    if as_collection:
        export_as_image_collection(result, assetpath, acquisition_times, aoi)
    else:
        export_as_image(result, assetpath, acquisition_times, aoi)
    
    
if __name__ == '__main__': 
    main()
            
    