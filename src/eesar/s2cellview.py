
import sys
import getopt
import time
import math
import numpy as np
import matplotlib.pyplot as plt    
from ast import literal_eval 

import ee    
ee.Initialize() 

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
            .filterDate('2021-09-01', '2022-03-30') \
            .select('label').mosaic()
watermask = dyn.gt(ee.Image.constant(0))        

def fetchS2CellGeometry(index):
    s2cells = ee.FeatureCollection('users/djq/bbox_s2_l14').geometry()
    result = ee.Geometry(s2cells.geometries().get(index))
    return result

# stone creek aoi
# def fetchS2CellGeometry(index):
#     coords = [
#           [
#             [
#               -95.71924209594727,
#               29.50860362786118
#             ],
#             [
#               -95.71971416473389,
#               29.506325347350128
#             ],
#             [
#               -95.71756839752197,
#               29.50535426138495
#             ],
#             [
#               -95.71898460388184,
#               29.502440947606576
#             ],
#             [
#               -95.71443557739258,
#               29.500498691853778
#             ],
#             [
#               -95.71263313293457,
#               29.50468196564529
#             ],
#             [
#               -95.7149076461792,
#               29.508192794181117
#             ],
#             [
#               -95.71924209594727,
#               29.50860362786118
#             ]
#           ]
#         ]
#     return ee.Geometry.Polygon(coords)

def viewS2cell(imgcoll, cell_index, t1 = '2016-01-01', t2 = '2023-01-01'):         
    ''' View change fractions in an s2geometry cell '''     
    def plot_iter(current,prev):
        current = ee.Image.constant(current)
        plots = ee.List(prev) 
        res = bmap.multiply(0) \
                  .where(bmap.eq(current),1) \
                  .reduceRegion(ee.Reducer.mean(),scale=10,maxPixels=10e10)
        return ee.List(plots.add(res))
    try:         
        print('Change fraction plots ...') 
        poly = fetchS2CellGeometry(cell_index)     
        bmap = ee.ImageCollection(imgcoll) \
                 .filterDate(ee.Date(t1), ee.Date(t2)) \
                 .filterBounds(poly) \
                 .filter(ee.Filter.contains(rightValue=poly,leftField='.geo')) \
                 .toBands() \
                 .clip(poly) \
                 .updateMask(watermask)      
        k = bmap.bandNames().length().getInfo()                 
        plots = ee.List(ee.List([1,2,3]).iterate(plot_iter,ee.List([]))).getInfo()           
        bns = np.array(list([s[3:9] for s in list(plots[0].keys())])) 
        x = range(1,k+1)  
        _ = plt.figure(figsize=(10,5))
        plt.plot(x,list(plots[0].values()),'ro-',label='posdef')
        plt.plot(x,list(plots[1].values()),'co-',label='negdef')
        plt.plot(x,list(plots[2].values()),'yo-',label='indef')        
        ticks = range(0,k+2)
        labels = [str(i) for i in range(0,k+2)]
        labels[0] = ' '
        labels[-1] = ' '
#        labels[1:-1] = bns 
        if k>50:
            for i in range(1,k+1,2):
                labels[i] = ''
        plt.xticks(ticks,labels,rotation=90)
        plt.legend()
        plt.show()
#            print('Saved to ~/%s'%fn)
    except Exception as e:
        print('Error: %s'%e)   

def main():
    
    
     
    usage = '''Usage:
    ---------------------------------------------
    Filter and display an s2Geometry cell 
    
    python %s  [OPTIONS] assetname
    
    Options:
      -h                  this help
      -i <int>            s2cell index
      -t <list>           time interval (default: ['2016-01-01','2023-01-01'])
    '''%sys.argv[0]
    

    options, args = getopt.getopt(sys.argv[1:], 'hi:t')
    
    t1, t2 =  ['2016-01-01','2023-01-01']
    cell_index = 1
    
    for option, value in options: 
        if option == '-h':
            print(usage)
            sys.exit(0)
        elif option == '-i':
            cell_index = literal_eval(value)
        elif option == '-t':
            t1, t2 = literal_eval(value)        
            
#    assetname = 'projects/sentinel-change-detection/assets/test/' + args[0]
    assetname = args[0]
    
    viewS2cell(assetname, cell_index, t1, t2)
    
if __name__ == '__main__': 
    main()    