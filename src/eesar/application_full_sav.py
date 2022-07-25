# application_ful.py
# widget interface for SAR sequential change detection, full scale version

import ee
ee.Initialize

from eesar.sarseqalgorithm import change_maps 

import time
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, gamma, f, chi2
from collections import Counter
import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        MeasureControl,
                        FullScreenControl,
                        basemaps,basemap_to_tiles,
                        LayersControl)
from geopy.geocoders import Nominatim

'''
 ********************
 The widget interface
 ********************
'''
poly = None

geolocator = Nominatim(timeout=10,user_agent='application_full.ipynb')

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
                    .filterDate('2021-09-01', '2022-03-30') \
                    .select('label').mosaic()

def maskDynamicWorld(image): 
    return image.eq(ee.Image.constant(6)) # Built 

watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)

w_location = widgets.Text(
    layout = widgets.Layout(width='150px'),
    value='Jülich',
    placeholder=' ',
    description='',
    disabled=False
)
w_orbitpass = widgets.RadioButtons(
    layout = widgets.Layout(width='200px'),
    options=['ASCENDING','DESCENDING'],
    value='ASCENDING',
    description='Pass:',
    disabled=False
)
w_changemap = widgets.RadioButtons(
    options=['Bitemporal','First','Last','Frequency'],
    value='First',
    layout = widgets.Layout(width='200px'),
    disabled=False
)
w_interval = widgets.BoundedIntText(
    min=1,
    value=1,
    layout = widgets.Layout(width='150px'),
    description='Interval:',
    disabled=True
)
w_maxfreq = widgets.BoundedIntText(
    min=1,
    value=20,
    layout = widgets.Layout(width='150px'),
    description='MaxFreq:',
    disabled=True
)
w_minfreq = widgets.BoundedIntText(
    min=1,
    value=1,
    layout = widgets.Layout(width='150px'),
    description='MinFreq:',
    disabled=True
)
w_platform = widgets.RadioButtons(
    layout = widgets.Layout(width='200px'),
    options=['Both','A','B'],
     value='A',
    description='Platform:',
    disabled=False
)
w_relativeorbitnumber = widgets.IntText(
    value='0',
    layout = widgets.Layout(width='150px'),
    description='RelOrbit:',
    disabled=False
)
w_exportassetsname = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='projects/<user>/assets/<path>',
    placeholder=' ',
    disabled=False
)
w_exportdrivename = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='<path>',
    placeholder=' ',
    disabled=False
)
w_exportscale = widgets.FloatText(
    value=10,
    placeholder=' ',
    description='Scale ',
    disabled=False
)
w_startdate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2019-01-01',
    placeholder=' ',
    description='StartDate:',
    disabled=False
)
w_enddate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2020-01-01',
    placeholder=' ',
    description='EndDate:',
    disabled=False
)
w_median = widgets.Checkbox(
    layout = widgets.Layout(width='200px'),
    value=True,
    description='MedianFilter',
    disabled=False
)
w_quick = widgets.Checkbox(
    value=True,
    description='QuickPreview',
    disabled=False
)
w_significance = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='0.01',
    min=0.0001,
    max=0.05,
    step=0.001,
    description='Signif:',
    disabled=False
)
w_maskchange = widgets.Checkbox(
    value=True,
    description='NCMask',
    disabled=False
)
w_maskwater = widgets.Checkbox(
    value=True,
    description='WaterMask',
    disabled=False
)
w_opacity = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='1.0',
    min=0.0,
    max=1.0,
    step=0.1,
    description='Opacity:',
    disabled=False
)
w_dw = widgets.Checkbox(
    value=False,
    description='DW Mask',
    disabled=False
)
w_out = widgets.Output(
    layout=widgets.Layout(width='700px',border='1px solid black')
)

w_collect = widgets.Button(description="Collect",disabled=True)
w_preview = widgets.Button(description="Preview",disabled=True)
w_reset = widgets.Button(description='Reset',disabled=False)
w_review = widgets.Button(description="ReviewAsset",disabled=False)
w_goto = widgets.Button(description='GoTo',disabled=False)
w_export_ass = widgets.Button(description='ExportToAssets',disabled=True)
w_export_drv = widgets.Button(description='ExportToDrive',disabled=True)
w_plot = widgets.Button(description='PlotAsset',disabled=False)

w_masks = widgets.VBox([w_maskchange,w_maskwater,w_dw,w_quick])
w_dates = widgets.VBox([w_startdate,w_enddate])
w_assets = widgets.VBox([w_review,w_plot,w_reset])
w_bmap = widgets.VBox([w_interval,w_maxfreq,w_minfreq])
w_export = widgets.VBox([widgets.HBox([w_export_ass,w_exportassetsname]),
                         widgets.HBox([w_export_drv,w_exportdrivename])])
w_signif = widgets.VBox([w_significance,w_median])


#Assemble the interface
row1 = widgets.HBox([w_platform,w_orbitpass,w_relativeorbitnumber,w_dates])
row2 = widgets.HBox([w_collect,w_signif,w_opacity,w_export])
row3 = widgets.HBox([w_preview,w_changemap,w_bmap,w_masks,w_assets])
row4 = widgets.HBox([w_out,w_goto,w_location])

box = widgets.VBox([row1,row2,row3,row4])

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    w_export_drv.disabled = True
    
def on_changemap_widget_change(b):   
    if b['new']=='Bitemporal':
        w_interval.disabled=False
    else:
        w_interval.disabled=True    
    if b['new']=='Frequency':
        w_maxfreq.disabled=False
        w_minfreq.disabled=False
    else:
        w_maxfreq.disabled=True 
        w_minfreq.disabled=True   

#These widget changes require a new collect
w_orbitpass.observe(on_widget_change,names='value')
w_platform.observe(on_widget_change,names='value')
w_relativeorbitnumber.observe(on_widget_change,names='value')
w_startdate.observe(on_widget_change,names='value')
w_enddate.observe(on_widget_change,names='value')
w_median.observe(on_widget_change,names='value')
w_significance.observe(on_widget_change,names='value')
w_changemap.observe(on_changemap_widget_change,names='value')  

        
def clear_layers():
    for i in range(20,2,-1): 
        if len(m.layers)>i:
            m.remove_layer(m.layers[i])    


def on_reset_button_clicked(b):
    try:
        clear_layers()
        w_out.clear_output()
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_reset.on_click(on_reset_button_clicked)         

def on_goto_button_clicked(b):
    try:
        location = geolocator.geocode(w_location.value)
        m.center = (location.latitude,location.longitude)
        m.zoom = 11
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_goto.on_click(on_goto_button_clicked)


def GetTileLayerUrl(image):
    map_id = ee.Image(image).getMapId()
    return map_id["tile_fetcher"].url_format        

def handle_draw(self, action, geo_json):
    global poly
    coords =  geo_json['geometry']['coordinates']
    if action == 'created':
        poly = ee.Geometry.Polygon(coords)
        w_preview.disabled = True
        w_export_ass.disabled = True
        w_export_drv.disabled = True 
        w_collect.disabled = False
    elif action == 'deleted':
        poly = None
        w_collect.disabled = True  
        w_preview.disabled = True    
        w_export_ass.disabled = True
        w_export_drv.disabled = True      

def getS1collection():
    return ee.ImageCollection('COPERNICUS/S1_GRD') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) \
                      .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) \
                      .filter(ee.Filter.eq('resolution_meters', 10)) \
                      .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                      .filter(ee.Filter.eq('orbitProperties_pass', w_orbitpass.value)) 
            
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
    
def on_collect_button_clicked(b):
    ''' Collect a time series from the archive '''
    global cmaps, bmaps, count, archive_crs
    
    def minimum(a, b):      
        if a <= b:
            return a
        else:
            return b
    
    def trim_list(current): 
        current = ee.Image(current)
        bns = current.bandNames().slice(0,count-1)
        return current.select(ee.List.sequence(0,count-2),bns)
    
    with w_out:
        try:
            w_out.clear_output()
            clear_layers()
            print('running on GEE archive COPERNICUS/S1_GRD')
            
            collection = getS1collection()          
            if w_platform.value != 'Both':
                collection = collection.filter(ee.Filter.eq('platform_number', w_platform.value))
            if w_relativeorbitnumber.value > 0:
                collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(w_relativeorbitnumber.value)))    
            count = collection.size().getInfo()     
            if count==0:
                raise ValueError('No images found') 
            print('Images found: %i, platform: %s'%(count,w_platform.value))
            
            collection = collection.sort('system:time_start')                                                                      
            archive_crs = ee.Image(collection.first()).select(0).projection().crs().getInfo()
            
            relativeorbitnumbers = map( int, ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo() )
            rons = list(set(relativeorbitnumbers))
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
                
                path_lengths = [ctr[uts] for uts in uniquetimestamps]
                
                # only process paths with unambiguous time series
                if len(set(path_lengths)) == 1:
                
                    mosaic_lengths = [ctr[timestamp] for timestamp in uniquetimestamps]
                 
                    cList = collection_ron.map(get_vvvh).toList(500)   
          
                    # make list of combined (mosaicked) images along orbit path 
                    mLen = ee.List(mosaic_lengths) 
                    first = ee.Dictionary({'plist': ee.List([]), 'clist': cList})  
                    pList = ee.List(ee.Dictionary(mLen.iterate(make_mosaics, first)).get('plist'))
     
                    first = ee.Dictionary({'imlist':ee.List([]),'enl':ee.Number(4.4),'poly':poly})        
                    imList = ee.List(ee.Dictionary(pList.iterate(clipList, first)).get('imlist'))   
                    
                    count = minimum(imList.size().getInfo(),count) 
                   
                    #Run the algorithm ************************************************
                    result = change_maps(imList, w_median.value, w_significance.value)
                    #******************************************************************
                    smap = ee.Image(result.get('smap')).byte()
                    cmap = ee.Image(result.get('cmap')).byte()
                    fmap = ee.Image(result.get('fmap')).byte() 
                    bmap = ee.Image(result.get('bmap')).byte()              
                    cmap_list = cmap_list.add( ee.Image.cat(cmap,smap,fmap).rename(['cmap','smap','fmap']) ) 
                    bmap_list = bmap_list.add( bmap ) 
                else:
                    print('relative orbit %i excluded: includes differing series lengths'%ron)
                    
            # mosaic smap,cmap,fmap images   
            cmaps = ee.ImageCollection(cmap_list).mosaic()
            # truncate bitemporal maps to length of shortest series
            bmap_list = bmap_list.map(trim_list)
            # mosaic bitemporal maps over orbit paths
            bmaps = ee.ImageCollection(bmap_list).mosaic().rename(uniquetimestamps[1:count])
             
            w_preview.disabled = False
            w_export_ass.disabled = False
            w_export_drv.disabled = False
            #Display S1 mosaic 
            if len(rons)>0:
                print( 'please wait for raster overlay ...' )
                clear_layers()
                S1 = collection.mosaic().select(0).visualize(min=-15, max=4)
                m.add_layer(TileLayer(url=GetTileLayerUrl(S1),name='S1'))                          
        except Exception as e:
            print('Error: %s'%e) 

w_collect.on_click(on_collect_button_clicked)                  

def on_preview_button_clicked(b):
    ''' Preview change maps '''
    with w_out:  
        try:       
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            palette = jet
            w_out.clear_output()
            print('Shortest orbit path series length: %i images, previewing (please wait for raster overlay) ...'%count)
            if w_changemap.value=='First':
                mp = ee.Image(cmaps.select('smap')).byte()
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))
                
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp=ee.Image(cmaps.select('cmap')).byte()
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = ee.Image(cmaps.select('fmap')).byte()
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))
                mx = w_maxfreq.value
                print('Change frequency :\n blue = few, red = many')
            elif w_changemap.value == 'Bitemporal':
                sel = int(w_interval.value)
                sel = min(sel,count-1)
                sel = max(sel,1)                               
                print('Bitemporal for interval ending: %s'%mp.bandNames().getInfo())
                print('red = positive definite, cyan = negative definite, yellow = indefinite')  
                mp = bmaps.select(sel-1)
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))
                palette = rcy
                mx = 3     
            if not w_quick.value:
                mp = mp.reproject(crs=archive_crs,scale=float(w_exportscale.value))
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask) 
            if w_maskchange.value==True:   
                if w_changemap.value=='Frequency':
                    mp = mp.updateMask(mp.gte(w_minfreq.value))
                else:
                    mp = mp.updateMask(mp.gt(0))    
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx,
                                                                   opacity=w_opacity.value, 
                                                                   palette=palette)),name=w_changemap.value))           
        except Exception as e:
            print('Error: %s'%e)

w_preview.on_click(on_preview_button_clicked)      

def on_review_button_clicked(b):
    ''' Examine change maps exported to user's assets ''' 
    with w_out:  
        try: 
#          test for existence of asset                  
            _ = ee.Image(w_exportassetsname.value).getInfo()
#          ---------------------------            
            asset = ee.Image(w_exportassetsname.value)
            poly = ee.Geometry.Polygon(ee.Geometry(asset.get('system:footprint')).coordinates())
            center = poly.centroid().coordinates().getInfo()
            center.reverse()
            m.center = center  
            bnames = asset.bandNames().getInfo()[3:]
            count = len(bnames)               
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            smap = asset.select('smap').byte()
            cmap = asset.select('cmap').byte()
            fmap = asset.select('fmap').byte()
            bmap = asset.select(list(range(3,count+3)),bnames).byte()      
            palette = jet
            w_out.clear_output()
            print('Series length: %i images, reviewing (please wait for raster overlay) ...'%(count+1))
            if w_changemap.value=='First':
                mp = smap
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp = cmap
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = fmap
                mx = w_maxfreq.value
                print('Change frequency :\n blue = few, red = many')
            elif w_changemap.value == 'Bitemporal':
                sel = int(w_interval.value)
                sel = min(sel,count-1)
                sel = max(sel,1)
                print('Bitemporal: interval %i'%w_interval.value)
                print('red = positive definite, cyan = negative definite, yellow = indefinite')  
                mp = bmap.select(sel-1)
                palette = rcy
                mx = 3     
            if len(m.layers)>6:
                m.remove_layer(m.layers[6])   
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask)
            if w_maskchange.value==True:    
                if w_changemap.value=='Frequency':
                    mp = mp.updateMask(mp.gte(w_minfreq.value))
                else:
                    mp = mp.updateMask(mp.gt(0))    
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx, opacity=w_opacity.value, palette=palette)),name=w_changemap.value))
        except Exception as e:
            print('Error: %s'%e)
    
w_review.on_click(on_review_button_clicked)   

def on_export_ass_button_clicked(b):
    ''' Export to assets '''
    try:         
        assexport = ee.batch.Export.image.toAsset(ee.Image.cat(cmaps,bmaps).byte().clip(poly),
                                    description='assetExportTask', 
                                    pyramidingPolicy={".default": 'sample'},
                                    assetId=w_exportassetsname.value,scale=10,maxPixels=1e10)      
        assexport.start()
        with w_out: 
            w_out.clear_output() 
            print('Exporting change maps to %s\n task id: %s'%(w_exportassetsname.value,str(assexport.id)))
    except Exception as e:
        with w_out:
            print('Error: %s'%e)                                          
    
w_export_ass.on_click(on_export_ass_button_clicked)  

def on_export_drv_button_clicked(b):
    ''' Export to Google Drive '''
    try:     
        cmaps = ee.Image.cat(cmaps,bmaps)  
        fileNamePrefix=w_exportdrivename.value.replace('/','-')            
        gdexport = ee.batch.Export.image.toDrive(ee.Image.cat(cmaps,bmaps).byte().clip(poly),
                                    description='driveExportTask', 
                                    folder = 'gee',
                                    fileNamePrefix=fileNamePrefix,scale=10,maxPixels=1e10)   
        gdexport.start()
        with w_out:
            w_out.clear_output()
            print('Exporting change maps to Drive/gee/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id))) 
    except Exception as e:
        with w_out:
            print('Error: %s'%e) 

w_export_drv.on_click(on_export_drv_button_clicked)            

def on_plot_button_clicked(b):          
    ''' plot change fractions form asset '''       
    global bmap1 
    def plot_iter(current,prev):
        current = ee.Image.constant(current)
        plots = ee.List(prev) 
        res = bmap1.multiply(0) \
                  .where(bmap1.eq(current),1) \
                  .reduceRegion(ee.Reducer.mean(),scale=10,maxPixels=10e10)
        return ee.List(plots.add(res))
    with w_out:
        try:
            w_out.clear_output()            
            print('Change fraction plots ...')                  
            assetImage = ee.Image(w_exportassetsname.value)
            k = assetImage.bandNames().length().subtract(3).getInfo()            
            bmap1 = assetImage.select(ee.List.sequence(3,k+2)).clip(poly)            
            if w_maskwater.value:
                bmap1 = bmap1.updateMask(watermask) 
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
            labels[1:-1] = bns 
            if k>50:
                for i in range(1,k+1,2):
                    labels[i] = ''
            plt.xticks(ticks,labels,rotation=90)
            plt.legend()
            fn = w_exportassetsname.value.replace('/','-')+'.png'
            plt.savefig(fn,bbox_inches='tight') 
            w_out.clear_output()
            plt.show()
            print('Saved to ~/%s'%fn)
        except Exception as e:
            print('Error: %s'%e)               
    
w_plot.on_click(on_plot_button_clicked)

def run():
    ''' Run the interface '''
    global m, center
    center = [51.0,6.4]
    osm = basemap_to_tiles(basemaps.OpenStreetMap.Mapnik)
    ews = basemap_to_tiles(basemaps.Esri.WorldStreetMap)
    ewi = basemap_to_tiles(basemaps.Esri.WorldImagery)
  
    dc = DrawControl(polyline={},circlemarker={})
    dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.polygon = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}

    dc.on_draw(handle_draw)
    
    lc = LayersControl(position='topright')
    fs = FullScreenControl()
 
    m = Map(center=center, 
                    zoom=6, 
                    layout={'height':'500px','width':'800px'},
                    layers=(ewi,ews,osm),
                    controls=(dc,lc,fs))
    with w_out:
        w_out.clear_output()
        print('Algorithm output') 
    display(m) 
    
    return box     
 