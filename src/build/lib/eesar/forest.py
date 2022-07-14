# forest.py
# widget interface for SAR forest mapping with omnibus change detection

import ee
ee.Initialize
import time, math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, gamma, f, chi2
import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        MeasureControl,
                        FullScreenControl,
                        Polygon,
                        basemaps,basemap_to_tiles,
                        LayersControl)
from geopy.geocoders import Nominatim

# *****************************************
# The sequental change detection algorithm
# *****************************************

def chi2cdf(chi2, df):
    """Calculates Chi square cumulative distribution function for
       df degrees of freedom using the built-in incomplete gamma
       function gammainc().
    """
    return ee.Image(chi2.divide(2)).gammainc(df.divide(2))

def det(im):
    """Calculates the  determinant of 2x2 diagonal covariance matrix."""
    return ee.Image(im).expression('b(0)*b(1)')

def logdet(im):
    """Calculates the log of the determinant of 2x2 diagonal covariance matrix."""
    return ee.Image(im).expression('b(0)*b(1)').log()

def log_det_sum(im_list):
    """Returns log of determinant of the sum of the  images in im_list."""
    im_list = ee.List(im_list)
    suml = ee.ImageCollection(im_list).reduce(ee.Reducer.sum())
    return ee.Image(det(suml)).log()

def omnibus(im_list, k):
    """Calculates the omnibus test statistic, bivariate case.""" 
    m = 4.4
    im_list = ee.List(im_list)
    k2logk = k.multiply(k.log()).multiply(2)
    k2logk = ee.Image.constant(k2logk)
    sumlogdets = ee.ImageCollection(im_list.map(logdet)).reduce(ee.Reducer.sum())
    logdetsum = log_det_sum(im_list)
    return k2logk.add(sumlogdets).subtract(logdetsum.multiply(k)).multiply(-2*m)

def change_map(im_list, median, alpha):
    '''Returns omnibus change map for im_list'''
    k = im_list.length()
    p_value = ee.Image.constant(1).subtract(chi2cdf(omnibus(im_list, k), k.subtract(1).multiply(2)))
    if median:
        p_value = p_value.focal_median(1.5)
    return p_value.multiply(0).where(p_value.lt(alpha), 1)
    
# ********************
# The widget interface
# ********************

poly = None
geolocator = Nominatim(timeout=10,user_agent='tutorial-pt-4.ipynb')

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
    options=['Omnibus','EWM','PALSAR','Hansen'],
    value='Omnibus',
    layout = widgets.Layout(width='150px'),
    disabled=False
)
w_platform = widgets.RadioButtons(
    layout = widgets.Layout(width='200px'),
    options=['Both','A','B'],
     value='Both',
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
    value='forest/',
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
    value='2018-04-01',
    placeholder=' ',
    description='StartDate:',
    disabled=False
)
w_enddate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2018-11-01',
    placeholder=' ',
    description='EndDate:',
    disabled=False
)
w_stride = widgets.BoundedIntText(
    value=1,
    min=1,
    description='Stride:',
    layout = widgets.Layout(width='150px'),
    disabled=False
)
w_median = widgets.Checkbox(
    layout = widgets.Layout(width='200px'),
    value=True,
    description='MedianFilter',
    disabled=False
)
w_significance = widgets.BoundedFloatText(
    layout = widgets.Layout(width='150px'),
    value='0.05',
    min=0.0001,
    max=0.05,
    step=0.001,
    description='Alpha:',
    disabled=False
)
w_maskwater = widgets.Checkbox(
    value=True,
    description='WaterMask',
    disabled=False
)
w_useshape = widgets.Checkbox(
    layout = widgets.Layout(width='180px'),
    value=False,
    description='Shape',
    disabled=False
)
w_county = widgets.BoundedIntText(
    value=45,
    min=1,
    description='County',
    layout = widgets.Layout(width='150px'),
    disabled=False
)
w_settlement = widgets.Checkbox(
    value=True,
    description='SettlementMask',
    disabled=False
)
w_forest = widgets.Checkbox(
    value=False,
    description='Forest bdries',
    disabled=False
)
w_out = widgets.Output(
    layout=widgets.Layout(width='700px',border='1px solid black')
)

w_collect = widgets.Button(description="Collect",disabled=False)
w_preview = widgets.Button(description="Preview",disabled=True,layout = widgets.Layout(width='200px'))
w_review = widgets.Button(description="ReviewAsset",disabled=False)
w_goto = widgets.Button(description='GoTo',disabled=False)
w_export_ass = widgets.Button(description='ExportToAssets',disabled=True)
w_export_drv = widgets.Button(description='ExportToDrive',disabled=True)
w_reset = widgets.Button(description='Reset',disabled=False)

w_masks = widgets.VBox([w_maskwater,w_settlement,w_forest])
w_dates = widgets.VBox([w_startdate,w_enddate])
w_export = widgets.VBox([widgets.HBox([w_useshape,w_county]),
                         widgets.HBox([w_export_ass,w_exportassetsname])])
w_signif = widgets.VBox([w_significance,w_median])

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    
def on_goto_button_clicked(b):
    try:
        location = geolocator.geocode(w_location.value)
        m.center = (location.latitude,location.longitude)
        m.zoom = 11
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_goto.on_click(on_goto_button_clicked)

#These widget changes require a new collect
w_orbitpass.observe(on_widget_change,names='value')
w_platform.observe(on_widget_change,names='value')
w_relativeorbitnumber.observe(on_widget_change,names='value')
w_startdate.observe(on_widget_change,names='value')
w_enddate.observe(on_widget_change,names='value')
w_stride.observe(on_widget_change,names='value')
w_median.observe(on_widget_change,names='value')
w_significance.observe(on_widget_change,names='value') 

row1 = widgets.HBox([w_platform,w_orbitpass,w_relativeorbitnumber,w_dates])
row2 = widgets.HBox([w_collect,w_signif,w_stride,w_export])
row3 = widgets.HBox([w_preview,w_changemap,w_masks,w_review,w_reset])
row4 = widgets.HBox([w_out,w_goto,w_location])

box = widgets.VBox([row1,row2,row3,row4])

#@title Collect

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
    elif action == 'deleted':
        poly = None
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
                      .filter(ee.Filter.eq('orbitProperties_pass', w_orbitpass.value)) \
                      .filter(ee.Filter.contains(rightValue=poly,leftField='.geo'))

def getS2collection():
    return ee.ImageCollection('COPERNICUS/S2_SR') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(w_startdate.value),ee.Date(w_enddate.value)) \
                      .sort('CLOUDY_PIXEL_PERCENTAGE',True) \
                      .filter(ee.Filter.contains(rightValue=poly,leftField='.geo'))             
                      
def getLS8collection():
    return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
                      .filterBounds(poly) \
                      .filterDate(ee.Date(w_startdate.value),ee.Date(w_enddate.value)) \
                      .sort('CLOUD_COVER_LAND',True) \
                      .filter(ee.Filter.contains(rightValue=poly,leftField='.geo'))     
                      
def getPALSAR():
    PALSAR = ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/FNF') \
                  .filterDate('2017-01-01', '2017-12-31')                   
    mp = PALSAR.first().clip(poly)
    prj = mp.projection()
    scale = prj.nominalScale()
    mp0 = mp.multiply(0) 
    mp = mp0.where(mp.eq(1),1).selfMask()          
    return mp.reproject(prj.atScale(scale)) 

def getEWM():
    mp = ee.ImageCollection("ESA/WorldCover/v100").first().clip(poly)             
    prj = mp.projection()
    scale = prj.nominalScale()
    mp0 = mp.multiply(0) 
    mp = mp0.where(mp.eq(10),1).selfMask()          
    return mp.reproject(prj.atScale(scale)) 
                                            
def getHansen(canopy=10,contiguous=9):
    gfc = ee.Image('UMD/hansen/global_forest_change_2020_v1_8')  
    treecover = gfc.select('treecover2000').clip(poly)
    lossyear = gfc.select('lossyear').clip(poly)
    loss = lossyear.where(lossyear.gte(1),1).selfMask()
    forest2000 = treecover.gte(ee.Number(canopy)).selfMask()               
    prj = forest2000.projection()
    scale = prj.nominalScale()
    mp = forest2000.connectedPixelCount(100).gte(ee.Number(contiguous)).selfMask()
    mp = mp.where(lossyear.gte(1),0).selfMask()
    return (mp.reproject(prj.atScale(scale)) , loss.reproject(prj.atScale(scale)))        
                  
def get_vvvh(image):   
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values '''
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()

def convert_timestamp_list(tsl):
    ''' make timestamps in YYYYMMDD format '''           
    tsl= [x.replace('/','') for x in tsl]  
    tsl = ['T20'+x[4:]+x[0:4] for x in tsl]         
    return tsl[::int(w_stride.value)]

def clipList(current,prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    poly = ee.Dictionary(prev).get('poly') 
    enl = ee.Number(ee.Dictionary(prev).get('enl')) 
    ctr = ee.Number(ee.Dictionary(prev).get('ctr'))   
    stride = ee.Number(ee.Dictionary(prev).get('stride'))
    imlist =  ee.Algorithms.If(ctr.mod(stride).eq(0),
        imlist.add(ee.Image(current).multiply(enl).clip(poly)),imlist)
    return ee.Dictionary({'imlist':imlist,'poly':poly,'enl':enl,'ctr':ctr.add(1),'stride':stride})  

def clear_layers():
    if len(m.layers)>7:
        m.remove_layer(m.layers[7])
    if len(m.layers)>6:
        m.remove_layer(m.layers[6])
    if len(m.layers)>5:
        m.remove_layer(m.layers[5])
    if len(m.layers)>4:
        m.remove_layer(m.layers[4])   
    if len(m.layers)>3:
        m.remove_layer(m.layers[3])
        
def on_reset_button_clicked(b):
    try:
        clear_layers()
        w_out.clear_output()
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_reset.on_click(on_reset_button_clicked)        

def on_collect_button_clicked(b):
    ''' Collect an S1 time series from the archive 
    '''
    global result, count, timestamplist, poly, prj
    with w_out:
        try:
            w_out.clear_output()
            clear_layers()        
            if w_useshape.value:
                polys = ee.FeatureCollection('projects/ee-mortcanty/assets/dvg1krs_nw').geometry()
                poly = ee.Geometry(polys.geometries().get(w_county.value))   
            collection = getS1collection()          
            if w_platform.value != 'Both':
                collection = collection.filter(ee.Filter.eq('platform_number', w_platform.value))
            if w_relativeorbitnumber.value > 0:
                collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(w_relativeorbitnumber.value)))    
            count = collection.size().getInfo()     
            if count<2:
                raise ValueError('Less than 2 images found') 
#            footprint = ee.Geometry.Polygon(collection.first().get('system:footprint').getInfo()['coordinates'])
            print('Images found: %i, platform: %s, stride: %i'%(count,w_platform.value,w_stride.value))
            acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()            
            prj = ee.Image(collection.first()).select(0).projection()            
            timestamplist = []
            for timestamp in acquisition_times:
                tmp = time.gmtime(int(timestamp)/1000)
                timestamplist.append(time.strftime('%x', tmp))  
            #Make timestamps in YYYYMMDD format       
            timestamplist = convert_timestamp_list(timestamplist) 
            print(timestamplist)    
            collection = collection.sort('system:time_start') 
            relativeorbitnumbers = map(int,ee.List(collection.aggregate_array('relativeOrbitNumber_start')).getInfo())
            rons = list(set(relativeorbitnumbers))
            print('Relative orbit numbers: '+str(rons))
            if len(rons)==2:
                print('Please select one relative orbit number before preview/export')
            pList = collection.map(get_vvvh).toList(500)      
            first = ee.Dictionary({'imlist':ee.List([]),'enl':ee.Number(4.4),'poly':poly,'ctr':ee.Number(0),'stride':ee.Number(int(w_stride.value))})          
            imList = ee.List(ee.Dictionary(pList.iterate(clipList,first)).get('imlist'))   
            count = imList.size().getInfo()           
            #Get a preview as collection mean                                                                                     
            S1 = collection.mean().select(0).visualize(min=-15,max=4)
            #Run the algorithm ************************************************
            result = change_map(imList, w_median.value, w_significance.value)
            #******************************************************************
            w_preview.disabled = False
            w_export_ass.disabled = False
            w_export_drv.disabled = False
            #Display preview 
            if len(rons)==1:
                print( 'please wait for raster overlay ...' )
                clear_layers()
                m.add_layer(TileLayer(url=GetTileLayerUrl(S1),name='S1'))    
                if w_forest.value:
                    def intersect(feature):
                        return feature.intersection(poly)                
                    geoms = ee.FeatureCollection('projects/ee-mortcanty/assets/veg02_f') \
                                              .filterBounds(poly).map(intersect) \
                                              .geometry().geometries() \
                                              .getInfo()                                         
                    polys = []
                    for dct in geoms:
                        if dct['type'] == 'Polygon':
                            polys.append(dct['coordinates'][0])      
                        for p in polys:
                            for l in p:
                                l.reverse()
                    layer = Polygon(locations=polys,
                                color="white", 
                                fill_color = 'black', 
                                name = 'Forest boundaries')    
                    m.add_layer(layer)    
                if w_useshape.value:
                    locations =  poly.coordinates().getInfo()[0]
                    locations = [tuple(list(reversed(i))) for i in locations]
                    layer = Polygon(locations=locations,
                                color="red", 
                                fill_color = 'black', 
                                name = 'County Nr. '+str(w_county.value))    
                    center = poly.centroid().coordinates().getInfo()
                    center.reverse() 
                    m.center = center              
                    m.add_layer(layer)                     
        except Exception as e:
            print('Error: %s'%e) 

w_collect.on_click(on_collect_button_clicked)                  

watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)
settlement = ee.Image("DLR/WSF/WSF2015/v1")

def on_preview_button_clicked(b):
    ''' Preview change maps
    '''
    global prj
    with w_out:  
        try:       
            ggg = 'black,green'
            yyy = 'black,yellow'
            ccc = 'black,cyan'
            rrr = 'black,red'
            bbb = 'black,blue'
            w_out.clear_output()
            print('Series length: %i images, previewing (please wait for raster overlay) ...'%count)
            if w_changemap.value=='Omnibus':
#              no-change map                
                palette = ggg
                mp = ee.Image(result).int().Not().selfMask()
#              minimum contiguous area requirement (0.5 hectare)                     
                contArea = mp.connectedPixelCount().selfMask()
                mp = contArea.gte(ee.Number(50)).selfMask() 
#              mask settled areas 
                if w_settlement.value:              
                    mp = mp.where(settlement.eq(255),0).selfMask()           
#              forest cover in hectares                
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
#              scale the map so not affected by zoom   
                scale = prj.nominalScale()
                mp = mp.reproject(prj.atScale(scale))   
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = poly,
                                    scale = scale,
                                    maxPixels = 1e13)    
#              dictionary of pixels counts (only on band in this case)                                   
                pixelCount = mp.reduceRegion( ee.Reducer.count(), geometry = poly ,scale = scale, maxPixels = 1e13 )
#              pixel size                
                onePixel = forestCover.getNumber('constant').divide(pixelCount.getNumber('constant'))
                minAreaUsed = onePixel.multiply(50).getInfo()
#              F1 score with Hansen as ground truth
                mph, _ = getHansen()
                TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly).get('constant')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly).get('constant')
                FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),poly).get('constant')
                P = TP.divide(TP.add(FP)).getInfo() 
                R = TP.divide(TP.add(FN)).getInfo()                
                F1 = 2.0*P*R/(P+R)                   
                print('Omnibus change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('constant').getInfo()))
                print('Minimum forest area used (ha) ', minAreaUsed)
                print('F1 score relative to Hansen ', F1)        
            elif w_changemap.value=='EWM':
                mp = getEWM()
                if w_settlement.value:              
                    mp = mp.where(settlement.eq(255),0).selfMask()     
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = poly,
                                    scale = 10,
                                    maxPixels = 1e13)         
                mph, _ = getHansen()
                TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('Map')
                TP = ee.Number(TP)              
                FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('Map')
                FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('Map')
                P = TP.divide(TP.add(FP)).getInfo() 
                R = TP.divide(TP.add(FN)).getInfo()                
                F1 = 2.0*P*R/(P+R)                   
                print('EWM change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('Map').getInfo()))
                print('F1 score relative to Hansen ', F1)                       
            elif w_changemap.value=='PALSAR':
                mp = getPALSAR()
                if w_settlement.value:              
                    mp = mp.where(settlement.eq(255),0).selfMask()     
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = poly,
                                    scale = 30,
                                    maxPixels = 1e13)          
                mph, _ = getHansen()
                TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly).get('fnf')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly).get('fnf')
                FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),poly).get('fnf')
                P = TP.divide(TP.add(FP)).getInfo() 
                R = TP.divide(TP.add(FN)).getInfo()                
                F1 = 2.0*P*R/(P+R)                         
                print('PALSAR change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('fnf').getInfo()))
                print('F1 score relative to Hansen ', F1)   
            elif w_changemap.value=='Hansen':                
                mp, loss = getHansen()       
                palette = ggg
                if w_settlement.value:              
                    mp = mp.where(settlement.eq(255),0).selfMask()        
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = poly,
                                    scale = 30,
                                    maxPixels = 1e13)                
                pixelCount = mp.reduceRegion( ee.Reducer.count(), geometry = poly ,scale = 30, maxPixels = 1e13 )
                onePixel = forestCover.getNumber('treecover2000').divide(pixelCount.getNumber('treecover2000'))
                minAreaUsed = onePixel.multiply(9).getInfo() 
                print('Hansen change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('treecover2000').getInfo()))
                print('Minimum forest area used (ha) ', minAreaUsed)  
                m.add_layer(TileLayer(url=GetTileLayerUrl(loss.visualize(min=0, max=1,  palette=rrr)), name=w_changemap.value+' loss'))            
            if len(m.layers)>8:
                m.remove_layer(m.layers[8])
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask)           
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=1,  palette=palette)), name=w_changemap.value))                
        except Exception as e:
            print('Error: %s'%e)

w_preview.on_click(on_preview_button_clicked)      

def on_review_button_clicked(b):
    ''' Examine change maps exported to user's assets
    ''' 
    global poly
    with w_out:  
        try: 
#           test for existence of asset                  
            _ = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value).getInfo()
#           ---------------------------           
            w_out.clear_output() 
            mp = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value).selfMask()
            poly = ee.Geometry.Polygon(ee.Geometry(mp.get('system:footprint')).coordinates())
            center = poly.centroid().coordinates().getInfo()
            center.reverse()
            m.center = center  
#           forest cover in hectares                
            pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
            forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = poly,
                                    scale = 10,
                                    maxPixels = 1e13)   
            mph, _ = getHansen()
            TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('constant') 
            TP = ee.Number(TP)
            FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('constant')
            FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),poly,maxPixels = 1e13).get('constant')
            P = TP.divide(TP.add(FP)).getInfo() 
            R = TP.divide(TP.add(FN)).getInfo()                
            F1 = 2.0*P*R/(P+R)                      
            print('Omnibus change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('constant').getInfo()))  
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=1, palette='black,green')),name='omnibus')) 
            print('F1 score relative to Hansen ', F1)          
        except Exception as e:
            print('Error: %s'%e)
    
w_review.on_click(on_review_button_clicked)   

def on_export_ass_button_clicked(b):
    ''' Export to assets
    '''
    try:
        mp = ee.Image(result).int()
        mp0 = mp.multiply(0)
        mp = mp0.where(mp.eq(0),1).selfMask()
#      minimum area requirement                      
        mp = mp.connectedPixelCount(150).gte(ee.Number(50)).selfMask()
#      mask settled areas 
        if w_settlement.value:              
            mp = mp.where(settlement.eq(255),0).selfMask()
        assexport = ee.batch.Export.image.toAsset(mp.byte().clip(poly),
                                description='assetExportTask', 
                                pyramidingPolicy={".default": 'sample'},
                                assetId='projects/ee-mortcanty/assets/'+w_exportassetsname.value, 
                                scale=10, maxPixels=1e13)      
        assexport.start()
        with w_out: 
            w_out.clear_output() 
            print('Exporting forest map to %s\n task id: %s'%('projects/ee-mortcanty/assets/'+w_exportassetsname.value, 
                                                              str(assexport.id)))             
    except Exception as e:
        with w_out:
            print('Error: %s'%e)                                          
    
w_export_ass.on_click(on_export_ass_button_clicked)  

#@title Run the interface
def run():
    global m, center
    center = [51.0,6.4]
    osm = basemap_to_tiles(basemaps.OpenStreetMap.Mapnik)
    ews = basemap_to_tiles(basemaps.Esri.WorldStreetMap)
    ewi = basemap_to_tiles(basemaps.Esri.WorldImagery)
    
    mc = MeasureControl(position='topright',primary_length_unit = 'kilometers')
    
    dc = DrawControl(polyline={},circlemarker={})
    dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.polygon = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}

    dc.on_draw(handle_draw)
    
    lc = LayersControl(position='topright')
    fs = FullScreenControl()
 
    m = Map(center=center, 
                    zoom=11, 
                    layout={'height':'500px','width':'800px'},
                    layers=(ewi,ews,osm),
                    controls=(dc,lc,mc,fs))
    with w_out:
        w_out.clear_output()
        print('Algorithm output') 
        
    display(m) 
    return box      