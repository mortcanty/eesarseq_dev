# application_ful.py
# widget interface for SAR sequential change detection, full scale version
# with classification of time series for crop masking

import ee
ee.Initialize

from eesar.eesarseq_class import assemble_and_run, minimum

import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, gamma, f, chi2
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
aoi = None

geolocator = Nominatim(timeout=10,user_agent='application_full.ipynb')

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
                    .filterDate('2021-09-01', '2022-12-31') \
                    .select('label').mosaic()

def maskDynamicWorld(image): 
    return image.eq(ee.Image.constant(6)) # Built 

#watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)
watermask = dyn.gt(ee.Image.constant(0))

def getS2():
    try:
        S2 = ee.ImageCollection('COPERNICUS/S2_SR') \
                          .filterBounds(aoi) \
                          .filterDate(ee.Date(w_startdate.value),ee.Date(w_enddate.value)) \
                          .sort('CLOUDY_PIXEL_PERCENTAGE',True) \
                          .filter(ee.Filter.contains(rightValue=aoi,leftField='.geo')) \
                          .first() \
                          .clip(aoi)
        S2_rgb = ee.Image.rgb(S2.select('B4'),S2.select('B3'),S2.select('B2')) 
        timestamp = S2.get('system:time_start').getInfo() 
        timestamp = time.gmtime(int(timestamp)/1000)
        timestamp = time.strftime('%x', timestamp).replace('/','')
        timestampS2 = 'S2: 20'+timestamp[4:]+timestamp[0:4]
        return (S2_rgb, timestampS2) 
    except:
        with w_out:
            print('Error: %s'%'S2 not found')
                         

w_location = widgets.Text(
    layout = widgets.Layout(width='150px'),
    value='Houston',
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
    options=['Bitemporal','First','Last','Frequency','ATSF','DynWorld','CropsMask','S2'],
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
w_maskcrops = widgets.Checkbox(
    value=False,
    description='CropsMask',
    disabled=True
)
w_dw = widgets.Checkbox(
    value=False,
    description='DW Mask',
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
w_exportmaskname = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='<path>',
    placeholder=' ',
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
w_plot = widgets.Button(description='PlotAsset',disabled=False)
w_classify = widgets.Button(description='ClassifyCrop',disabled=True)

w_masks = widgets.VBox([w_maskchange,w_maskwater,w_maskcrops,w_dw,w_quick])
w_dates = widgets.VBox([w_startdate,w_enddate])
w_assets = widgets.VBox([w_review,w_plot,w_reset])
w_bmap = widgets.VBox([w_interval,w_maxfreq,w_minfreq])
w_export = widgets.VBox([widgets.HBox([w_export_ass,w_exportassetsname]),
                         widgets.HBox([w_classify,w_exportmaskname])])
w_signif = widgets.VBox([w_significance,w_median])

#Assemble the interface

row1 = widgets.HBox([w_platform,w_orbitpass,w_relativeorbitnumber,w_dates])
row2 = widgets.HBox([w_collect,w_signif,w_opacity,w_export])
row3 = widgets.HBox([w_preview,w_changemap,w_bmap,w_masks,w_assets])
row4 = widgets.HBox([w_out,w_goto,w_location])

box = widgets.VBox([row1,row2,row3,row4])

#Event handers

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    
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
    global aoi
    coords =  geo_json['geometry']['coordinates']
    if action == 'created':
        aoi = ee.Geometry.Polygon(coords)
        w_preview.disabled = True
        w_export_ass.disabled = True
        w_collect.disabled = False
    elif action == 'deleted':
        aoi = None
        w_collect.disabled = True  
        w_preview.disabled = True    
        w_export_ass.disabled = True    
    
def on_collect_button_clicked(b):
    ''' Collect a time series from the archive '''
    global cmaps, bmaps, atsf, count, crs, time_sequence
    with w_out:
        w_out.clear_output()
        clear_layers()
        print('Running on GEE archive COPERNICUS/S1_GRD')
        #assemble time series and run the algorithm
        cmaps, bmaps, count, rons, collection, atsf, time_sequence = assemble_and_run(aoi, median = w_median.value, 
                                                  significance = w_significance.value, startdate=w_startdate.value, 
                                                  enddate=w_enddate.value, platform=w_platform.value, 
                                                  orbitpass=w_orbitpass.value, ron=w_relativeorbitnumber.value)
        crs = ee.Image(collection.first()).select(0).projection().crs().getInfo()        
        w_preview.disabled = False
        w_export_ass.disabled = False
        w_classify.disabled = False       
        if len(rons)>0:
            #display S1 mosaic 
            print('Shortest orbit path series length: %i images\n please wait for raster overlay ...'%count)
            clear_layers()
            S1 = collection.mean()
            m.add_layer(TileLayer(url=GetTileLayerUrl(S1.select(0).visualize(min=-15, max=4)),name='S1'))            
            # scale linear time sequence for classifier and rename bands
            bns = ['band%i'%i for i in range(2*count)]
            time_sequence = time_sequence.rename(bns)
            
            
w_collect.on_click(on_collect_button_clicked)        

def on_classify_button_clicked(b): 
    global crops
    with w_out:
        w_out.clear_output()
        print('Classifying time sequence of %i polarimetric bands ...'%time_sequence.bandNames().size().getInfo())  
        try:
            # landcover band with just 2 classes: crops, etc= 1, else 0
            classValues = [0,1,2,3,4,5,6,7,8]            
            remapValues = [0,0,0,0,1,0,0,0,0]
            label = 'lc'
            lc = dyn.remap(classValues,remapValues).rename(label).toByte()
            sample = time_sequence.multiply(10000).addBands(lc).clip(aoi).stratifiedSample(
                          numPoints = 10000,
                          classBand = label,
                          region = aoi,
                          tileScale = 4,
                          scale = 30 ) 
            # add a random value field to the sample and use it to approximately split 80%
            # of the features into a training set and 20% into a validation set.
            sample = sample.randomColumn()            
            trainingSample = sample.filter('random <= 0.8')
            validationSample = sample.filter('random > 0.8')            
            # train a classifier
            classifier = ee.Classifier.smileNaiveBayes()
            trainedClassifier = classifier.train(
                      features = trainingSample,
                      classProperty = label,
                      inputProperties = time_sequence.bandNames()
            )
            # get information about the trained classifier.
            print('Results of trained classifier', trainedClassifier.explain().getInfo())   
            # get a confusion matrix and overall accuracy for the training sample.
            trainAccuracy = trainedClassifier.confusionMatrix()
            print('Training error matrix', trainAccuracy.getInfo())
            print('Training overall accuracy', trainAccuracy.accuracy().getInfo())  
            
            # get a confusion matrix and overall accuracy for the validation sample.
            validationSample = validationSample.classify(trainedClassifier)
            validationAccuracy = validationSample.errorMatrix(label, 'classification')
            print('Validation error matrix', validationAccuracy.getInfo())
            print('Validation accuracy', validationAccuracy.accuracy().getInfo())

            
            
            
            # classify the time sequence
            crops = time_sequence.multiply(10000).clip(aoi).classify(trainedClassifier)  
            # minimum contiguous area requirement (4.5 hectare)    
            crops = crops.int().selfMask()                 
            contArea = crops.connectedPixelCount(50).selfMask() 
            crops = contArea.gte(ee.Number(50)).selfMask()
            w_maskcrops.disabled = False 
            assetID = 'projects/ee-mortcanty/assets/'+w_exportmaskname.value
            assexport = ee.batch.Export.image.toAsset(crops,
                                    description='assetExportTask', 
                                    pyramidingPolicy={".default": 'sample'},
                                    assetId=assetID, scale=30, maxPixels=1e10)      
            assexport.start()   
            print('Exporting crop map to %s\n task id: %s'%(w_exportmaskname.value,str(assexport.id)))
                
        except Exception as e:
            print('Error: %s'%e)
            
w_classify.on_click(on_classify_button_clicked)           

def on_preview_button_clicked(b):
    ''' Preview change maps '''
    with w_out:  
        try:       
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            bw = 'black,white'
            palette = jet
            rgb = False
            w_out.clear_output()
            print('Shortest orbit path series length: %i images\n previewing please wait for raster overlay ...'%count)
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
                mp = bmaps.select(sel-1)   
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                      
                print('Bitemporal for interval ending: %s'%mp.bandNames().getInfo())
                print('red = positive definite, cyan = negative definite, yellow = indefinite')                 
                palette = rcy
                mx = 3   
            elif w_changemap.value=='ATSF':
                img = ee.Image(atsf).log10().multiply(10).add(15)
                mp = img.select(0)               
                mx = 19 
                palette = bw 
                print( 'ATSF VV' )   
            elif w_changemap.value=='DynWorld':
                mp = dyn.clip(aoi)              
                mx = 8 
                palette = ['419BDF','397D49','88B053','7A87C6','E49635','DFC35A','C4281B','A59B8F','B39FE1']
                print( 'Dynamic World Map' )                   
            elif w_changemap.value=='CropsMask':
                mp = crops = ee.Image('projects/ee-mortcanty/assets/'+w_exportmaskname.value).clip(aoi)
                mx = 1
                palette = bw
                print('Crops')
                w_maskcrops.disabled = False 
            elif w_changemap.value=='S2':
                mp, timestamp = getS2()  
                rgb = True 
                mx = 2000           
                print( timestamp )     
            if not w_quick.value:
                mp = mp.reproject(crs=crs,scale=float(w_exportscale.value))
            if w_maskcrops.value==True:
                crops = ee.Image('projects/ee-mortcanty/assets/'+w_exportmaskname.value)
                mp = mp.updateMask(crops.unmask().neq(1))    
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask) 
            if w_maskchange.value==True:   
                if w_changemap.value=='Frequency':
                    mp = mp.updateMask(mp.gte(w_minfreq.value))
                else:
                    mp = mp.updateMask(mp.gt(0))    
            if rgb: 
                m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=mx,opacity=w_opacity.value)),
                                      name=timestamp))
            else:
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
            aoi = ee.Geometry.Polygon(ee.Geometry(asset.get('system:footprint')).coordinates())
            center = aoi.centroid().coordinates().getInfo()
            center.reverse()
#            m.center = center  
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
    global aoi
    try:       
        assexport = ee.batch.Export.image.toAsset(ee.Image.cat(cmaps,bmaps).byte().clip(aoi),
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
            bmap1 = assetImage.select(ee.List.sequence(3,k+2)).clip(aoi)            
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
#            fn = w_exportassetsname.value.replace('/','-')+'.png'
#            plt.savefig(fn,bbox_inches='tight') 
            w_out.clear_output()
            plt.show()
#            print('Saved to ~/%s'%fn)
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
  
    dc = DrawControl(aoiline={},circlemarker={})
    dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.polygon = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    
    mc = MeasureControl(position='topright',primary_length_unit = 'kilometers')

    dc.on_draw(handle_draw)
    
    lc = LayersControl(position='topright')
    fs = FullScreenControl()
 
    m = Map(center=center, 
                    zoom=6, 
                    layout={'height':'500px','width':'800px'},
                    layers=(ewi,ews,osm),
                    controls=(dc,mc,lc,fs))
    with w_out:
        w_out.clear_output()
        print('Algorithm output') 
    display(m) 
    
    return box     
 