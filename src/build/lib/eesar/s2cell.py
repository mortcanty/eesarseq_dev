# forest.py
# widget interface for SAR forest mapping with omnibus change detection

import ee
ee.Initialize
import time, csv
from tqdm.notebook import tqdm
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

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
            .filterDate('2021-09-01', '2022-03-30') \
            .select('label').mosaic()
watermask = dyn.gt(ee.Image.constant(0)) 

def maskDynamicWorld(image): 
    return image.eq(ee.Image.constant(6)) # Built 

#  Houston AOI
coords = [[[-96.06878489327627, 29.823701939611126],
           [-96.06878489327627, 29.1423007492751],
           [-95.00860422921377, 29.1423007492751],
           [-95.00860422921377, 29.823701939611126]]]
poly_houston = ee.Geometry.Polygon(coords)

geolocator = Nominatim(timeout=10,user_agent='s2cell.ipynb')
        
watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)
settlement = ee.Image("DLR/WSF/WSF2015/v1")

def GetTileLayerUrl(image):
    map_id = ee.Image(image).getMapId()
    return map_id["tile_fetcher"].url_format        

def handle_draw(self, action, geo_json):
    coords =  geo_json['geometry']['coordinates']  
    if action == 'created':
        with w_out:
            point = ee.Geometry.Point(coords)
            s2cell = ee.FeatureCollection('users/djq/bbox_s2_l14') \
                       .filterBounds(poly_houston) \
                       .filter(ee.Filter.contains(rightValue=point,leftField='.geo')) \
                       .geometry()
            coordinates =  s2cell.geometries().getInfo()[0]['coordinates']            
            poly = ee.Geometry.Polygon(coordinates)            
            viewS2cell('users/djq/demo', poly)           
            name = 's2cell'+str(len(m.layers)-3)
            
            locations = [tuple(list(reversed(i))) for i in coordinates[0]]      
            layer = Polygon(locations=locations,
                        color="white", 
                        fill_color = 'black', 
                        name = name)   
            m.add_layer(layer)
    elif action == 'deleted':
        point = None

def get_timestamp_list(collection):
    ''' make timestamps from image collection in YYYYMMDD format '''   
    acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()
    tsl = []
    for timestamp in acquisition_times:
        tmp = time.gmtime(int(timestamp)/1000)
        tsl.append(time.strftime('%x', tmp))        
    tsl= [x.replace('/','') for x in tsl]  
    tsl = ['20'+x[4:]+x[0:4] for x in tsl]         
    return tsl

def viewS2cell(imgcoll, poly):         
    ''' View change fractions in an s2geometry cell '''     
    def plot_iter(current,prev):
        current = ee.Image.constant(current)
        plots = ee.List(prev) 
        res = bmap.multiply(0) \
                  .where(bmap.eq(current),1) \
                  .reduceRegion(ee.Reducer.mean(),scale=10,maxPixels=10e10)
        return ee.List(plots.add(res))
    with w_out:
        try:          
            collection = ee.ImageCollection(imgcoll) \
                           .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) 
            tsl = get_timestamp_list(collection)
            bmap = collection.toBands() \
                             .clip(poly) \
                             .updateMask(watermask)                        
            bmap = bmap.mask(bmap.mask().And(maskDynamicWorld(dyn)))                     
                     
            k = bmap.bandNames().length().getInfo()                 
            plots = ee.List(ee.List([1,2,3]).iterate(plot_iter,ee.List([]))).getInfo()           
            x = range(1,k+1)  
            _ = plt.figure(figsize=(10,4),)
            plt.plot(x,list(plots[0].values()),'ro-',label='posdef')
            plt.plot(x,list(plots[1].values()),'co-',label='negdef')
            plt.plot(x,list(plots[2].values()),'yo-',label='indef')        
            ticks = range(0,k+2)
            labels = [str(i) for i in range(0,k+2)]        
            labels[0] = ' '
            labels[-1] = ' '
            labels[1:-1] = tsl 
            if k > 60:
                for i in range(1,k+1):
                    if i % 4 != 0:
                        labels[i] = ''
            plt.xticks(ticks,labels,rotation=90)
            plt.ylim([0,0.4])
            plt.legend()
            plt.show()
    #            print('Saved to ~/%s'%fn)
        except Exception as e:
            print('Error: %s'%e)   
            
# ********************
# The widget interface
# ********************
w_location = widgets.Text(
    layout = widgets.Layout(width='150px'),
    value='Houston',
    placeholder=' ',
    description='',
    disabled=False
)

w_changemapassetsname = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='users/djq/demo',
    placeholder=' ',
    disabled=False
)    

w_startdate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2018-01-01',
    placeholder=' ',
    description='StartDate:',
    disabled=False
)
w_enddate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2018-04-01',
    placeholder=' ',
    description='EndDate:',
    disabled=False
)

w_out = widgets.Output(
    layout=widgets.Layout(width='700px',border='1px solid black')
)

w_goto = widgets.Button(description='GoTo',disabled=False)
w_reset = widgets.Button(description='Reset',disabled=False)
w_scan = widgets.Button(description='Scan All',disabled=False)
w_dates = widgets.HBox([w_startdate,w_enddate])    
row1 = widgets.HBox([w_dates,w_reset,w_goto,w_location])
row2 = widgets.HBox([w_out,w_scan])
box = widgets.VBox([row1,row2])

def clear_layers():
    if len(m.layers)>9:
        m.remove_layer(m.layers[9])    
    if len(m.layers)>8:
        m.remove_layer(m.layers[8])
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

def on_scan_button_clicked(b):
        
    with w_out:
        try:
            clear_layers()
            w_out.clear_output()
            print('Scanning ...')
            # get the collection of bitemporal change images
            # filter, recast to a single image, rename bands with start times
            # and update masks
            collection = ee.ImageCollection('users/djq/demo') \
                           .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value))
            tsl = get_timestamp_list(collection)   
            bmap = collection.toBands().updateMask(watermask)                        
            bmap = bmap.mask(bmap.mask().And(maskDynamicWorld(dyn)))  
            bmap = bmap.rename(tsl)            
            # get the S2 cells within the Houston AOI        
            s2cells = ee.FeatureCollection('users/djq/bbox_s2_l14') \
                     .filter(ee.Filter.isContained(rightValue=poly_houston,leftField='.geo'))
            print('Number of s2 cells contained in AOI: %s'%s2cells.size().getInfo())
            # recast to a list of cells and iterate over them to create csv
            cells = s2cells.toList(50000).getInfo()           
            with open('s2cell_features.csv', mode='w') as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                fields = ['s2cell id', 'total image pixel count']
                for item in tsl:
                    fields.append(item)            
                writer.writerow(fields)
                for cell in tqdm(cells[:100]):    
                    cell_id = cell['id']       
                    cell_coords = cell['geometry']['coordinates']
                    poly = ee.Geometry.Polygon(cell_coords)      
                    pixels_in_cell = bmap.unmask().select(0) \
                                         .clip(poly) \
                                         .reduceRegion(ee.Reducer.count(),scale=10,maxPixels=10e10) \
                                         .getInfo()[tsl[0]]        
                    changes_per_period = bmap.multiply(0) \
                                             .where(bmap.gte(1),1) \
                                             .clip(poly) \
                                             .reduceRegion(ee.Reducer.sum(),scale=10,maxPixels=10e10) \
                                             .getInfo() 
                    cpp= list(changes_per_period.values())   
                    row = [cell_id, pixels_in_cell]
                    for item in cpp:
                        row.append(int(item))                                                     
                    writer.writerow(row)
            print('done')               
        except Exception as e:
            print('Error: %s'%e)

w_scan.on_click(on_scan_button_clicked)        

def on_goto_button_clicked(b):
    try:
        location = geolocator.geocode(w_location.value)
        m.center = (location.latitude,location.longitude)
        m.zoom = 11
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_goto.on_click(on_goto_button_clicked)
            
def run():
    global m, center
    center = poly_houston.centroid().coordinates().getInfo()[::-1]
    osm = basemap_to_tiles(basemaps.OpenStreetMap.Mapnik)
    ews = basemap_to_tiles(basemaps.Esri.WorldStreetMap)
    ewi = basemap_to_tiles(basemaps.Esri.WorldImagery)
    
    mc = MeasureControl(position='topright',primary_length_unit = 'kilometers')
    
    dc = DrawControl(polyline={})
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
        
    display(m) 
    return box      
