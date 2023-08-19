# s2cell.py
# widget interface for scanning s2geometry cells

import ee
ee.Initialize
import time, csv, numbers
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from operator import add
from scipy.stats import norm, gamma, f, chi2
import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        MeasureControl,
                        FullScreenControl,
                        Polygon,
                        Popup,
                        basemaps,basemap_to_tiles,
                        LayersControl)
from geopy.geocoders import Nominatim

#  Houston AOI
coords = [[[-96.06878489327627, 29.823701939611126],
           [-96.06878489327627, 29.1423007492751],
           [-95.00860422921377, 29.1423007492751],
           [-95.00860422921377, 29.823701939611126]]]
poly_houston = ee.Geometry.Polygon(coords)

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
            .filterDate('2021-09-01', '2022-03-30') \
            .select('label').mosaic() 

def maskDynamicWorld(image): 
    classes = {'built': 6, 'trees': 1, 'grass': 2, 'crops': 4, 'bare': 7, 'scrub': 5}
    include = [classes[v] for v in w_include.value]
    result = ee.Image.constant(0)
    for i in include:
        result = result.Or(image.eq(ee.Image.constant(i)))     
    return result

geolocator = Nominatim(timeout=10,user_agent='s2cell.ipynb')

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
                       .filter(ee.Filter.contains(rightValue=point,leftField='.geo')) 
                       
            cellid = s2cell.first().getInfo()['id']            
            cellcoords = (point.getInfo())['coordinates']
                       
            coordinates =  s2cell.geometry().geometries().getInfo()[0]['coordinates']            
            poly = ee.Geometry.Polygon(coordinates)            
            k, fmap = viewS2cell('projects/sentinel-change-detection/assets/houston_aoi', poly) # 0.001
            fmap = fmap.updateMask(fmap.gt(0))           
            name = 's2cell'+str(len(m.layers)-3)
            
            locations = [tuple(list(reversed(i))) for i in coordinates[0]]      
            layer = Polygon(locations=locations,
                        color="white", 
                        fill_color = 'black', 
                        name = name)   
            
            cellid = s2cell.first().getInfo()['id']            
            cellcoords = (point.getInfo())['coordinates'][::-1]
            message = widgets.HTML()
            message.value = 's2cell: '+ cellid[-4:]
            popup = Popup(location = cellcoords, child = message, name = 'cellid')
            m.add_layer(popup)
            m.add_layer(layer)
            
            m.add_layer(TileLayer(url=GetTileLayerUrl(fmap.visualize(min=0, max=k/4, 
                                                                     palette='black,blue,cyan,yellow,red')),
                                                                     name='fmap'+cellid[-4:]))           
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
                  .divide(ee.Image.constant(pixels_in_cell)) \
                  .reduceRegion(ee.Reducer.sum(),scale=10,maxPixels=10e10) 
                  
        return ee.List(plots.add(res))
    def mapped(x):
        if isinstance(x,numbers.Number):
            return x
        return 0
    with w_out:
        try:                
            collection = ee.ImageCollection(imgcoll) \
                           .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) \
                           .sort('system:time_start')
            tsl = get_timestamp_list(collection)
            bmap = collection.toBands() \
                             .clip(poly)
            pixels_in_cell = bmap.unmask().select(0) \
                                         .reduceRegion(ee.Reducer.count(),scale=10,maxPixels=10e10) \
                                         .values() \
                                         .getInfo()                                                                       
            bmap = bmap.mask(bmap.mask().And(maskDynamicWorld(dyn)))   
            
            fmap = bmap.multiply(0).where(bmap.gte(1),1).reduce(ee.Reducer.sum())      
                     
            k = bmap.bandNames().length().getInfo()                 
            plots = ee.List(ee.List([1,2,3]).iterate(plot_iter,ee.List([]))).getInfo()           
            x = range(1,k+1)  
            _ = plt.figure(figsize=(10,4),)
            posdef = list(map(mapped, plots[0].values()))
            negdef = list(map(mapped, plots[1].values()))
            indef = list(map(mapped, plots[2].values()))           
            alldef = list(map(add,posdef,negdef))
            alldef = list(map(add,alldef,indef))
            plt.plot(x,posdef,'r-',label='posdef')
            plt.plot(x,negdef,'c-',label='negdef')
            plt.plot(x,indef,'y-',label='indef') 
            plt.plot(x,alldef,'ko-',label='all')        
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
            plt.ylim([0,w_yaxis.value])
            plt.legend()
            plt.show()
            return k, fmap
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
    value='2020-01-01',
    placeholder=' ',
    description='StartDate:',
    disabled=False
)
w_enddate = widgets.Text(
    layout = widgets.Layout(width='200px'),
    value='2022-12-31',
    placeholder=' ',
    description='EndDate:',
    disabled=False
)
w_yaxis = widgets.BoundedFloatText(
    layout = widgets.Layout(width='200px'),
    value='0.2',
    min=0.05,
    max=0.8,
    step=0.05,
    description='y-axis:',
    disabled=False
)
w_include = widgets.SelectMultiple(
    layout = widgets.Layout(height='120px'),
    options=['built','trees','grass','bare','scrub','crops'],
    value=['built'],
    description='Include class',
    disabled=False
)
w_out = widgets.Output(
    layout=widgets.Layout(width='700px',border='1px solid black')
)

w_goto = widgets.Button(description='GoTo',disabled=False)
w_reset = widgets.Button(description='Reset',disabled=False)
w_scan = widgets.Button(description='Scan All',disabled=False)
w_dates = widgets.HBox([w_startdate,w_enddate])    
w_auxil = widgets.VBox([w_include,w_yaxis])
row1 = widgets.HBox([w_dates,w_reset,w_scan])
row2 = widgets.HBox([w_out,w_auxil])
box = widgets.VBox([row1,row2])        

def clear_layers():
    for i in range(20,3,-1): 
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

def on_scan_button_clicked(b):
        
    with w_out:
        try:
            clear_layers()
            w_out.clear_output()
            print('Scanning ...')
            # get the collection of bitemporal change images
            # filter, recast to a single image, rename bands with start times
            # and update masks
            collection = ee.ImageCollection('projects/sentinel-change-detection/assets/houston_aoi') \
                           .filterDate(ee.Date(w_startdate.value), ee.Date(w_enddate.value)) \
                           .sort('system:time_start')
            tsl = get_timestamp_list(collection)   
            bmap = collection.toBands()                        
            bmap = bmap.mask(bmap.mask().And(maskDynamicWorld(dyn)))  
            bmap = bmap.rename(tsl)            
            # get the S2 cells within the Houston AOI        
            s2cells = ee.FeatureCollection('users/djq/bbox_s2_l14') \
                     .filter(ee.Filter.isContained(rightValue=poly_houston,leftField='.geo'))
            print('Number of s2 cells contained in AOI: %s'%s2cells.size().getInfo())
            # recast to a list of cells and iterate over them to create csv
            # cells = s2cells.toList(50000).getInfo()
            cells = s2cells.toList(50000).slice(10000,10100).getInfo()           
            with open('s2cell_features.csv', mode='w') as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                fields = ['s2cell id', 'total image pixel count']
                for item in tsl:
                    fields.append(item)            
                writer.writerow(fields)
                for cell in tqdm(cells):    
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
                    cpp = list(changes_per_period.values())   
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
                    scroll_wheel_zoom=True,
                    layout={'height':'500px','width':'900px'},
                    layers=(ewi,ews,osm),
                    controls=(dc,lc,mc,fs))
    palette = ['419BDF','397D49','88B053','7A87C6','E49635','DFC35A','C4281B','A59B8F','B39FE1']
    m.add_layer( TileLayer(url=GetTileLayerUrl(dyn.visualize(min=0,max=8,
                                  palette=palette)),name= 'Dynamic World') )
    with w_out:
        w_out.clear_output()
        
    display(m) 
    return box      
