# application.py
# widget interface for SAR sequential change detection
# colab version

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
    return ee.Image(chi2.divide(2)).gammainc(ee.Number(df).divide(2))

def det(im):
    """Calculates determinant of 2x2 diagonal covariance matrix."""
    return im.expression('b(0)*b(1)')

def log_det_sum(im_list, j):
    """Returns log of determinant of the sum of the first j images in im_list."""
    im_ist = ee.List(im_list)
    sumj = ee.ImageCollection(im_list.slice(0, j)).reduce(ee.Reducer.sum())
    return ee.Image(det(sumj)).log()

def log_det(im_list, j):
    """Returns log of the determinant of the jth image in im_list."""
    im = ee.Image(ee.List(im_list).get(j.subtract(1)))
    return ee.Image(det(im)).log()

def pval(im_list, j, m=4.4):
    """Calculates -2logRj for im_list and returns P value and -2mlogRj."""
    im_list = ee.List(im_list)
    j = ee.Number(j)
    

    m2logRj = (log_det_sum(im_list, j.subtract(1))
               .multiply(j.subtract(1))
               .add(log_det(im_list, j))
               .add(ee.Number(2).multiply(j).multiply(j.log()))
               .subtract(ee.Number(2).multiply(j.subtract(1))
               .multiply(j.subtract(1).log()))
               .subtract(log_det_sum(im_list,j).multiply(j))
               .multiply(-2).multiply(m))
    
    # correction to simple Wilks approximation
    # pv = ee.Image.constant(1).subtract(chi2cdf(m2logRj, 2))    
    one= ee.Number(1)
    rhoj = one.subtract(one.add(one.divide(j.multiply(j.subtract(one)))).divide(6).divide(m))
    omega2j = one.subtract(one.divide(rhoj)).pow(2.0).divide(-2)    
    rhojm2logRj = m2logRj.multiply(rhoj)
    pv = ee.Image.constant(1).subtract(chi2cdf(rhojm2logRj,2).add(chi2cdf(rhojm2logRj,6).multiply(omega2j)).subtract(chi2cdf(rhojm2logRj,2).multiply(omega2j)))   

    return (pv, m2logRj)

def p_values(im_list,m=4.4):
    """Pre-calculates the P-value array for a list of images."""
    im_list = ee.List(im_list)
    k = im_list.length()

    def ells_map(ell):
        """Arranges calculation of pval for combinations of k and j."""
        ell = ee.Number(ell)
        # Slice the series from k-l+1 to k (image indices start from 0).
        im_list_ell = im_list.slice(k.subtract(ell), k)

        def js_map(j):
            """Applies pval calculation for combinations of k and j."""
            j = ee.Number(j)
            pv1, m2logRj1 = pval(im_list_ell, j)
            return ee.Feature(None, {'pv': pv1, 'm2logRj': m2logRj1})

        # Map over j=2,3,...,l.
        js = ee.List.sequence(2, ell)
        pv_m2logRj = ee.FeatureCollection(js.map(js_map))

        # Calculate m2logQl from collection of m2logRj images.
        m2logQl = ee.ImageCollection(pv_m2logRj.aggregate_array('m2logRj')).sum()

        # correction to simple Wilks approximation
        # pvQl = ee.Image.constant(1).subtract(chi2cdf(m2logQl, ell.subtract(1).multiply(2)))        
        one = ee.Number(1)
        f = ell.subtract(1).multiply(2)
        rho = one.subtract(ell.divide(m).subtract(one.divide(ell.multiply(m))).divide(f).divide(3))
        omega2 = f.multiply(one.subtract(one.divide(rho)).pow(2)).divide(-4)
        rhom2logQl = m2logQl.multiply(rho)   
        pvQl = ee.Image.constant(1).subtract(chi2cdf(rhom2logQl,f).add(chi2cdf(rhom2logQl,f.add(4)).multiply(omega2)).subtract(chi2cdf(rhom2logQl,f).multiply(omega2)))                     
        
        pvs = ee.List(pv_m2logRj.aggregate_array('pv')).add(pvQl)
        return pvs

    # Map over l = k to 2.
    ells = ee.List.sequence(k, 2, -1)
    pv_arr = ells.map(ells_map)

    # Return the P value array ell = k,...,2, j = 2,...,l.
    return pv_arr

def filter_j(current, prev):
    """Calculates change maps; iterates over j indices of pv_arr."""
    pv = ee.Image(current)
    prev = ee.Dictionary(prev)
    pvQ = ee.Image(prev.get('pvQ'))
    i = ee.Number(prev.get('i'))
    cmap = ee.Image(prev.get('cmap'))
    smap = ee.Image(prev.get('smap'))
    fmap = ee.Image(prev.get('fmap'))
    bmap = ee.Image(prev.get('bmap'))
    alpha = ee.Image(prev.get('alpha'))
    j = ee.Number(prev.get('j'))
    cmapj = cmap.multiply(0).add(i.add(j).subtract(1))
    # Check      Rj?            Ql?                  Row i?
    tst = pv.lt(alpha).And(pvQ.lt(alpha)).And(cmap.eq(i.subtract(1)))
    # Then update cmap...
    cmap = cmap.where(tst, cmapj)
    # ...and fmap...
    fmap = fmap.where(tst, fmap.add(1))
    # ...and smap only if in first row.
    smap = ee.Algorithms.If(i.eq(1), smap.where(tst, cmapj), smap)
    # Create bmap band and add it to bmap image.
    idx = i.add(j).subtract(2)
    tmp = bmap.select(idx)
    bname = bmap.bandNames().get(idx)
    tmp = tmp.where(tst, 1)
    tmp = tmp.rename([bname])
    bmap = bmap.addBands(tmp, [bname], True)
    return ee.Dictionary({'i': i, 'j': j.add(1), 'alpha': alpha, 'pvQ': pvQ,
                          'cmap': cmap, 'smap': smap, 'fmap': fmap, 'bmap':bmap})

def filter_i(current, prev):
    """Arranges calculation of change maps; iterates over row-indices of pv_arr."""
    current = ee.List(current)
    pvs = current.slice(0, -1 )
    pvQ = ee.Image(current.get(-1))
    prev = ee.Dictionary(prev)
    i = ee.Number(prev.get('i'))
    alpha = ee.Image(prev.get('alpha'))
    median = prev.get('median')
    # Filter Ql p value if desired.
    pvQ = ee.Algorithms.If(median, pvQ.focal_median(1.5), pvQ)
    cmap = prev.get('cmap')
    smap = prev.get('smap')
    fmap = prev.get('fmap')
    bmap = prev.get('bmap')
    first = ee.Dictionary({'i': i, 'j': 1, 'alpha': alpha ,'pvQ': pvQ,
                           'cmap': cmap, 'smap': smap, 'fmap': fmap, 'bmap': bmap})
    result = ee.Dictionary(ee.List(pvs).iterate(filter_j, first))
    return ee.Dictionary({'i': i.add(1), 'alpha': alpha, 'median': median,
                          'cmap': result.get('cmap'), 'smap': result.get('smap'),
                          'fmap': result.get('fmap'), 'bmap': result.get('bmap')})
    
def dmap_iter(current, prev):
    """Reclassifies values in directional change maps."""
    prev = ee.Dictionary(prev)
    j = ee.Number(prev.get('j'))
    image = ee.Image(current)
    avimg = ee.Image(prev.get('avimg'))
    diff = image.subtract(avimg)
    # Get positive/negative definiteness.
    posd = ee.Image(diff.select(0).gt(0).And(det(diff).gt(0)))
    negd = ee.Image(diff.select(0).lt(0).And(det(diff).gt(0)))
    bmap = ee.Image(prev.get('bmap'))
    bmapj = bmap.select(j)
    dmap = ee.Image.constant(ee.List.sequence(1, 3))
    bmapj = bmapj.where(bmapj, dmap.select(2))
    bmapj = bmapj.where(bmapj.And(posd), dmap.select(0))
    bmapj = bmapj.where(bmapj.And(negd), dmap.select(1))
    bmap = bmap.addBands(bmapj, overwrite=True)
    # Update avimg with provisional means.
    i = ee.Image(prev.get('i')).add(1)
    avimg = avimg.add(image.subtract(avimg).divide(i))
    # Reset avimg to current image and set i=1 if change occurred.
    avimg = avimg.where(bmapj, image)
    i = i.where(bmapj, 1)
    return ee.Dictionary({'avimg': avimg, 'bmap': bmap, 'j': j.add(1), 'i': i})

def change_maps(im_list, median=False, alpha=0.01):
    """Calculates thematic change maps."""
    k = im_list.length()
    # Pre-calculate the P value array.
    pv_arr = ee.List(p_values(im_list))
    # Filter P values for change maps.
    cmap = ee.Image(im_list.get(0)).select(0).multiply(0)
    bmap = ee.Image.constant(ee.List.repeat(0,k.subtract(1))).add(cmap)
    alpha = ee.Image.constant(alpha)
    first = ee.Dictionary({'i': 1, 'alpha': alpha, 'median': median,
                           'cmap': cmap, 'smap': cmap, 'fmap': cmap, 'bmap': bmap})
    result = ee.Dictionary(pv_arr.iterate(filter_i, first))
    # Post-process bmap for change direction.
    bmap =  ee.Image(result.get('bmap'))
    avimg = ee.Image(im_list.get(0))
    j = ee.Number(0)
    i = ee.Image.constant(1)
    first = ee.Dictionary({'avimg': avimg, 'bmap': bmap, 'j': j, 'i': i})
    dmap = ee.Dictionary(im_list.slice(1).iterate(dmap_iter, first)).get('bmap')
    return ee.Dictionary(result.set('bmap', dmap))

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
    value='contract/',
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
    layout = widgets.Layout(width='200px'),
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
    value=False,
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
w_review = widgets.Button(description="ReviewAsset",disabled=False)
w_goto = widgets.Button(description='GoTo',disabled=False)
w_export_ass = widgets.Button(description='ExportToAssets',disabled=True)
w_export_drv = widgets.Button(description='ExportToDrive',disabled=True)
w_plot = widgets.Button(description='PlotAsset',disabled=False)
w_reset = widgets.Button(description='Reset',disabled=False)


w_masks = widgets.VBox([w_maskchange,w_maskwater,w_dw,w_quick])
w_dates = widgets.VBox([w_startdate,w_enddate])
w_assets = widgets.VBox([w_review,w_plot,w_reset])
w_bmap = widgets.VBox([w_interval,w_maxfreq,w_minfreq])
w_export = widgets.VBox([widgets.HBox([w_export_ass,w_exportassetsname]),widgets.HBox([w_export_drv,w_exportdrivename])])
w_signif = widgets.VBox([w_significance,w_median])

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
w_changemap.observe(on_changemap_widget_change,names='value')

w_stride_opac = widgets.VBox([w_stride,w_opacity])  

row1 = widgets.HBox([w_platform,w_orbitpass,w_relativeorbitnumber,w_dates])
row2 = widgets.HBox([w_collect,w_signif,w_stride_opac,w_export])
row3 = widgets.HBox([w_preview,w_changemap,w_bmap,w_masks,w_assets])
row4 = widgets.HBox([w_out,w_goto,w_location])

box = widgets.VBox([row1,row2,row3,row4])

dyn = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
                    .filterDate('2021-09-01', '2022-03-30') \
                    .select('label').mosaic()

    
def maskDynamicWorld(image): 
    return image.eq(ee.Image.constant(6)) # Built 

#watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)
watermask = dyn.gt(ee.Image.constant(0))

def GetTileLayerUrl(image):
    map_id = ee.Image(image).getMapId()
    return map_id["tile_fetcher"].url_format 


def handle_draw(self, action, geo_json):
    global poly
    coords =  geo_json['geometry']['coordinates']  
#     coords = [[[-96.06878489327627, 29.823701939611126],
#           [-96.06878489327627, 29.1423007492751],
#           [-95.00860422921377, 29.1423007492751],
#           [-95.00860422921377, 29.823701939611126]]]
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
    return ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/FNF') \
                  .filterDate('2017-01-01', '2017-12-31') 
                  
def get_vvvh(image):   
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values '''
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()

def convert_timestamp_list(tsl):
    # Make timestamps in YYYYMMDD format            
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
        
def on_reset_button_clicked(b):
    try:
        clear_layers()
        w_out.clear_output()
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_reset.on_click(on_reset_button_clicked)              

#@title Collect
def on_collect_button_clicked(b):
    ''' Collect a time series from the archive 
    '''
    global result, count, timestamplist, archive_crs
    with w_out:
        try:
            w_out.clear_output()
            print('running on GEE archive COPERNICUS/S1_GRD')
            
            collection = getS1collection()          
            if w_platform.value != 'Both':
                collection = collection.filter(ee.Filter.eq('platform_number', w_platform.value))
            if w_relativeorbitnumber.value > 0:
                collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(w_relativeorbitnumber.value)))  
            collection = collection.sort('system:time_start')  
            count = collection.size().getInfo()      
            if count<2: 
                raise ValueError('Less than 2 images found') 
#            footprint = ee.Geometry.Polygon(collection.first().get('system:footprint').getInfo()['coordinates'])
            print('Images found: %i, platform: %s, stride: %i'%(count,w_platform.value,w_stride.value))
            acquisition_times = ee.List(collection.aggregate_array('system:time_start')).getInfo()
            archive_crs = ee.Image(collection.first()).select(0).projection().crs().getInfo()
            timestamplist = []
            for timestamp in acquisition_times:
                tmp = time.gmtime(int(timestamp)/1000)
                timestamplist.append(time.strftime('%x', tmp))  
            #Make timestamps in YYYYMMDD format       
            timestamplist = convert_timestamp_list(timestamplist) 
            print(timestamplist)    
             
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
            S1 = collection.mean().select(0)                  
            #Run the algorithm ************************************************
            result = change_maps(imList, w_median.value, w_significance.value)
            #******************************************************************
            w_preview.disabled = False
            w_export_ass.disabled = False
            w_export_drv.disabled = False
            #Display preview 
            if len(rons)==1:
                print( 'please wait for raster overlay ...' )
                clear_layers()
                m.add_layer(TileLayer(url=GetTileLayerUrl(S1.visualize(min=-15, max=4)),name='S1'))   
            collectionS2 = getS2collection()      
            if collectionS2.size().getInfo() > 0:    
                s2_image =  ee.Image(collectionS2.first()).clip(poly) \
                                                            .select(['B4', 'B3', 'B2'])                                                              
                timestamp = s2_image.get('system:time_start').getInfo() 
                timestamp = time.gmtime(int(timestamp)/1000)
                timestamp = time.strftime('%x', timestamp).replace('/','')
                timestamps2 = '20'+timestamp[4:]+timestamp[0:4]
                print('Best Sentinel 2 from %s'%timestamps2)
                m.add_layer(TileLayer(url=GetTileLayerUrl(s2_image.visualize(min=0.0,max=5000)),name='S2 RGB'))
            else:
                print('No S2 image found')                                    
        except Exception as e:
            print('Error: %s'%e) 

w_collect.on_click(on_collect_button_clicked)                  

#@title Preview
def on_preview_button_clicked(b):
    ''' Preview change maps
    '''
    with w_out:  
        try:       
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            ggg = 'black,green'
            palette = jet
            w_out.clear_output()
            print('Series length: %i images, previewing (please wait for raster overlay) ...'%count)
            if w_changemap.value=='First':
                mp = ee.Image(result.get('smap')).byte()
                
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))
                
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp=ee.Image(result.get('cmap')).byte()
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = ee.Image(result.get('fmap')).byte()
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                mx = w_maxfreq.value
                print('Change frequency :\n blue = few, red = many')
            elif w_changemap.value == 'Bitemporal':
                sel = int(w_interval.value)
                sel = min(sel,count-1)
                sel = max(sel,1)
                print('Bitemporal: %s-->%s'%(timestamplist[sel-1],timestamplist[sel]))
                print('red = positive definite, cyan = negative definite, yellow = indefinite')  
                bmap = ee.Image(result.get('bmap')).byte()   
                mp = bmap.select(sel-1)
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                palette = rcy
                mx = 3     
            if len(m.layers)>9:
                m.remove_layer(m.layers[9])
            if not w_quick.value:
                mp = mp.reproject(crs=archive_crs,scale=float(w_exportscale.value))
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

w_preview.on_click(on_preview_button_clicked)      

def on_review_button_clicked(b):
    ''' Examine change maps exported to user's assets
    ''' 
    with w_out:  
        try: 
#          test for existence of asset                  
            _ = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value).getInfo()
#          ---------------------------            
            asset = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value)
            poly = ee.Geometry.Polygon(ee.Geometry(asset.get('system:footprint')).coordinates())
            center = poly.centroid().coordinates().getInfo()
            center.reverse()
#            m.center = center  
            bnames = asset.bandNames().getInfo()[3:]
            count = len(bnames)               
            jet = 'black,blue,cyan,yellow,red'
            rcy = 'black,red,cyan,yellow'
            ggg = 'black,green'
            smap = asset.select('smap').byte()
            cmap = asset.select('cmap').byte()
            fmap = asset.select('fmap').byte()
            bmap = asset.select(list(range(3,count+3)),bnames).byte()      
            palette = jet
            w_out.clear_output()
            print('Series length: %i images, reviewing (please wait for raster overlay) ...'%(count+1))
            if w_changemap.value=='First':
                mp = smap
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                mx = count
                print('Interval of first change:\n blue = early, red = late')
            elif w_changemap.value=='Last':
                mp = cmap
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                mx = count
                print('Interval of last change:\n blue = early, red = late')
            elif w_changemap.value=='Frequency':
                mp = fmap
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                mx = w_maxfreq.value
                print('Change frequency :\n blue = few, red = many')
            elif w_changemap.value == 'Bitemporal':
                sel = int(w_interval.value)
                sel = min(sel,count-1)
                sel = max(sel,1)
                print('Bitemporal: interval %i'%w_interval.value)
                print('red = positive definite, cyan = negative definite, yellow = indefinite')  
                mp = bmap.select(sel-1)
                if w_dw.value:
                    mp = mp.mask(mp.mask().And(maskDynamicWorld(dyn)))                
                palette = rcy
                mx = 3     
            if len(m.layers)>9:
                m.remove_layer(m.layers[9])   
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
    ''' Export to assets
    '''
    try:
        smap = ee.Image(result.get('smap')).byte()
        cmap = ee.Image(result.get('cmap')).byte()
        fmap = ee.Image(result.get('fmap')).byte() 
        bmap = ee.Image(result.get('bmap')).byte() 
#      in case of duplicates                  
        timestamplist1 = [timestamplist[i] + '_' + str(i+1) for i in range(len(timestamplist))]             
        cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist1[1:])  
        assexport = ee.batch.Export.image.toAsset(cmaps.byte().clip(poly),
                                    description='assetExportTask', 
                                    pyramidingPolicy={".default": 'sample'},
                                    assetId='projects/ee-mortcanty/assets/'+w_exportassetsname.value,scale=10,maxPixels=1e9)      
        assexport.start()
        with w_out: 
            w_out.clear_output() 
            print('Exporting change maps to %s\n task id: %s'%('projects/ee-mortcanty/assets/'+w_exportassetsname.value,str(assexport.id)))
    except Exception as e:
        with w_out:
            print('Error: %s'%e)                                          
    
w_export_ass.on_click(on_export_ass_button_clicked)  

def on_export_drv_button_clicked(b):
    ''' Export to Google Drive
    '''
    try:
        smap = ee.Image(result.get('smap')).byte()
        cmap = ee.Image(result.get('cmap')).byte()
        fmap = ee.Image(result.get('fmap')).byte() 
        bmap = ee.Image(result.get('bmap')).byte()            
        cmaps = ee.Image.cat(cmap,smap,fmap,bmap).rename(['cmap','smap','fmap']+timestamplist[1:])  
        fileNamePrefix=w_exportdrivename.value.replace('/','-')            
        gdexport = ee.batch.Export.image.toDrive(cmaps.byte().clip(poly),
                                    description='driveExportTask', 
                                    folder = 'gee',
                                    fileNamePrefix=fileNamePrefix,scale=10,maxPixels=1e9)   
        gdexport.start()
        with w_out:
            w_out.clear_output()
            print('Exporting change maps to Drive/gee/%s\n task id: %s'%(fileNamePrefix,str(gdexport.id))) 
    except Exception as e:
        with w_out:
            print('Error: %s'%e) 

w_export_drv.on_click(on_export_drv_button_clicked)            

def on_plot_button_clicked(b):          
#  plot change fractions        
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
            assetImage = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value)
            k = assetImage.bandNames().length().subtract(3).getInfo()            
            bmap1 = assetImage.select(ee.List.sequence(3,k+2)).clip(poly)            
            if w_maskwater.value:
                bmap1 = bmap1.updateMask(watermask) 
            if w_dw.value:
                bmap1 = bmap1.mask(bmap1.mask().And(maskDynamicWorld(dyn))) 
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
#            fn = 'projects/ee-mortcanty/assets/'+w_exportassetsname.value.replace('/','-')+'.png'
#            plt.savefig(fn,bbox_inches='tight') 
            w_out.clear_output()
            plt.show()
#            print('Saved to ~/%s'%fn)
        except Exception as e:
            print('Error: %s'%e)               
    
w_plot.on_click(on_plot_button_clicked)

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

    palette = ['419BDF','397D49','88B053','7A87C6','E49635','DFC35A','C4281B','A59B8F','B39FE1']
    m.add_layer( TileLayer(url=GetTileLayerUrl(dyn.visualize(min=0,max=8,
                                  palette=palette)),name= 'Dynamic World') )        
        
    display(m) 
    return box      