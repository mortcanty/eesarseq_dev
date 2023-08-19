# forest.py
# widget interface for SAR forest mapping with omnibus change detection

import ee
ee.Initialize
import time, math

from eesar.eesarseq import assemble_and_run

import ipywidgets as widgets
from IPython.display import display
from ipyleaflet import (Map,DrawControl,TileLayer,
                        MeasureControl,
                        FullScreenControl,
                        Polygon,
                        basemaps,basemap_to_tiles,
                        LayersControl)
from geopy.geocoders import Nominatim

# ********************
# The widget interface
# ********************

aoi = None
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
    options=['Omnibus','EWM','PALSAR','DW','Hansen'],
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
w_compareassetname = widgets.Text(
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
w_F1 = widgets.Checkbox(
    layout = widgets.Layout(width='180px'),
    value=False,
    description='F1 Scores',
    disabled=False
)
w_stride = widgets.BoundedIntText(
    value=1,
    min=1,
    description='Stride:',
    layout = widgets.Layout(width='150px'),
    disabled=False
)
w_shape = widgets.BoundedIntText(
    value=45,
    min=1,
    description='Index',
    layout = widgets.Layout(width='150px'),
    disabled=False
)
w_settlement = widgets.Checkbox(
    value=True,
    description='SettlementMask',
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
w_calcloss = widgets.Button(description='CalculateLoss',disabled=True)
w_export_drv = widgets.Button(description='ExportToDrive',disabled=True)
w_reset = widgets.Button(description='Reset',disabled=False)

w_masks = widgets.VBox([w_maskwater,w_settlement])
w_dates = widgets.VBox([w_startdate,w_enddate])
w_export = widgets.VBox([widgets.HBox([w_useshape,w_shape]),
                         widgets.HBox([w_export_ass,w_exportassetsname]),
                         widgets.HBox([w_calcloss,w_compareassetname])])
w_signif = widgets.VBox([w_significance,w_median])

def on_widget_change(b):
    w_preview.disabled = True
    w_export_ass.disabled = True
    w_calcloss.disabled = True


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
w_median.observe(on_widget_change,names='value')
w_significance.observe(on_widget_change,names='value')

row1 = widgets.HBox([w_platform,w_orbitpass,w_relativeorbitnumber,w_dates])
row2 = widgets.HBox([w_collect,w_signif,widgets.VBox([w_stride,w_F1]),w_export])
row3 = widgets.HBox([w_preview,w_changemap,w_masks,w_review,w_reset])
row4 = widgets.HBox([w_out,w_goto,w_location])

box = widgets.VBox([row1,row2,row3,row4])

#@title Collect

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
        w_calcloss.disabled = True
        w_export_drv.disabled = True
    elif action == 'deleted':
        aoi = None
        w_preview.disabled = True
        w_export_ass.disabled = True
        w_calcloss.disabled = True
        w_export_drv.disabled = True


def getPALSAR():
    PALSAR = ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/FNF') \
                  .filterDate('2017-01-01', '2017-12-31')
    mp = PALSAR.first().clip(aoi)
    prj = mp.projection()
    scale = prj.nominalScale()
    mp0 = mp.multiply(0)
    mp = mp0.where(mp.eq(1),1)#.selfMask()
    return mp.reproject(prj.atScale(scale))

def getDW():
    DW = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
                  .filterDate('2020-01-01', '2023-01-01') \
                  .filterBounds(aoi) \
                  .select('label') \
                  .reduce(ee.Reducer.mode())
    return DW.eq(1).clip(aoi)

def getEWM():
    mp = ee.ImageCollection("ESA/WorldCover/v100").first().clip(aoi)
    prj = mp.projection()
    scale = prj.nominalScale()
    mp0 = mp.multiply(0)
    mp = mp0.where(mp.eq(10),1)#.selfMask()
    return mp.reproject(prj.atScale(scale))

def getHansen(canopy=10,contiguous=9,loss_year=1):
    gfc = ee.Image('UMD/hansen/global_forest_change_2021_v1_9')
    treecover = gfc.select('treecover2000').clip(aoi)
    lossyear = gfc.select('lossyear').clip(aoi)
    loss = lossyear.multiply(0).where(lossyear.gte(loss_year),1).selfMask()
    forest2000 = treecover.gte(ee.Number(canopy)).selfMask()
    prj = forest2000.projection()
    scale = prj.nominalScale()
    mp = forest2000.connectedPixelCount(100).gte(ee.Number(contiguous)).selfMask()
    mp = mp.where(lossyear.gte(loss_year),0).selfMask()
    return (mp.reproject(prj.atScale(scale)) , loss.reproject(prj.atScale(scale)))

def clear_layers():
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

def on_collect_button_clicked(b):
    ''' Collect an S1 time series from the archive
    '''
    global omaps, count, timestamplist, aoi, prj
    with w_out:
        try:
            w_out.clear_output()
            clear_layers()
            if w_useshape.value:
                if w_shape.value > 10:
#              NRW Landkreise
                    print('NRW')
                    aois = ee.FeatureCollection('projects/ee-mortcanty/assets/dvg1krs_nw').geometry()
                    aoi = ee.Geometry(aois.geometries().get(w_shape.value))
                elif w_shape.value == 1:
#              Colombia
                    print('Colombia')
                    aois = ee.FeatureCollection('projects/ee-mortcanty/assets/colombia_admin0').geometry()
                    aoi = ee.Geometry(aois.geometries().get(-1))
                locations =  aoi.coordinates().getInfo()[0]
                locations = [tuple(list(reversed(i))) for i in locations]
                layer = Polygon(locations=locations,
                            color="white",
                            fill_color = 'black',
                            name = 'Shapefile')
                center = aoi.centroid().coordinates().getInfo()
                center.reverse()
                m.center = center
                m.add_layer(layer)
            #Assemble and run
            _, _, count, _, collection, _, omaps, _ = assemble_and_run(
                            aoi, median = w_median.value, significance = w_significance.value, stride = w_stride.value, startdate = w_startdate.value,
                            enddate = w_enddate.value, platform=w_platform.value, orbitpass=w_orbitpass.value, ron=w_relativeorbitnumber.value)
            #Get a preview as collection mean
            prj = ee.Image(collection.first()).select(0).projection()
            S1 = collection.mean().select(0).visualize(min=-15,max=4)
            w_preview.disabled = False
            w_export_ass.disabled = False
            w_calcloss.disabled = False
            w_export_drv.disabled = False
            #Display S1
            print( 'Length of time series: %i\nplease wait for raster overlay ...'%count )
            m.add_layer(TileLayer(url=GetTileLayerUrl(S1),name='S1'))
        except Exception as e:
            print('Error: %s'%e)

w_collect.on_click(on_collect_button_clicked)

watermask = ee.Image('UMD/hansen/global_forest_change_2015').select('datamask').eq(1)
settlement = ee.Image("DLR/WSF/WSF2015/v1")

mpo = None

def on_preview_button_clicked(b):
    ''' Preview change maps
    '''
    global prj
    with w_out:
        try:
            ggg = 'black,green'
            rrr = 'black,red'
            w_out.clear_output()
            print('Series length: %i images, previewing (please wait for raster overlay) ...'%count)
            if w_changemap.value=='Omnibus':
#              no-change map
                palette = ggg
                mp = ee.Image(omaps).int().Not().selfMask()
#              minimum contiguous area requirement (0.5 hectare)
                contArea = mp.connectedPixelCount().selfMask()
                mp = contArea.gte(ee.Number(50)).selfMask()
#              mask settled areas
                if w_settlement.value:
                    mp = mp.where(settlement.eq(255),0).selfMask()
#              forest cover in hectares
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
#              scale the map so not affected by zoom
                scale = 10
                mp = mp.reproject(prj.atScale(scale))
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = scale,
                                    maxPixels = 1e13)
#              dictionary of pixels counts (only one band in this case)
                pixelCount = mp.reduceRegion( ee.Reducer.count(), geometry = aoi ,scale = scale, maxPixels = 1e13 )
#              pixel size
                onePixel = forestCover.getNumber('constant').divide(pixelCount.getNumber('constant'))
                minAreaUsed = onePixel.multiply(50).getInfo()
                print('Omnibus change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('constant').getInfo()))
                print('Minimum forest area used (ha) ', minAreaUsed)
                if w_F1.value:
#              F1 score with Hansen as ground truth
                    mp = mp.unmask()
                    mph, _ = getHansen()
                    mph = mph.clip(aoi).unmask()
                    TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1h = 2.0*P*R/(P+R)
#              F1 score with PALSAR as ground truth
                    mpp = getPALSAR().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpp.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpp.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    FN = mp.multiply(0).where(mpp.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1p = 2.0*P*R/(P+R)
#              F1 score with EWM as ground truth
                    mpe = getEWM().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpe.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpe.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    FN = mp.multiply(0).where(mpe.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1e = 2.0*P*R/(P+R)
#              F1 score with DW as ground truth
                    mpd = getDW().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpd.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpd.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    FN = mp.multiply(0).where(mpd.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1d = 2.0*P*R/(P+R)
                    print('F1 score relative to Hansen ', F1h)
                    print('F1 score relative to PALSAR ', F1p)
                    print('F1 score relative to EWM ', F1e)
                    print('F1 score relative to DW ', F1d)

            elif w_changemap.value=='EWM':
                mp = getEWM().clip(aoi).unmask()
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = 10,
                                    maxPixels = 1e13)
                print('EWM change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('Map').getInfo()))
                if w_F1.value:
#              F1 with hansen as ground truth
                    mph, _ = getHansen()
                    mph = mph.clip(aoi).unmask()
                    TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mph.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    FN = mp.multiply(0).where(mph.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1h = 2.0*P*R/(P+R)

#              F1 with PALSAR as ground truth
                    mpp = getPALSAR().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpp.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpp.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    FN = mp.multiply(0).where(mpp.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1p = 2.0*P*R/(P+R)

#              F1 with DW as ground truth
                    mpd = getDW().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpd.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpd.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    FN = mp.multiply(0).where(mpd.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('Map')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1d = 2.0*P*R/(P+R)

                    print('F1 score relative to Hansen ', F1h)
                    print('F1 score relative to PALSAR ', F1p)
                    print('F1 score relative to DW', F1d)

            elif w_changemap.value=='PALSAR':
                mp = getPALSAR().clip(aoi).unmask()
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = 30,
                                    maxPixels = 1e13)
                print('PALSAR change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('fnf').getInfo()))
                if w_F1.value:
#              F1 with Hansen as ground truth
                    mph, _ = getHansen()
                    mph = mph.clip(aoi).unmask()
                    TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mph.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    FN = mp.multiply(0).where(mph.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1h = 2.0*P*R/(P+R)

#              F1 with EWM as ground truth
                    mpe = getEWM().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpe.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpe.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    FN = mp.multiply(0).where(mpe.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1e = 2.0*P*R/(P+R)

#              F1 with DW as ground truth
                    mpd = getDW().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpd.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpd.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    FN = mp.multiply(0).where(mpd.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('fnf')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1d = 2.0*P*R/(P+R)

                    print('F1 score relative to Hansen ', F1h)
                    print('F1 score relative to EWM ', F1e)
                    print('F1 score relative to DW ', F1d)
            elif w_changemap.value=='DW':
                mp = getDW().clip(aoi).unmask()
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = 10,
                                    maxPixels = 1e13)
                print('DW change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('label_mode').getInfo()))
                if w_F1.value:
                    mph, _ = getHansen()
                    mph = mph.clip(aoi).unmask()
                    TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mph.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    FN = mp.multiply(0).where(mph.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1 = 2.0*P*R/(P+R)

                    mpp = getPALSAR().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpp.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpp.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    FN = mp.multiply(0).where(mpp.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1p = 2.0*P*R/(P+R)

                    mpe = getEWM().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpe.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpe.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    FN = mp.multiply(0).where(mpe.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('label_mode')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1e = 2.0*P*R/(P+R)

                    print('F1 score relative to Hansen ', F1)
                    print('F1 score relative to PALSAR ', F1p)
                    print('F1 score relative to EWM ', F1e)
            elif w_changemap.value=='Hansen':
                mp, loss = getHansen(loss_year=1)
                mp = mp.clip(aoi).unmask()
                loss = loss.clip(aoi).unmask()
                palette = ggg
                pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
                forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = 30,
                                    maxPixels = 1e13)
                pixelCount = mp.reduceRegion( ee.Reducer.count(), geometry = aoi ,scale = 30, maxPixels = 1e13 )
                onePixel = forestCover.getNumber('treecover2000').divide(pixelCount.getNumber('treecover2000'))
                minAreaUsed = onePixel.multiply(9).getInfo()
                print('Hansen change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('treecover2000').getInfo()))
                loss = loss.updateMask(loss.gt(0))
                m.add_layer(TileLayer(url=GetTileLayerUrl(loss.visualize(min=0, max=1,  palette=rrr)), name=w_changemap.value+' loss'))
                if w_F1.value:
                    mpd = getDW().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpd.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpd.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    FN = mp.multiply(0).where(mpd.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1d = 2.0*P*R/(P+R)

                    mpp = getPALSAR().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpp.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpp.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    FN = mp.multiply(0).where(mpp.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1p = 2.0*P*R/(P+R)

                    mpe = getEWM().clip(aoi).unmask()
                    TP = mp.multiply(0).where(mpe.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    TP = ee.Number(TP)
                    FP = mp.multiply(0).where(mpe.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    FN = mp.multiply(0).where(mpe.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('treecover2000')
                    P = TP.divide(TP.add(FP)).getInfo()
                    R = TP.divide(TP.add(FN)).getInfo()
                    F1e = 2.0*P*R/(P+R)

                    print('Minimum forest area used (ha) ', minAreaUsed)
                    print('F1 score relative to PALSAR ', F1p)
                    print('F1 score relative to EWM ', F1e)
                    print('F1 score relative to DW ', F1d)

            if len(m.layers)>8:
                m.remove_layer(m.layers[8])
            if w_maskwater.value==True:
                mp = mp.updateMask(watermask)
            mp = mp.updateMask(mp.gt(0))
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.visualize(min=0, max=1,  palette=palette)), name=w_changemap.value))
        except Exception as e:
            print('Error: %s'%e)

w_preview.on_click(on_preview_button_clicked)

def on_review_button_clicked(b):
    ''' Examine change maps exported to user's assets
    '''
    global aoi
    with w_out:
        try:
#           test for existence of asset
            _ = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value).getInfo()
#           ---------------------------
            w_out.clear_output()
            mp = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value).unmask()
            aoi = ee.Geometry.Polygon(ee.Geometry(mp.get('system:footprint')).coordinates())
            center = aoi.centroid().coordinates().getInfo()
            center.reverse()
            m.center = center
#           forest cover in hectares
            pixelArea = mp.multiply(ee.Image.pixelArea()).divide(10000)
            scale = 10
            forestCover = pixelArea.reduceRegion(
                                    reducer = ee.Reducer.sum(),
                                    geometry = aoi,
                                    scale = scale,
                                    maxPixels = 1e13)
            print('Omnibus change map\nForest Cover (ha): %i'%math.trunc(forestCover.get('constant').getInfo()))
            m.add_layer(TileLayer(url=GetTileLayerUrl(mp.updateMask(mp.gt(0)).visualize(min=0, max=1, palette='black,green')),name='omnibus'))
            if w_F1.value:
                mph, _ = getHansen()
                mph = mph.unmask()
                TP = mp.multiply(0).where(mph.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mph.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')

                print(FP.getInfo())

                FN = mp.multiply(0).where(mph.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')

                print(FN.getInfo())

                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1h = 2.0*P*R/(P+R)

                FP = mp.multiply(0).where(mph.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')

                print(FP.getInfo())

                FN = mp.multiply(0).where(mph.eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')

                print(FN.getInfo())

                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1ho = 2.0*P*R/(P+R)

                mpe = getEWM()
                TP = mp.multiply(0).where(mpe.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('constant')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mpe.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpe.eq(1).And(mp.eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1e = 2.0*P*R/(P+R)

                FP = mp.multiply(0).where(mpe.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpe.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,scale=10,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1eo = 2.0*P*R/(P+R)

                mpd = getDW()
                TP = mp.multiply(0).where(mpd.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mpd.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpd.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1d = 2.0*P*R/(P+R)

                FP = mp.multiply(0).where(mpd.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpd.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1do = 2.0*P*R/(P+R)

                mpp = getPALSAR()
                TP = mp.multiply(0).where(mpp.eq(1).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                TP = ee.Number(TP)
                FP = mp.multiply(0).where(mpp.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpp.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1p = 2.0*P*R/(P+R)

                FP = mp.multiply(0).where(mpp.eq(1).And(mp.unmask().eq(0)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                FN = mp.multiply(0).where(mpp.unmask().eq(0).And(mp.eq(1)),1).reduceRegion(ee.Reducer.sum(),aoi,maxPixels = 1e13).get('constant')
                P = TP.divide(TP.add(FP)).getInfo()
                R = TP.divide(TP.add(FN)).getInfo()
                F1po = 2.0*P*R/(P+R)

                print('F1 score relative to Hansen ', F1h)
                print('F1 score relative to PALSAR ', F1p)
                print('F1 score relative to EWM ', F1e)
                print('F1 score relative to DW ', F1d)

                print('F1 score Hansen relative to Omnibus', F1ho)
                print('F1 score PALSAR relative to Omnibus', F1po)
                print('F1 score EWM relative to Omnibus', F1eo)
                print('F1 score DW relative to Omnibus', F1do)
        except Exception as e:
            print('Error: %s'%e)

w_review.on_click(on_review_button_clicked)

def on_export_ass_button_clicked(b):
    ''' Export to assets
    '''
    try:
        mp = ee.Image(omaps).int().Not().selfMask()
#      minimum contiguous area requirement (0.5 hectare)
        contArea = mp.connectedPixelCount().selfMask()
        mp = contArea.gte(ee.Number(50)).selfMask()
#      mask settled areas
        if w_settlement.value:
            mp = mp.where(settlement.eq(255),0).selfMask()
        assexport = ee.batch.Export.image.toAsset(mp.byte().clip(aoi),
                                description='assetExportTask',
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

def on_calcloss_button_clicked(b):
    ''' Calculate loss
    '''
    try:
        compare_map = ee.Image('projects/ee-mortcanty/assets/'+w_compareassetname.value)
        current_map = ee.Image('projects/ee-mortcanty/assets/'+w_exportassetsname.value)
        w_out.clear_output()
        loss_map = compare_map.multiply(0).where(compare_map.eq(1).And(current_map.unmask().eq(0)),1).selfMask()
        with w_out:
            print('Omnibus loss map')
        m.add_layer(TileLayer(url=GetTileLayerUrl(compare_map.visualize(min=0, max=1, palette='black,green')),name='compare'))
        m.add_layer(TileLayer(url=GetTileLayerUrl(current_map.visualize(min=0, max=1, palette='black,green')),name='current'))
        m.add_layer(TileLayer(url=GetTileLayerUrl(loss_map.visualize(min=0, max=1, palette='black,red')),name='loss'))
    except Exception as e:
        with w_out:
            print('Error: %s'%e)

w_calcloss.on_click(on_calcloss_button_clicked)



#@title Run the interface
def run():
    global m, center
    center = [51.0,6.4]
    osm = basemap_to_tiles(basemaps.OpenStreetMap.Mapnik)
    ews = basemap_to_tiles(basemaps.Esri.WorldStreetMap)
    ewi = basemap_to_tiles(basemaps.Esri.WorldImagery)

    mc = MeasureControl(position='topright',primary_length_unit = 'kilometers')

    dc = DrawControl(aoiline={},circlemarker={})
    dc.rectangle = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}
    dc.aoigon = {"shapeOptions": {"fillColor": "#0000ff","color": "#0000ff","fillOpacity": 0.05}}

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
