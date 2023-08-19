# eesarseq.py
# Utilities to assemble SAR time series and
# run the sequential change detection algorithm
# (c) Mort Canty, 2022

import time
import math
from collections import Counter

import ee
import concurrent
ee.Initialize()

# *****************************************
# The sequential change detection algorithm
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
    pv = ee.Image.constant(1).subtract(chi2cdf(rhojm2logRj,2) \
                             .add(chi2cdf(rhojm2logRj,6).multiply(omega2j)) \
                             .subtract(chi2cdf(rhojm2logRj,2).multiply(omega2j)))
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
        pvQl = ee.Image.constant(1).subtract(chi2cdf(rhom2logQl,f) \
                                             .add(chi2cdf(rhom2logQl,f.add(4)).multiply(omega2)) \
                                             .subtract(chi2cdf(rhom2logQl,f).multiply(omega2)))

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
    pvQ = ee.Algorithms.If(median, pvQ.focal_median(2.5), pvQ)
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
    atsf = ee.Image(prev.get('atsf'))
    diff = image.subtract(atsf)
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
    # Update atsf with provisional means.
    i = ee.Image(prev.get('i')).add(1)
    atsf = atsf.add(image.subtract(atsf).divide(i))
    # Reset atsf to current image and set i=1 if change occurred.
    atsf = atsf.where(bmapj, image)
    i = i.where(bmapj, 1)
    return ee.Dictionary({'atsf': atsf, 'bmap': bmap, 'j': j.add(1), 'i': i})

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
    atsf = ee.Image(im_list.get(0))
    j = ee.Number(0)
    i = ee.Image.constant(1)
    first = ee.Dictionary({'atsf': atsf, 'bmap': bmap, 'j': j, 'i': i})
    d = ee.Dictionary(im_list.slice(1).iterate(dmap_iter, first))
    dmap = d.get('bmap')
    atsf = ee.Image(d.get('atsf')).divide(4.4)
    return ee.Dictionary(result.set('bmap', dmap)).combine({'atsf': atsf})

# ****************************
# Assemble time series and run
# ****************************
def getS1collection(aoi, orbitpass, startdate, enddate):
    return ee.ImageCollection('COPERNICUS/S1_GRD') \
                      .filterBounds(aoi) \
                      .filterDate(ee.Date(startdate), ee.Date(enddate)) \
                      .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH'])) \
                      .filter(ee.Filter.eq('resolution_meters', 10)) \
                      .filter(ee.Filter.eq('instrumentMode', 'IW')) \
                      .filter(ee.Filter.eq('orbitProperties_pass', orbitpass))

def get_vvvh(image):
    ''' get 'VV' and 'VH' bands from sentinel-1 imageCollection
    and restore linear signal from db-values, keeping system:footprint '''   
    fp = image.get('system:footprint')
    return image.select('VV','VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp().set({'system:footprint': fp})

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
    tsl = ['T20' + x[4:] + x[0:4] for x in tsl]
    return tsl

def clipList(current, prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    aoi = ee.Dictionary(prev).get('aoi')
    enl = ee.Number(ee.Dictionary(prev).get('enl'))
    imlist = imlist.add(ee.Image(current).multiply(enl).clip(aoi))
    return ee.Dictionary({'imlist': imlist,'aoi': aoi,'enl': enl})

def make_mosaics(current, prev):
    ''' return equitemporal mosaicked images in plist '''
    mLen = ee.Number(current)
    prev = ee.Dictionary(prev)
    pList = ee.List(prev.get('plist'))
    cList = ee.List(prev.get('clist'))   
    # mosaic of all equitemporal images along orbit path    
    images_on_path = cList.slice(0, mLen)
    mosaic = ee.ImageCollection(images_on_path).mosaic()    
    pList = pList.add( mosaic )
    return ee.Dictionary({'plist':pList, 'clist':cList.slice(mLen)})

def assemble_and_run(aoi, median = True, significance = 0.01, startdate = '20180101',
                     enddate = '20190101', platform='A', orbitpass='DESCENDING', ron=0):
    '''
    *** Collect a time series from all intersecting orbit paths and invoke algorithm ***
    Input:
    aoi          <list>     polygon region of interest
    median       <boolean>  if true apply 3x3 median filter to p-values (default: true)
    significance <float>    likelihood ratio test significance (default: 0.01)
    startdate    <string>   start of time series (default: 20180101)
    enddate      <string>   end of time series (default: 20190101)
    platform     <string>   sentinel-1 satellite (default: A)
    orbitpass    <string>   ascending of descending node (default: DESCENDING)
    ron          <int>      rel orbit number, if 0 use all intersecting orbits (default: 0)
    Output:
    cmaps        <ee.Image><byte>           3 bands: cmap, smap, fmap
    bmaps        <ee.Image><byte>           count-1 bands: bmap
    count        <int>                      length of time series
    rons         <list>                      relative orbit numbers used
    collection   <ee.ImageCollection><float> the filtered image collection
    atsf         <ee.Image><float>           ATSF  
    '''
    def trim_list(current):
        current = ee.Image(current)
        bns = current.bandNames().slice(0, count-1)
        return current.select(ee.List.sequence(0, count-2), bns)
    
    try:
        collection = getS1collection(aoi, orbitpass, startdate, enddate)   
        if platform != 'Both':
            collection = collection.filter(ee.Filter.eq('platform_number', platform))
        if ron > 0:
            collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', int(ron)))
        count = collection.size().getInfo()
        if count==0:
            raise ValueError('No images found')
        collection = collection.sort('system:time_start')       
        rons = map( int, ee.List(collection.aggregate_array('relativeOrbitNumber_start')) \
                    .getInfo() )
        rons = list(set(rons))
        rons.sort()
        cmap_list = ee.List([])
        bmap_list = ee.List([])
        atsf_list = ee.List([])
        count = 500
        for ron in rons:
            # filter for relative orbit number
            collection_ron = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', ron))    
            # to eliminate edge effects, get the bounding geometry 
            # of the images in the ron and shrink by 5km
            geo = collection_ron.geometry().dissolve().buffer(-5000)        
            # determine timestamps for the orbit
            timestamplist = get_timestamp_list(collection_ron)
            ctr = Counter(timestamplist)
            uniquetimestamps = list(set(timestamplist))
            uniquetimestamps.sort()
            print('Orbit number %i'%ron)
            # get orbit path lengths for each observation time (should all be same)
            # the path length determines the number of swaths which are to be mosaicked
            orbit_lengths = ee.List( [ctr[timestamp] for timestamp in uniquetimestamps] )
            print('Number of images in orbit path:')
            print(list(zip(uniquetimestamps,orbit_lengths.getInfo())))
            # get list of all images for the present orbit number
            cList = collection_ron.map(get_vvvh).toList(500)    
            # make list of combined (mosaicked) images along orbit path          
            first = ee.Dictionary({'plist': ee.List([]), 'clist': cList})
            pList = ee.List(ee.Dictionary(orbit_lengths.iterate(make_mosaics, first)).get('plist'))                      
            # clip time series to aoi and multiply by enl
            first = ee.Dictionary({'imlist':ee.List([]),'enl':ee.Number(4.4),'aoi':aoi})
            imList = ee.List(ee.Dictionary(pList.iterate(clipList, first)).get('imlist'))           
            # length of shortest time series so far
            count = minimum(imList.size().getInfo(), count)
            # run the algorithm for this relative orbit
            result = change_maps(imList, median, significance)            
            smap = ee.Image(result.get('smap')).byte()        
            cmap = ee.Image(result.get('cmap')).byte()
            fmap = ee.Image(result.get('fmap')).byte()
            bmap = ee.Image(result.get('bmap')).byte()
            atsf = ee.Image(result.get('atsf')).float()
            # combine to lists and clip to reduced geometry
            cmap_list = cmap_list.add( ee.Image.cat(cmap,smap,fmap).rename(['cmap','smap','fmap']).clip(geo) )
            bmap_list = bmap_list.add( bmap.clip(geo) )
            atsf_list = atsf_list.add( atsf.clip(geo) )
        # mosaic cmap, smap, fmap images
        cmaps = ee.ImageCollection.fromImages(cmap_list).mosaic()
        # truncate bitemporal maps to length of shortest series
        bmap_list = bmap_list.map(trim_list)
        # mosaick the bitemporal maps labeling bands with most recent timestamp list 
        bmaps = ee.ImageCollection(bmap_list).mosaic().rename(uniquetimestamps[1:count])
        # mosaic ATSF images
        atsf = ee.ImageCollection.fromImages(atsf_list).mosaic()
        return (cmaps, bmaps, count, rons, collection, atsf)
    except Exception as e:
        print('Error: %s'%e)

if __name__ == '__main__':
    pass
