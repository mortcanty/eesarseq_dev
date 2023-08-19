# eesarseq.py
# Utilities to assemble SAR time series and
# run the non-sequential omnibus and 
# sequential change detection algorithms
#
# Mort Canty
# August, 2022

import time
import math
import numpy as np
from collections import Counter

import ee
ee.Initialize()

# *******************************************************
# The non-sequential (omnibus) change detection algorithm
# *******************************************************

def det(im):
    """Calculates the  determinant of 2x2 diagonal covariance matrix."""
    return ee.Image(im).expression('b(0)*b(1)')

def chi2cdf(chi2, df):
    """Calculates Chi square cumulative distribution function for
       df degrees of freedom using the built-in incomplete gamma
       function gammainc().
    """
    return ee.Image(chi2.divide(2)).gammainc(ee.Number(df).divide(2))

def logdet(im):
    """Calculates the log of the determinant of 2x2 diagonal covariance matrix."""
    return ee.Image(im).expression('b(0)*b(1)').log()

def log_det_sum_all(im_list):
    """Returns log of determinant of the sum of the  images in im_list."""
    im_list = ee.List(im_list)
    suml = ee.ImageCollection(im_list).reduce(ee.Reducer.sum())
    return ee.Image(det(suml)).log()

def omnibus(im_list):
    """
    Calculates the omnibus test statistic Q, bivariate case.
    k = series length
    f = 2*(k-1)
    m = ENL
    rho = 1 - (k**2 - 1)/3mfk
    omega = f*(1-1/rho)**2/4 
    Returns -2rho*logQ, omega
    """ 
    m = 4.4
    im_list = ee.List(im_list)
    k = im_list.length().getInfo()
    f = 2*(k-1)
    rho = 1 - (k**2-1)/(3*m*f*k)
    omega = ee.Number((f*(1-1/rho)**2)/4)       
#    k2logk = k.multiply(k.log()).multiply(2)
    
    k2logk = 2*k*np.log(k) 
    
    k2logk = ee.Image.constant(k2logk)
    sumlogdets = ee.ImageCollection(im_list.map(logdet)).reduce(ee.Reducer.sum())
    logdetsum = log_det_sum_all(im_list)   
    return ( k2logk.add(sumlogdets).subtract(logdetsum.multiply(k)).multiply(2*m).multiply(-rho), omega )

def change_map(im_list, median, alpha):
    '''Returns omnibus change map for im_list''' 
    k = im_list.length()
    f = k.subtract(1).multiply(2)
    f4 = f.add(4)   
    m2rhologQ, omega = omnibus(im_list) 
    c2 = chi2cdf(m2rhologQ, f)
    c4 = chi2cdf(m2rhologQ, f4)   
    p_value = ee.Image.constant(1).subtract(c2).subtract(c4.subtract(c2).multiply(omega))  
    if median:
        p_value = p_value.focal_median(1.5)
    return p_value.multiply(0).where(p_value.lt(alpha), 1)

# *****************************************
# The sequential change detection algorithm
# *****************************************

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

# ******************
# Refined Lee Filter
# ******************

def rl(img):
    '''
    Refined Lee Speckle Filter for S1 images only
    Created on 03.03.2020
    Transcribed from Guido Lemoine's 2017 JavaScript 
    implementation on the GEE
    '''
#   img must be in natural units, and single band
#   Set up 3x3 kernels 
    weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
    kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

    mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
    variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

#  Use a sample of the 3x3 windows inside a 7x7 window to determine gradients and directions
    sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], 
            [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

    sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False)

#  Calculate mean and variance for the sampled windows and store as 9 bands
    sample_mean = mean3.neighborhoodToBands(sample_kernel) 
    sample_var = variance3.neighborhoodToBands(sample_kernel)

#  Determine the 4 gradients for the sampled windows
    gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
    gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
    gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
    gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

#  And find the maximum gradient amongst gradient bands
    max_gradient = gradients.reduce(ee.Reducer.max())

#  Create a mask for band pixels that are the maximum gradient
    gradmask = gradients.eq(max_gradient)

#  duplicate gradmask bands: each gradient represents 2 directions
    gradmask = gradmask.addBands(gradmask)

#  Determine the 8 directions
    directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
    directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
    directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
    directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
#  The next 4 are the not() of the previous 4
    directions = directions.addBands(directions.select(0).Not().multiply(5))
    directions = directions.addBands(directions.select(1).Not().multiply(6))
    directions = directions.addBands(directions.select(2).Not().multiply(7))
    directions = directions.addBands(directions.select(3).Not().multiply(8))

#  Mask all values that are not 1-8
    directions = directions.updateMask(gradmask)

#  "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
    directions = directions.reduce(ee.Reducer.sum())

#  var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
#  Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

    sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

#  Calculate localNoiseVariance
    sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

#  Set up the 7*7 kernels for directional statistics
    rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

    diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

    rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
    diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

#  Create stacks for mean and variance using the original kernels. Mask with relevant direction.
    dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
    dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

#  and add the bands for rotated kernels
    for i in range(1,4):
        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
        dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
        dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))


#  "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
    dir_mean = dir_mean.reduce(ee.Reducer.sum())
    dir_var = dir_var.reduce(ee.Reducer.sum())

#  And finally generate the filtered value
    varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))

    b = varX.divide(dir_var)

    result = dir_mean.add(b.multiply(img.subtract(dir_mean)))
    
    return result.arrayFlatten([['sum']])

def refinedLee(img):
    return ee.Image.cat(rl(img.select(0)), rl(img.select(1)))

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
    and restore linear signal from db-values '''   
    return image.select('VV', 'VH').multiply(ee.Image.constant(math.log(10.0)/10.0)).exp()

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
    tsl= [x.replace('/', '') for x in tsl]
    tsl = ['T20' + x[4:] + x[0:4] for x in tsl]
    return tsl

def clipList(current, prev):
    ''' clip a list of images and multiply by ENL'''
    imlist = ee.List(ee.Dictionary(prev).get('imlist'))
    aoi = ee.Dictionary(prev).get('aoi')
    enl = ee.Number(ee.Dictionary(prev).get('enl'))
    stride = ee.Number(ee.Dictionary(prev).get('stride'))
    ctr = ee.Number(ee.Dictionary(prev).get('ctr'))
    imlist =  ee.Algorithms.If(ctr.mod(stride).eq(0),
        imlist.add(ee.Image(current).multiply(enl).clip(aoi)),imlist)
    return ee.Dictionary({'imlist':imlist,'aoi':aoi,'enl':enl,'ctr':ctr.add(1),'stride':stride})  

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

def assemble_and_run(aoi, median = True, significance = 0.01, stride = 1, startdate = '20180101',
                     enddate = '20190101', platform='A', orbitpass='DESCENDING', ron=0):
    '''
    *** Collect a time series from all intersecting orbit paths and invoke algorithm ***
    Input:
    aoi          <list>     polygon region of interest
    median       <boolean>  if true apply 3x3 median filter to p-values (default: true)
    significance <float>    likelihood ratio test significance (default: 0.01)
    stride       <int>      image series stride (default: 1)
    startdate    <string>   start of time series (default: 20180101)
    enddate      <string>   end of time series (default: 20190101)
    platform     <string>   sentinel-1 satellite (default: A)
    orbitpass    <string>   ascending or descending node (default: DESCENDING)
    ron          <int>      rel orbit number, if 0 use all intersecting orbits (default: 0)
    Output:
    cmaps        <ee.Image><byte>            3 bands: cmap, smap, fmap
    bmaps        <ee.Image><byte>            count-1 bands: bmap
    count        <int>                       length of time series
    rons         <list>                      relative orbit numbers used
    collection   <ee.ImageCollection><float> the filtered image collection
    atsf         <ee.Image><float>           ATSF 
    sequence     <ee.Image><float>           Lee-filtered time sequence, 2*count bands 
    '''
    def trim_list(current):
        current = ee.Image(current)
        bns = current.bandNames().slice(0, count-1)
        return current.select(ee.List.sequence(0, count-2), bns)
    
    def trim_sequence(current):
        current = ee.Image(current)
        bns = current.bandNames().slice(0, 2*count)
        return current.select(bns) 
    
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
        omap_list = ee.List([])
        sequence_list = ee.List([])
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
            print(list(zip(uniquetimestamps, orbit_lengths.getInfo())))
            # get list of all images for the present orbit number and convert from decibels
            cList = collection_ron.map(get_vvvh).toList(500)    
            # make list of combined (mosaicked) images along orbit path          
            first = ee.Dictionary({'plist': ee.List([]), 'clist': cList})
            pList = ee.List(ee.Dictionary(orbit_lengths.iterate(make_mosaics, first)).get('plist'))         
            # construct Lee-filtered time series for this orbit
            sequence = ee.ImageCollection(pList) \
                            .map(refinedLee) \
                            .toBands() \
                            .clip(aoi) \
                            .float()                        
            sequence_list = sequence_list.add(sequence)                                
            # clip time series to aoi and multiply by enl
            first = ee.Dictionary({'imlist': ee.List([]), 'enl': ee.Number(4.4),
                                   'aoi': aoi, 'ctr': ee.Number(0), 'stride': ee.Number(int(stride))})
            imList = ee.List(ee.Dictionary(pList.iterate(clipList, first)).get('imlist'))           
            # length of shortest time series so far
            count = minimum(imList.size().getInfo(), count)
            # run the algorithm for this relative orbit #####
            result = change_maps(imList, median, significance) # sequential
            omap = change_map(imList, median, significance).byte() # non-sequential
            #################################################
            smap = ee.Image(result.get('smap')).byte()        
            cmap = ee.Image(result.get('cmap')).byte()
            fmap = ee.Image(result.get('fmap')).byte()
            bmap = ee.Image(result.get('bmap')).byte()
            atsf = ee.Image(result.get('atsf')).float()
            # combine to lists and clip to reduced geometry
            cmap_list = cmap_list.add(ee.Image.cat(cmap, smap, fmap)
                                        .rename(['cmap', 'smap', 'fmap']).clip(geo))
            bmap_list = bmap_list.add(bmap.clip(geo))
            atsf_list = atsf_list.add(atsf.clip(geo))
            omap_list = omap_list.add(omap.clip(geo))
        # mosaic cmap, smap, fmap images
        cmaps = ee.ImageCollection.fromImages(cmap_list).mosaic()
        # truncate bitemporal maps to length of shortest series
        bmap_list = bmap_list.map(trim_list)
        # mosaick the bitemporal maps labeling bands with most recent timestamp list 
        bmaps = ee.ImageCollection(bmap_list).mosaic().rename(uniquetimestamps[stride:count*stride:stride])
        # mosaic ATSF images
        atsf = ee.ImageCollection.fromImages(atsf_list).mosaic()
        # mosaic omnibus change maps
        omaps = ee.ImageCollection.fromImages(omap_list).mosaic()
        # mosaic time sequences  
        # truncate sequences to length of shortest sequence
        sequence_list = sequence_list.map(trim_sequence)   
        time_sequence = ee.ImageCollection.fromImages(sequence_list).mosaic()     
        # return results
        return (cmaps, bmaps, count, rons, collection, atsf, omaps, time_sequence)
    except Exception as e:
        print('Error: %s'%e)

if __name__ == '__main__':
    pass
