"""Proof of concept code to run Mort's algo on EE

Need to manually authenticate beforehand using `earthengine authenticate`
"""

import sys
import getopt
import time
import math

from absl import app
from absl import flags
import IPython

from google3.geo.gestalt.client.python import ee

FLAGS = flags.FLAGS

flags.DEFINE_boolean('debug', False, 'Interactive shell.')

# *****************************************
# The sequental change detection algorithm
# *****************************************


def chi2cdf(chi2, df):
  """Calculates Chi square cumulative distribution function.

  Calculates this for df degrees of freedom using the built-in incomplete gamma
  function gammainc().

  Args:
    chi2:
    df:

  Returns:
    cumulative distribution function
  """
  return ee.Image(chi2.divide(2)).gammainc(ee.Number(df).divide(2))


def det(im):
  """Calculates determinant of 2x2 diagonal covariance matrix."""
  return im.expression('b(0)*b(1)')


def log_det_sum(im_list, j):
  """Returns log of determinant of the sum of the first j images in im_list."""
  im_list = ee.List(im_list)
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
  m2logRj = (
      log_det_sum(im_list, j.subtract(1)).multiply(j.subtract(1)).add(
          log_det(im_list, j)).add(ee.Number(2).multiply(j).multiply(
              j.log())).subtract(
                  ee.Number(2).multiply(j.subtract(1)).multiply(
                      j.subtract(1).log())).subtract(
                          log_det_sum(im_list,
                                      j).multiply(j)).multiply(-2).multiply(m))

  # correction to simple Wilks approximation
  # pv = ee.Image.constant(1).subtract(chi2cdf(m2logRj, 2))
  one = ee.Number(1)
  rhoj = one.subtract(
      one.add(one.divide(j.multiply(j.subtract(one)))).divide(6).divide(m))
  omega2j = one.subtract(one.divide(rhoj)).pow(2.0).divide(-2)
  rhojm2logRj = m2logRj.multiply(rhoj)
  pv = ee.Image.constant(1).subtract(
      chi2cdf(rhojm2logRj,
              2).add(chi2cdf(rhojm2logRj, 6).multiply(omega2j)).subtract(
                  chi2cdf(rhojm2logRj, 2).multiply(omega2j)))

  return (pv, m2logRj)


def p_values(im_list, m=4.4):
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
    # pvQl = ee.Image.constant(1).subtract(
    #     chi2cdf(m2logQl, ell.subtract(1).multiply(2)))
    one = ee.Number(1)
    f = ell.subtract(1).multiply(2)
    rho = one.subtract(
        ell.divide(m).subtract(one.divide(ell.multiply(m))).divide(f).divide(3))
    omega2 = f.multiply(one.subtract(one.divide(rho)).pow(2)).divide(-4)
    rhom2logQl = m2logQl.multiply(rho)
    pvQl = ee.Image.constant(1).subtract(
        chi2cdf(rhom2logQl,
                f).add(chi2cdf(rhom2logQl, f.add(4)).multiply(omega2)).subtract(
                    chi2cdf(rhom2logQl, f).multiply(omega2)))

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
  return ee.Dictionary({
      'i': i,
      'j': j.add(1),
      'alpha': alpha,
      'pvQ': pvQ,
      'cmap': cmap,
      'smap': smap,
      'fmap': fmap,
      'bmap': bmap
  })


def filter_i(current, prev):
  """Arranges calculation of change maps; iterates over row-indices of pv_arr."""
  current = ee.List(current)
  pvs = current.slice(0, -1)
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
  first = ee.Dictionary({
      'i': i,
      'j': 1,
      'alpha': alpha,
      'pvQ': pvQ,
      'cmap': cmap,
      'smap': smap,
      'fmap': fmap,
      'bmap': bmap
  })
  result = ee.Dictionary(ee.List(pvs).iterate(filter_j, first))
  return ee.Dictionary({
      'i': i.add(1),
      'alpha': alpha,
      'median': median,
      'cmap': result.get('cmap'),
      'smap': result.get('smap'),
      'fmap': result.get('fmap'),
      'bmap': result.get('bmap')
  })


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
  bmap = ee.Image.constant(ee.List.repeat(0, k.subtract(1))).add(cmap)
  alpha = ee.Image.constant(alpha)
  first = ee.Dictionary({
      'i': 1,
      'alpha': alpha,
      'median': median,
      'cmap': cmap,
      'smap': cmap,
      'fmap': cmap,
      'bmap': bmap
  })
  result = ee.Dictionary(pv_arr.iterate(filter_i, first))
  # Post-process bmap for change direction.
  bmap = ee.Image(result.get('bmap'))
  avimg = ee.Image(im_list.get(0))
  j = ee.Number(0)
  i = ee.Image.constant(1)
  first = ee.Dictionary({'avimg': avimg, 'bmap': bmap, 'j': j, 'i': i})
  dmap = ee.Dictionary(im_list.slice(1).iterate(dmap_iter, first)).get('bmap')
  return ee.Dictionary(result.set('bmap', dmap))


#   --------------------------
#    Auxiliary functions
#   -------------------------


def getS1collection(t1, t2, poly, node):
  return ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(poly).filterDate(
      ee.Date(t1), ee.Date(t2)).filter(
          ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH'])).filter(
              ee.Filter.eq('resolution_meters', 10)).filter(
                  ee.Filter.eq('instrumentMode', 'IW')).filter(
                      ee.Filter.eq('orbitProperties_pass', node)).filter(
                          ee.Filter.contains(rightValue=poly, leftField='.geo'))


def get_vvvh(image):
  """Get 'VV' and 'VH' bands from sentinel-1 imageCollection and restore linear signal from db-values."""
  return image.select('VV', 'VH').multiply(
      ee.Image.constant(math.log(10.0) / 10.0)).exp()


def convert_timestamp_list(tsl):
  """Make timestamps in YYYYMMDD format."""
  tsl = [x.replace('/', '') for x in tsl]
  tsl = ['T20' + x[4:] + x[0:4] for x in tsl]
  return tsl


def clipList(current, prev):
  """Clip a list of images and multiply by ENL."""
  imlist = ee.List(ee.Dictionary(prev).get('imlist'))
  poly = ee.Dictionary(prev).get('poly')
  ctr = ee.Number(ee.Dictionary(prev).get('ctr'))
  imlist = imlist.add(ee.Image(current).multiply(4.4).clip(poly))
  return ee.Dictionary({'imlist': imlist, 'poly': poly, 'ctr': ctr.add(1)})


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  # To debug, pass in flag `--debug=True`
  if FLAGS.debug:
    IPython.embed(user_ns=dict(globals(), **locals()))

  usage = """Usage:
  ---------------------------------------------
  Perform sequential SAR change detection

  python %s  [OPTIONS] assetname

  Options:
    -h                  this help
    -c <list>           polygon coordinates of area of interest (default: Houston AOI)
    -t <list>           time interval (default: ['2018-01-01','2019-01-01'])
    -a <float>          false positive rate (default: 0.01)
    -n <string>         node (default: 'DESCENDING')
    -o <int>            relative orbit (default: 143 )

  For a sequence of k observations, change maps are exported to GEE assets as k+3-band image:
  band 1: cmap  interval of last recorded change (0 ... k-1) byte
  band 2: fmap  interval of first recorded change (0 ... k-1) byte
  band 3: fmap  number of changes in all (0, ... k-1) byte
  bands 4 through k+3:  bmap  loewner order of changes in each interval (0, 1, 2, 3) byte     
  """ % sys.argv[0]
  options, args = getopt.getopt(sys.argv[1:], 'hc:t:a:n:o:')

  #  defaults
  coords = [[[-96.06878489327627, 29.823701939611126],
             [-96.06878489327627, 29.1423007492751],
             [-95.00860422921377, 29.1423007492751],
             [-95.00860422921377, 29.823701939611126]]]
  t1, t2 = ['2018-01-01', '2019-01-01']
  alpha = 0.01
  node = 'DESCENDING'
  orbitnum = 143

  for option, value in options:
    if option == '-h':
      print(usage)
      print(value)
      sys.exit(0)
    # djq - lint does not like eval
    # elif option == '-c':
    #   coords = eval(value)
    # elif option == '-t':
    #   t1, t2 = eval(value)
    # elif option == '-a':
    #   alpha = eval(value)
    # elif option == '-n':
    #   node = value
    # elif option == '-o':
    #   orbitnum = eval(value)

  assetpath = 'projects/ee-mortcanty/assets/contract/' + args[0]

  ee.Initialize()

  print('-------------------------------')
  print('Sequential SAR Change Detection')
  print('-------------------------------')

  try:
    poly = ee.Geometry.Polygon(coords)
    collection = getS1collection(t1, t2, poly, node)
    collection = collection.filter(
        ee.Filter.eq('relativeOrbitNumber_start',
                     orbitnum)).sort('system:time_start').filter(
                         ee.Filter.eq('platform_number', 'A'))
    count = collection.size().getInfo()
    if count < 2:
      raise ValueError('Less than 2 images found')
    print('Images found: %i' % count)
    acquisition_times = ee.List(
        collection.aggregate_array('system:time_start')).getInfo()
    timestamplist = []
    for timestamp in acquisition_times:
      tmp = time.gmtime(int(timestamp) / 1000)
      timestamplist.append(time.strftime('%x', tmp))
    timestamplist = convert_timestamp_list(timestamplist)
    print(timestamplist)
    pList = collection.map(get_vvvh).toList(500)
    first = ee.Dictionary({
        'imlist': ee.List([]),
        'poly': poly,
        'ctr': ee.Number(0)
    })
    imList = ee.List(
        ee.Dictionary(pList.iterate(clipList, first)).get('imlist'))
    # Run the algorithm ************************************************
    result = change_maps(imList, True, alpha)
    #******************************************************************
    smap = ee.Image(result.get('smap')).byte()
    cmap = ee.Image(result.get('cmap')).byte()
    fmap = ee.Image(result.get('fmap')).byte()
    bmap = ee.Image(result.get('bmap')).byte()
    #      in case of duplicates
    timestamplist1 = [
        timestamplist[i] + '_' + str(i + 1) for i in range(len(timestamplist))
    ]
    cmaps = ee.Image.cat(cmap, smap, fmap,
                         bmap).rename(['cmap', 'smap', 'fmap'] +
                                      timestamplist1[1:])
    assexport = ee.batch.Export.image.toAsset(
        cmaps.byte().clip(poly),
        description='assetExportTask',
        pyramidingPolicy={'.default': 'sample'},
        assetId=assetpath,
        scale=10,
        maxPixels=1e11)
    assexport.start()
    print('Exporting change maps to %s\n task id: %s' %
          (assetpath, str(assexport.id)))
  except Exception as e:
    print('Error: %s' % e)


if __name__ == '__main__':
  app.run(main)
