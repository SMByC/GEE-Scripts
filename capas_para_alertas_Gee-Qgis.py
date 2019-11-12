import ee
ee.Initialize()

from ee_plugin import Map
from ee_plugin import utils

### datos para modificar por corrida:
# menu de ingreso de fechas a comparar
startDate1='2019-07-15'
endDate1= '2019-08-27'

startDate2='2019-08-28'
endDate2='2019-09-05'

### Otros parametros
bandNamesLandsatTOA = ee.List(['blue','green','red','nir','swir1','swir2','pixel_qa'])
sensorBandDictLandsatTOA = ee.Dictionary({'L8': ee.List([1,2,3,4,5,6,10]), 'L7': ee.List([0,1,2,3,4,6,9])})

vizParams = {"bands":["B8","B11","B12"],"max":0.4,"gamma":1.6}
imageVisParam = {"opacity":1,"bands":["VV_last"],"min":-25,"max":0,"gamma":1}

# sentinel
s2Bands = ["B2", "B3", "B4", "B8", "B11", "B12", "QA60", "date"]

### funtions

# landsat 8RT cloudmask
def mascara_l8(image):
    mask_8 = image.select('BQA')
    cloudAndShadowb = mask_8.neq(6896).And(mask_8.neq(2800))
    return image.updateMask(image.mask()).And(cloudAndShadowb)

# landsat 7RT cloudmask
def mascara_l7(image):
    mask_7 = image.select('BQA')
    cloudAndShadow = mask_7.neq(760).And(mask_7.neq(764))
    return image.updateMask(image.mask().And(cloudAndShadow))

# landsat 8SR cloudmask
def maskL8SR(image):
    # Bits 3 and 5 are cloud shadow and cloud, respectively.
    cloudShadowBitMask = ee.Number(2).pow(3).int()
    cloudsBitMask = ee.Number(2).pow(5).int()
    # Get the QA band.
    qa = image.select('pixel_qa')
    # Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))
    date = ee.Number.parse(image.date().format("YYYYMMdd"))
    return image\
     .addBands(ee.Image.constant(date).rename("date").double()).updateMask(mask)
    
# landsat 7SR cloudmask
def maskL57SR(image):
    qa = image.select('pixel_qa')
    # Second bit must be zero, meaning none to low cloud confidence.
    mask1 = qa.bitwiseAnd(ee.Number(2).pow(7).int()).eq(0).And( \
        qa.bitwiseAnd(ee.Number(2).pow(3).int()).lte(0)) # cloud shadow
    # This gets rid of irritating fixed-pattern noise at the edge of the images.
    # mask2 = image.select('B.*').gt(0).reduce('min')
    date = ee.Number.parse(image.date().format("YYYYMMdd"))
    return image \
        .addBands(ee.Image.constant(date).rename("date").double()) \
        .updateMask(mask1) #.and(mask2))

def maskS2clouds(image):
  qa = image.select('QA60')
  # Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudBitMask).eq(0)
  #.and(
  #           qa.bitwiseAnd(cirrusBitMask).eq(0));
  date = ee.Number.parse(image.date().format("YYYYMMdd"))

  return image.updateMask(mask).divide(10000).float() \
        .addBands(ee.Image.constant(date).rename("date").double())

def maskS2clouds2(image2):
  qa2 = image2.select('QA60')
  # Both flags should be set to zero, indicating clear conditions.
  mask2 = qa2.bitwiseAnd(cloudBitMask).eq(0)
  #.and(
  #           qa.bitwiseAnd(cirrusBitMask).eq(0));
  date = ee.Number.parse(image2.date().format("YYYYMMdd"))

  return image2.updateMask(mask2).divide(10000).float() \
        .addBands(ee.Image.constant(date).rename("date").double())

### sentinel

sent_antes = ee.ImageCollection('COPERNICUS/S2') \
    .filterDate(startDate1,endDate1)

sent_despues = ee.ImageCollection('COPERNICUS/S2') \
    .filterDate(startDate2,endDate2)

# Bits 10 and 11 are clouds and cirrus, respectively.
cloudBitMask = ee.Number(2).pow(10).int()
cirrusBitMask = ee.Number(2).pow(11).int()

# compuestos de mediana

sentinel_antes = sent_antes.map(maskS2clouds)
sentinel_despues = sent_despues.map(maskS2clouds2)

median2 = sentinel_antes.median()
sentinel_antes_med = median2 #.clip(AOI)

mediand2 = sentinel_despues.median()
sentinel_despues_med = mediand2 #.clip(AOI)

###

l8s = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') \
    .filterDate(startDate1,endDate1) \
    .select(sensorBandDictLandsatTOA.get('L8'),bandNamesLandsatTOA) \
    .map(maskL8SR)

l7s = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR') \
    .filterDate(startDate1,endDate1) \
    .select(sensorBandDictLandsatTOA.get('L7'),bandNamesLandsatTOA) \
    .map(maskL57SR)
          
ls_antes = ee.ImageCollection(l7s.merge(l8s))

# modulo de fusion de compuestos landsat
landsat1 = ls_antes.median()
LT_previo_fusion = landsat1 #.clip(AOI)

### colecciones RT_TOA individuales enmascaradas
# landsat 8 RT_TOA enmascarada
lRt8rt_toa = ee.ImageCollection('LANDSAT/LC08/C01/T1_RT_TOA') \
          .filterDate(startDate2,endDate2) \
          .map(mascara_l8)
lRt8rt_toa_median = lRt8rt_toa.reduce(ee.Reducer.median())

# landsat 7 RT_TOA enmascarada
lRt7rt_toa = ee.ImageCollection('LANDSAT/LE07/C01/T1_RT_TOA') \
          .filterDate(startDate2,endDate2) \
          .map(mascara_l7)
lRt7rt_toa_median = lRt7rt_toa.reduce(ee.Reducer.median())

### colecciones RT individuales sin enmascarar

# landsat 8
lRt8rt = ee.ImageCollection('LANDSAT/LC08/C01/T1_RT').filterDate(startDate1,endDate1)
lRt8rt_t1_median = lRt8rt.reduce(ee.Reducer.median())

# landsat 7
lRt7rt = ee.ImageCollection('LANDSAT/LE07/C01/T1_RT') \
          .filterDate(startDate2,endDate2)
lRt7rt_t1_median = lRt7rt.reduce(ee.Reducer.median())

### sentinel 1

# Load the Sentinel-1 ImageCollection.
sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')

# compuesto antes
vv = sentinel1 \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
    .filter(ee.Filter.eq('instrumentMode', 'IW')) \
    .filterDate(startDate1, endDate1)
  
# filter vv
vvAscending = vv.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
Ascending_antes = vvAscending.reduce(ee.Reducer.last())

Map.addLayer(Ascending_antes, imageVisParam, 'composite_antes_ascending', shown=False)

vv2 = sentinel1 \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
    .filter(ee.Filter.eq('instrumentMode', 'IW')) \
    .filterDate(startDate2, endDate2)

vvAscending2 = vv2.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))

Ascending_despues = vvAscending2.reduce(ee.Reducer.last())

Map.addLayer(Ascending_despues, imageVisParam, 'composite_despues_ascending', shown=False)

vv3 = sentinel1 \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
    .filter(ee.Filter.eq('instrumentMode', 'IW')) \
    .filterDate(startDate2, endDate2)

vvDescending = vv3.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
Descending = vvDescending.reduce(ee.Reducer.last()) 

Map.addLayer(Descending, imageVisParam, 'composite_despues_descending', shown=False)

vv4 = sentinel1 \
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
  .filter(ee.Filter.eq('instrumentMode', 'IW')) \
  .filterDate(startDate1, endDate1)

vvDescending2 = vv4.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
Descending_antes = vvDescending2.reduce(ee.Reducer.last())

Map.addLayer(Descending_antes, imageVisParam, 'composite_antes_descending', shown=False)

### mosaicos

Map.addLayer(lRt8rt_toa_median, {'bands': ["B5_median", "B6_median", "B4_median"], 'min':0, 'max':1, 'gamma': 1.6}, 'LC8_TOA_posterior_enmask', shown=False)  # landsat TOA enmascarado
Map.addLayer(lRt8rt_t1_median, {'bands': ["B5_median", "B6_median", "B4_median"], 'min':0, 'max':30000}, 'LC8_RT_posterior', shown=False)  # landsat RT sin enmascarar

Map.addLayer(lRt7rt_toa_median, {'bands': ["B4_median", "B5_median", "B3_median"], 'min':0, 'max':1, 'gamma': 1.6}, 'LE7_RT_TOA_posterior_enmask', shown=False)  # landsat TOA enmascarado
Map.addLayer(lRt7rt_t1_median, {'bands': ["B4_median", "B5_median", "B3_median"], 'min':0, 'max':255, 'gamma': 1.6}, 'LE7_RT_posterior', shown=False)  # landsat RT sin enmascarar

Map.addLayer(LT_previo_fusion, {'bands': ['nir', 'swir1', 'red'], 'min':0, 'max':4000, 'gamma': 1.6}, 'mediana_landsat_antes', shown=False)  # mediana landsat previo

Map.addLayer(sentinel_antes_med, {'bands': ['B8', 'B11', 'B12'], 'max': 0.6, 'gamma': 1.6}, 'mediana_sentinel_antes', shown=False)  # mediana sentinel2 previo
Map.addLayer(sentinel_despues_med, {'bands': ['B8', 'B11', 'B12'], 'max': 0.6, 'gamma': 1.6}, 'mediana_sentinel_despues', shown=False)  # mediana sentinel2 posterior




