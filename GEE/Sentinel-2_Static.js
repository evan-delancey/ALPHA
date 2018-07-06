//Sentinel-2_Static
/////Description: Get Alberta wide Sentinel-2 vegetation indices and bands
/////Written by: Evan R. DeLancey - Spatial Data Scientist ABMI
/////email: edelance@ualberta.ca
/////date: 2018-07-05
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------
//load packages
//-----------------------------------------------------------------------------------------
var colorbrewer = require('users/gena/packages:colorbrewer');
//-----------------------------------------------------------------------------------------



//Cloud masking functions
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//parameters
var cloudThresh = 15;
var cloudHeights = ee.List.sequence(200,10000,250);//Height of clouds to use to project cloud shadows
var irSumThresh =0.35;//Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
var dilatePixels = 2; //Pixels to dilate around clouds
var contractPixels = 1;//Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors
//Functions
// var vizParams = {bands: ['red', 'green', 'blue'], min: 0, max: 0.3};
//////////////////////////////////////////////////////////////////////////
var rescale = function(img, exp, thresholds) {
    return img.expression(exp, {img: img})
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  };
  
  var getNotWaterClusterID = function(clusterizedImage){
  var ID = clusterizedImage.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:NotWaterPoint,
    scale:30
  });
  ID = ID.get('cluster');
  return ee.Number.parse(ID);
}
var getWaterClusterID = function(clusterizedImage){
  var ID = clusterizedImage.reduceRegion({
    reducer:ee.Reducer.mean(),
    geometry:WaterPoint,
    scale:30
  });
  ID = ID.get('cluster');
  return ee.Number.parse(ID);
}
////////////////////////////////////////
////////////////////////////////////////
// Cloud masking algorithm for Sentinel2
//Built on ideas from Landsat cloudScore algorithm
//Currently in beta and may need tweaking for individual study areas
function sentinelCloudScore(img) {
  

  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);
  
  // Clouds are reasonably bright in the blue and cirrus bands.
  score = score.min(rescale(img, 'img.blue', [0.1, 0.5]));
  score = score.min(rescale(img, 'img.cb', [0.1, 0.3]));
  score = score.min(rescale(img, 'img.cb + img.cirrus', [0.15, 0.2]));
  
  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]));

  
  //Clouds are moist
  var ndmi = img.normalizedDifference(['nir','swir1']);
  score=score.min(rescale(ndmi, 'img', [-0.1, 0.1]));
  
  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['green', 'swir1']);
  score=score.min(rescale(ndsi, 'img', [0.8, 0.6]));
  
  score = score.multiply(100).byte();
 
  return img.addBands(score.rename('cloudScore'));
}
//////////////////////////////////////////////////////////////////////////
// Function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60').int16();
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Return the masked and scaled data.
  return image.updateMask(mask);
}
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Function for finding dark outliers in time series
//Masks pixels that are dark, and dark outliers
function simpleTDOM2(c){
  var shadowSumBands = ['nir','swir1'];
  var irSumThresh = 0.4;
  var zShadowThresh = -1.2;
  //Get some pixel-wise stats for the time series
  var irStdDev = c.select(shadowSumBands).reduce(ee.Reducer.stdDev());
  var irMean = c.select(shadowSumBands).mean();
  var bandNames = ee.Image(c.first()).bandNames();
  
  //Mask out dark dark outliers
  c = c.map(function(img){
    var z = img.select(shadowSumBands).subtract(irMean).divide(irStdDev);
    var irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum());
    var m = z.lt(zShadowThresh).reduce(ee.Reducer.sum()).eq(2).and(irSum.lt(irSumThresh)).not();
    
    return img.updateMask(img.mask().and(m));
  });
  
  return c.select(bandNames);
}
////////////////////////////////////////////////////////
/////////////////////////////////////////////
/***
 * Implementation of Basic cloud shadow shift
 * 
 * Author: Gennadii Donchyts
 * License: Apache 2.0
 */
function projectShadows(cloudMask,image,cloudHeights){
  var meanAzimuth = image.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var meanZenith = image.get('MEAN_SOLAR_ZENITH_ANGLE');
  ///////////////////////////////////////////////////////
  // print('a',meanAzimuth);
  // print('z',meanZenith)
  
  //Find dark pixels
  var darkPixels = image.select(['nir','swir1','swir2']).reduce(ee.Reducer.sum()).lt(irSumThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)
  ;//.gte(1);
  
  
  //Get scale of image
  var nominalScale = cloudMask.projection().nominalScale();
  //Find where cloud shadows should be based on solar geometry
  //Convert to radians
  var azR =ee.Number(meanAzimuth).add(180).multiply(Math.PI).divide(180.0);
  var zenR  =ee.Number(meanZenith).multiply(Math.PI).divide(180.0);
  
  
 
  //Find the shadows
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    
    var shadowCastedDistance = zenR.tan().multiply(cloudHeight);//Distance shadow is cast
    var x = azR.sin().multiply(shadowCastedDistance).divide(nominalScale);//X distance of shadow
    var y = azR.cos().multiply(shadowCastedDistance).divide(nominalScale);//Y distance of shadow
    // print(x,y)
   
    return cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y));
    
    
  });
  
  
  var shadowMask = ee.ImageCollection.fromImages(shadows).max();
  // Map.addLayer(cloudMask.updateMask(cloudMask),{'min':1,'max':1,'palette':'88F'},'Cloud mask');
  // Map.addLayer(shadowMask.updateMask(shadowMask),{'min':1,'max':1,'palette':'880'},'Shadow mask');
  
  //Create shadow mask
  shadowMask = shadowMask.and(cloudMask.not());
  shadowMask = shadowMask.and(darkPixels).focal_min(contractPixels).focal_max(dilatePixels);
  
  var cloudShadowMask = shadowMask.or(cloudMask);
  
  image = image.updateMask(cloudShadowMask.not()).addBands(shadowMask.rename(['cloudShadowMask']));
  return image;
}
//////////////////////////////////////////////////////
//Function to bust clouds from S2 image
function bustClouds(img){
  img = sentinelCloudScore(img);
  img = img.updateMask(img.select(['cloudScore']).gt(cloudThresh).focal_min(contractPixels).focal_max(dilatePixels).not());
  return img;
}
//////////////////////////////////////////////////////
//Function for wrapping the entire process to be applied across collection
function wrapIt(img){
  img = sentinelCloudScore(img);
  var cloudMask = img.select(['cloudScore']).gt(cloudThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)

  img = projectShadows(cloudMask,img,cloudHeights);

  return img.clip(geometry3);
}
//////////////////////////////////////////////////////
//Function to find unique values of a field in a collection
function uniqueValues(collection,field){
    var values  =ee.Dictionary(collection.reduceColumns(ee.Reducer.frequencyHistogram(),[field]).get('histogram')).keys();
    
    return values;
  }
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
//Define Alberta shpae and projection file
//-----------------------------------------------------------------------------------------
var PrjFile = ee.Image('users/abmigc/ProcessingUnits');
var prj = PrjFile.projection();
var AB = ee.FeatureCollection('users/abmigc/Alberta_prjGEE');
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//Get static Sentinel-2 image collections
//-----------------------------------------------------------------------------------------
var S2_1 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2016-05-15', '2016-08-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2_2 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2017-05-15', '2017-08-31')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);

var S2_3 = ee.ImageCollection('COPERNICUS/S2')
  .filterDate('2018-05-15', '2018-06-29')
  .filterBounds(AB)
  .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 50)
  .map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
    return out;
    })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);  
    
var S2 = ee.ImageCollection(S2_1.merge(S2_2));
var S2 = ee.ImageCollection(S2.merge(S2_3));


//-----------------------------------------------------------------------------------------
//Apply cloud masking algorithms
//-----------------------------------------------------------------------------------------
//Bust clouds using BQA method
var S2QA = S2.map(maskS2clouds);
//Bust clouds using cloudScore method
var S2 = S2.map(bustClouds);
//Bust clouds using cloudScore and shadows using TDOM
var S2 = simpleTDOM2(S2);


//get median of all bands
var S2 = S2.median();



//-----------------------------------------------------------------------------------------
//get bands as varibles
//-----------------------------------------------------------------------------------------
var B2 = S2.select(['blue']);
var B3 = S2.select(['green']);
var B4 = S2.select(['red']);
var B5 = S2.select(['re1']);
var B6 = S2.select(['re2']);
var B7 = S2.select(['re3']);
var B8 = S2.select(['nir']);
var B8A = S2.select(['nir2']);
var B11 = S2.select(['swir1']);
var B12 = S2.select(['swir2']);
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//calcuate vegetation indices
//-----------------------------------------------------------------------------------------
var NDVI = S2.normalizedDifference(['nir', 'red']);
var NDVI705 = S2.normalizedDifference(['re2', 're1']);
var NDWI1 = S2.normalizedDifference(['green', 'nir']);
var NDWI2 = S2.normalizedDifference(['nir', 'swir1']);
var EVI = S2.expression(
  '((2.5*(B8 - B4))/((B8 + (6*B4)) - (7.5*B2) + 1))', {
    'B8': S2.select(['nir']),
    'B4': S2.select(['red']),
    'B2': S2.select(['blue'])
  }  
);
var ARI = S2.expression(
      '(B8 / B2) - (B8 / B3)', {
        'B8': S2.select(['nir']),
        'B2': S2.select(['blue']),
        'B3': S2.select(['green'])
      }
);
var REIP = S2.expression(
  '705 + 35*((((RED + RE3)/2) - RE1) / (RE2 - RE1))', {
      'RE1': S2.select(['re1']),
      'RE2': S2.select(['re2']),
      'RE3': S2.select(['re3']),
      'RED' : S2.select(['red'])
  }
);

var NBR = S2.normalizedDifference(['nir', 'swir2']);
var IRECI = S2.expression(
  '(B7 - B4)/(B5/B6)', {
    'B7': S2.select(['re3']),
    'B4': S2.select(['red']),
    'B5': S2.select(['re1']),
    'B6': S2.select(['re2'])
  }  
);// From Frampton et al 2013
var TCBgreen = S2.expression(
  '(-0.2848*B2)+(-0.2435*B3)+(-0.5436*B4)+(0.7243*B8)+(0.0840*B11)+(-0.1800*B12)',  {
    'B2': S2.select(['blue']),
    'B3': S2.select(['green']),
    'B4': S2.select(['red']),
    'B8': S2.select(['nir']),
    'B11': S2.select(['swir1']),
    'B12': S2.select(['swir2']),
  }
);
var TCBwet = S2.expression(
  '(0.1509*B2)+(0.1973*B3)+(0.3279*B4)+(0.3406*B8)+(-0.7112*B11)+(-0.4572*B12)',  {
    'B2': S2.select(['blue']),
    'B3': S2.select(['green']),
    'B4': S2.select(['red']),
    'B8': S2.select(['nir']),
    'B11': S2.select(['swir1']),
    'B12': S2.select(['swir2']),
  }
);
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//Map vegetation indices
//-----------------------------------------------------------------------------------------
/*
Map.addLayer(NDVI, {min:-0.7, max:0.75, palette: colorbrewer.Palettes.RdYlGn[11]}, 'NDVI');
Map.addLayer(NDVI705, {min:-0.7, max:0.75, palette: colorbrewer.Palettes.RdYlGn[11]}, 'NDVI705');
Map.addLayer(NDWI1, {min:-0.7, max:0.75, palette: colorbrewer.Palettes.Blues[9]}, 'NDWI1');
Map.addLayer(NDWI2, {min:-0.7, max:0.75, palette: colorbrewer.Palettes.Blues[9]}, 'NDWI2');
Map.addLayer(EVI, {min:-0.7, max:0.75, palette: colorbrewer.Palettes.RdYlGn[9]}, 'EVI');
Map.addLayer(ARI, {min:-1.5, max:0.5, palette: colorbrewer.Palettes.PRGn[11]}, 'ARI');
Map.addLayer(REIP, {min:715, max:735, palette: colorbrewer.Palettes.RdBu[10]}, 'REIP');
Map.addLayer(NBR, {min:-1, max:1, palette: colorbrewer.Palettes.BrBG[10]}, 'NBR');
Map.addLayer(IRECI, {min:-1, max:1, palette: colorbrewer.Palettes.RdYlGn[10]}, 'IRECI');
*/
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//reproject indices and bands
//-----------------------------------------------------------------------------------------
var NDVI = NDVI.reproject(prj, null, 10);
var NDVI705 = NDVI705.reproject(prj, null, 10);
var NDWI1 = NDWI1.reproject(prj, null, 10);
var NDWI2 = NDWI2.reproject(prj, null, 10);
var EVI = EVI.reproject(prj, null, 10);
var ARI = ARI.reproject(prj, null, 10);
var REIP = REIP.reproject(prj, null, 10);
var NBR = NBR.reproject(prj, null, 10);
var IRECI = IRECI.reproject(prj, null, 10);
var TCBgreen = TCBgreen.reproject(prj, null, 10);
var TCBwet = TCBwet.reproject(prj, null, 10);
var B2 = B2.reproject(prj, null, 10);
var B3 = B3.reproject(prj, null, 10);
var B4 = B4.reproject(prj, null, 10);
var B5 = B5.reproject(prj, null, 10);
var B6 = B6.reproject(prj, null, 10);
var B7 = B7.reproject(prj, null, 10);
var B8 = B8.reproject(prj, null, 10);
var B8A = B8A.reproject(prj, null, 10);
var B11 = B11.reproject(prj, null, 10);
var B12 = B12.reproject(prj, null, 10);
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//clip indices and bands to Alberta
//-----------------------------------------------------------------------------------------
var NDVI = NDVI.clip(AB);
var NDVI705 = NDVI705.clip(AB);
var NDWI1 = NDWI1.clip(AB);
var NDWI2 = NDWI2.clip(AB);
var EVI = EVI.clip(AB);
var ARI = ARI.clip(AB);
var REIP = REIP.clip(AB);
var NBR = NBR.clip(AB);
var IRECI = IRECI.clip(AB);
var TCBgreen = TCBgreen.clip(AB);
var TCBwet = TCBwet.clip(AB);
var B2 = B2.clip(AB);
var B3 = B3.clip(AB);
var B4 = B4.clip(AB);
var B5 = B5.clip(AB);
var B6 = B6.clip(AB);
var B6 = B6.clip(AB);
var B7 = B7.clip(AB);
var B8 = B8.clip(AB);
var B8A = B8A.clip(AB);
var B11 = B11.clip(AB);
var B12 = B12.clip(AB);
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//clip high and low outlier values
//-----------------------------------------------------------------------------------------
var ARI = ARI.where(ARI.gt(2), 2);
var ARI = ARI.where(ARI.lt(-2.5), -2.5);
var REIP = REIP.where(REIP.gt(735), 735);
var REIP = REIP.where(REIP.lt(715), 715);
var EVI = EVI.where(EVI.lt(-1), -1);
var EVI = EVI.where(EVI.gt(1), 1);
//-----------------------------------------------------------------------------------------

//Histogram of indices
//print(ui.Chart.image.histogram(TCBwet, geometry, 10));



//-----------------------------------------------------------------------------------------
//export vegetation indices and bands to drive
//-----------------------------------------------------------------------------------------
Export.image.toDrive({
  image: NDVI,
  description: 'NDVI',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: NDVI705,
  description: 'NDVI705',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: NDWI1,
  description: 'NDWI1',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: NDWI2,
  description: 'NDWI2',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: EVI,
  description: 'EVI',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: ARI,
  description: 'ARI',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: REIP,
  description: 'REIP',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: NBR,
  description: 'NBR',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: IRECI,
  description: 'IRECI',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: TCBgreen,
  description: 'TCBgreen',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: TCBwet,
  description: 'TCBwet',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B2,
  description: 'B2',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B3,
  description: 'B3',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B4,
  description: 'B4',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B5,
  description: 'B5',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B6,
  description: 'B6',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B7,
  description: 'B7',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});
Export.image.toDrive({
  image: B8,
  description: 'B8',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B8A,
  description: 'B8A',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B11,
  description: 'B11',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});

Export.image.toDrive({
  image: B12,
  description: 'B12',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 10E10
});
//-----------------------------------------------------------------------------------------
//END







