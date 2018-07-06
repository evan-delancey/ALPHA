/////////////////////////////////////////////////////////////////////
//-------------------Sentinel-1_Static------------------------------/
//Written by:                                                       /
////Evan R. DeLancey                                                /
////Spatial Data Scientist Alberta Biodiversity Monitoring Institute/
////2018-05-31                                                      /
////edelance@ualberta.ca                                            /
//                                                                  /
//Description:                                                      /
////-Get static Sentinel-1 Varibles for the province of Alberta     /
////-Use S1 multi temporal filtering and angle correction to get    /
//// VV, VH, and DPOL varibles                                      /
/////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//Define Functions
///////////////////////////////////////////////////////////////////////////////
function S1Prep(image) {
  var first = image.addBands(image.select('VV').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VV_Gamma'));
  var second = first.addBands(image.select('VH').subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VH_Gamma'));
  return second.addBands(second.select(['VH']).divide(second.select(['VV'])).rename('DPOL'))
}

function multitemporalDespeckle(images, radius, units, opt_timeWindow) {
  var timeWindow = opt_timeWindow || { before: -3, after: 3, units: 'month' }
  
  var bandNames = ee.Image(images.first()).bandNames()
  var bandNamesMean = bandNames.map(function(b) { return ee.String(b).cat('_mean') })
  var bandNamesRatio = bandNames.map(function(b) { return ee.String(b).cat('_ratio') })
  
  // compute space-average for all images
  var meanSpace = images.map(function(i) {
    var reducer = ee.Reducer.mean()
    var kernel = ee.Kernel.square(radius, units)
    
    var mean = i.reduceNeighborhood(reducer, kernel).rename(bandNamesMean)
    var ratio = i.divide(mean).rename(bandNamesRatio)

    return i.addBands(mean).addBands(ratio)
  })

  /***
   * computes a multi-temporal despeckle function for a single image
   */
  function multitemporalDespeckleSingle(image) {
    var t = image.date()
    var from = t.advance(ee.Number(timeWindow.before), timeWindow.units)
    var to = t.advance(ee.Number(timeWindow.after), timeWindow.units)
    
    var meanSpace2 = ee.ImageCollection(meanSpace).select(bandNamesRatio).filterDate(from, to)
      .filter(ee.Filter.eq('relativeOrbitNumber_start', image.get('relativeOrbitNumber_start'))) // use only images from the same cycle
    
    var b = image.select(bandNamesMean)

    return b.multiply(meanSpace2.sum()).divide(meanSpace2.count()).rename(bandNames)
  }
  
  return meanSpace.map(multitemporalDespeckleSingle).select(bandNames)
}

function maskAngle(image) {
  var ang = image.select(['angle']);
  var first = image.updateMask(ang.gt(30.63993));
  return first.updateMask(ang.lt(44.73993));
}

function maskEdge(img) {
  var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
  return img.updateMask(mask.select(0));  
}
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////

 

///////////////////////////////////////////////////////////////////////////////
//Define Alberta shpae and projection file
///////////////////////////////////////////////////////////////////////////////
var PrjFile = ee.Image('users/abmigc/ProcessingUnits');
var prj = PrjFile.projection();
var AB = ee.FeatureCollection('users/abmigc/Alberta_prjGEE');

///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
//Get static Sentinel-1 varibles (VV, VH, DPOL)
///////////////////////////////////////////////////////////////////////////////
var S1_1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2017-05-15', '2017-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
  
var S1_2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2018-05-15', '2018-06-20')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
  
var S1_3 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2016-05-15', '2016-08-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);
var S1 = ee.ImageCollection(S1_1.merge(S1_2));
var S1 = ee.ImageCollection(S1.merge(S1_3));

//Prep S1 images
var S1 = S1.map(maskEdge);
var S1 = S1.map(maskAngle);
var S1 = S1.map(S1Prep);
var VV = S1.select(['VV_Gamma']);
var VH = S1.select(['VH_Gamma']);

//Set parameters for multi temporal filtering (MTF)
var radius = 7; 
var units = 'pixels';

//Apply MTF
var VV = multitemporalDespeckle(VV, radius, units, { before: -12, after: 12, units: 'month'});
var VVsd = VV.reduce(ee.Reducer.stdDev());
var VV = VV.mean();
Map.addLayer(VV, {min:-25, max:-5}, 'VV');

var VH = multitemporalDespeckle(VH, radius, units, { before: -12, after: 12, units: 'month'});
var VHsd = VH.reduce(ee.Reducer.stdDev());
var VH = VH.mean();
Map.addLayer(VH, {min:-25, max:-5}, 'VH');

//Generate difference in polarization
var DPOL = VH.divide(VV);
Map.addLayer(DPOL, {min:1, max:5}, 'DPOL');
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////
//Export varibles
///////////////////////////////////////////////////////////////////////////////
var VV = VV.reproject(prj, null, 10);
var VV = VV.clip(AB);

var VVsd = VVsd.reproject(prj, null, 10);
var VVsd = VVsd.clip(AB);

var VHsd = VHsd.reproject(prj, null, 10);
var VHsd = VHsd.clip(AB);

var VH = VH.reproject(prj, null, 10);
var VH = VH.clip(AB);

var DPOL = DPOL.reproject(prj, null, 10);
var DPOL = DPOL.clip(AB);

Export.image.toDrive({
  image: VV,
  description: 'VV',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: VVsd,
  description: 'VVsd',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: VHsd,
  description: 'VHsd',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: VH,
  description: 'VH',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: DPOL,
  description: 'DPOL',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//Get Summer - fall Sentinel-1 variables
///////////////////////////////////////////////////////////////////////////////
var S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AB)
  .filterDate('2017-04-01', '2017-12-31')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV', 'VH'])
  .filterMetadata('resolution_meters', 'equals' , 10);



//Prep S1 images
var S1 = S1.map(maskAngle);
var S1 = S1.map(S1Prep);
var VV = S1.select(['VV_Gamma']);
var VH = S1.select(['VH_Gamma']);

//Get summer and fall VV and VH
var VV_summer = VV.filterDate('2017-06-15', '2017-08-01');
print(VV_summer);
var VV_fall = VV.filterDate('2017-09-15', '2017-10-31');
print(VV_fall);
var VH_summer = VH.filterDate('2017-06-15', '2017-08-01');
var VH_fall = VH.filterDate('2017-09-15', '2017-10-31');

//Set parameters for multi temporal filtering (MTF)
var radius = 7; 
var units = 'pixels';

//Apply MTF
var VV_summer = multitemporalDespeckle(VV_summer, radius, units, { before: -12, after: 12, units: 'month'});
var VV_summer = VV_summer.median();

var VV_fall = multitemporalDespeckle(VV_fall, radius, units, { before: -12, after: 12, units: 'month'});
var VV_fall = VV_fall.median();

var VH_summer = multitemporalDespeckle(VH_summer, radius, units, { before: -12, after: 12, units: 'month'});
var VH_summer = VH_summer.median();

var VH_fall = multitemporalDespeckle(VH_fall, radius, units, { before: -12, after: 12, units: 'month'});
var VH_fall = VH_fall.median();

//Generate summer - fall
var VVdiff = VV_summer.subtract(VV_fall);
var VHdiff = VH_summer.subtract(VH_fall);

var DPOL_summer = VH_summer.divide(VV_summer);
var DPOL_fall = VH_fall.divide(VV_fall);
var DPOLdiff = DPOL_summer.subtract(DPOL_fall);

Map.addLayer(VVdiff, {min:-10, max:10}, 'VV difference');
Map.addLayer(VHdiff, {min:-10, max:10}, 'VH difference');
Map.addLayer(DPOLdiff, {min:-2, max:2}, 'DPOL difference');


///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////
//Export varibles
///////////////////////////////////////////////////////////////////////////////
var VVdiff = VVdiff.reproject(prj, null, 10);
var VVdiff = VVdiff.clip(AB);

var VHdiff = VHdiff.reproject(prj, null, 10);
var VHdiff = VHdiff.clip(AB);

var DPOLdiff = DPOLdiff.reproject(prj, null, 10);
var DPOLdiff = DPOLdiff.clip(AB);

Export.image.toDrive({
  image: VVdiff,
  description: 'VVdiff',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: VHdiff,
  description: 'VHdiff',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});

Export.image.toDrive({
  image: DPOLdiff,
  description: 'DPOLdiff',
  scale: 10,
  region: AB,
  folder: 'Alberta_EO_Data',
  maxPixels: 3E10
});
///////////////////////////////////////////////////////////////////////////////
//END
///////////////////////////////////////////////////////////////////////////////
