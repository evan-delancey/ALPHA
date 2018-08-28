////////////////////////////////////////////////////////////////////
/------------------------------------------------------------------/
/------------------------------------------------------------------/
/WETLAND MAPPING IN GEE USING RANDON FOREST CLASSIFIER			   /
/Filename:" WetlandBinary_RF_Dataset7"				    		   /
/Written and developed by Alex O. Onojeghuo (PhD)				   /
/Alberta Biodiversity Monitoring Institute, Aug, 2 2018			   /
/------------------------------------------------------------------/
/------------------------------------------------------------------/
////////////////////////////////////////////////////////////////////

//1. Define Study area and clip area

var SA = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/SA'); 
var Clip1 = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/Clip1');
var Clip2 = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/Clip2');
var Clip3 = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/Clip3'); 

Map.addLayer(SA, {}, 'Study Area');   

//2. Function to define training layer - binary training data
var train = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/ABMI_HFI2016_Binary');
var train = train.select(['b1'])
var Apal = ['#eae8e6', '#6a82a6', '#faff67'];

//3. Specify input training data
var ABMI_Plots = ee.FeatureCollection('users/Alberta_Sentinel_analysis_v20/grassland_input/GrassParkland_Plots');

//4. Stack 14 input layers for RF image classification
var MNDWI = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/MNDWI');
var PC1 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/PC1');
var PC3_10 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/PC3_10');
var PC2 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/PC2');
var HAND1 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/HAND1_SRTM');
var SWI = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/SWI_LIDAR');
var VBF = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/VBF_LIDAR');
var VVsd = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/VVsd');
var ARI = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/ARI');
var REIP = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/REIP');
var TPI500 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/TPI500_LIDAR');
var Ent_VH = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/Entropy_VHmean');
var C10 = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/C10_LIDAR');
var TRI = ee.Image('users/Alberta_Sentinel_analysis_v20/grassland_input/TRI_LIDAR');

//print(train); 

//Image stack of onput layers
var ImgStack = ee.Image([MNDWI, PC1, PC3_10, PC2, TRI, VVsd, C10, HAND1, SWI, VBF, ARI]);

print(ImgStack);

//5. Function to create stratified sample points for image classification training / validation
var points = ImgStack.addBands(train).stratifiedSample({
  numPoints: 3000, 
  classBand: "b1_11", 
  region: ABMI_Plots, 
  scale: 10,
  tileScale: 3
}).randomColumn();  

var training = points.filter(ee.Filter.lt('random', 0.5));
var validation = points.filter(ee.Filter.gte('random', 0.5));

//print(ImgStack.bandNames());



//6. Apply Random Forest (RF) algorithim to train classifier
var classifier = ee.Classifier.randomForest({numberOfTrees: 500}).train(training, "b1_11", ImgStack.bandNames());

//7. Accuracy assessment
//Get confusion matrix representing training accuracy
var trainAccuracy = classifier.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());

// Get a confusion matrix representing expected accuracy.
var Validated = validation.classify(classifier);
var testAccuracy = Validated.errorMatrix('b1_11', 'classification');
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
print('Validation consumer accuracy: ', testAccuracy.consumersAccuracy());
print('Validation producer accuracy: ', testAccuracy.producersAccuracy());
print('Kappa Coefficient: ', testAccuracy.kappa());

//8. Clip classified output to area of interest
var pred = ImgStack.classify(classifier);
var pred_clip = pred.clip(SA);
var pred_clip1 = pred.clip(Clip1);
var pred_clip2 = pred.clip(Clip2);
var pred_clip3 = pred.clip(Clip3);

//9. Add layers to map
Map.addLayer(pred_clip, {min:0, max:2, palette:Apal}, 'Prediction');
Map.addLayer(train, {min:0, max:2, palette:Apal}, 'training data');

//---------------------------------------------------------------------------------------------//
// START EXPORT INPUT VARIABLES CREATED ABOVE
// -------------------------------------------------------------------------------------------//

Export.image.toDrive({ 
image: pred_clip1,
description: 'RF_Wetlands_Clip1',
scale: 10,
region: Clip1,
maxPixels: 3E10
});

Export.image.toDrive({ 
image: pred_clip2,
description: 'RF_Wetlands_Clip2',
scale: 10,
region: Clip2,
maxPixels: 3E10
});

Export.image.toDrive({ 
image: pred_clip3,
description: 'RF_Wetlands_Clip3',
scale: 10,
region: Clip3,
maxPixels: 3E10
});
