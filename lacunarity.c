#include "lacunarity.h"

#include "raster.h"
#include "gdal.h"


int lacunarity (char *input_raster, int band, 
        int binary, long binaryThreshold, int f3d,
        int gbox_min, int gbox_max, int gbox_step)
{
  
  long *data, *dataPtr;
  int rasterX, rasterY, i;
  int ok, g;
  double l;
  
  ok = raster_band_read_long(input_raster, band, &data, &rasterX, &rasterY);
  if (ok != 0) return 1;
  
  // Convert to binary if needed.
  if (binary){
    dataPtr = data;
    for (i = 0; i < (rasterX * rasterY); i++){
      if (*dataPtr >= binaryThreshold) *dataPtr = 1;
      else *dataPtr = 0;
      dataPtr++;
    }
  }
  
  fprintf(stdout, "Lacunarity index for %s:\n", input_raster);
  fprintf(stdout, "Gliding box size\tLacunarity index\n");
  
  for (g = gbox_min; g <= gbox_max; g += gbox_step){
    l = lacunarity_in_window(data, rasterX, rasterY, f3d, g, 0, 0, rasterX, rasterY);
    fprintf(stdout, "%i\t%f\n", g, l);
  }
  
  free(data);
  return 0;
}


int spatial_lacunarity (char *input_raster, int band, 
            int binary, long binaryThreshold, int f3d,
            int gbox, int mwin, 
            char *output_file, char *format)
{
  long *data, *dataPtr;         // The input data array.
  int rasterX, rasterY;         // The size of the input raster.
  int i, j;                     // x- and y-coordinates loop variables
  double *lacunarity;           // The lacunarity data array.
  double *lacunarityPtr;
  int outRasterX, outRasterY;   // The size of the output raster.
  double georeference[6];       // Georeference for output raster file.
  GDALDatasetH dataset;         // The GDAL dataset for the input raster file.
  int ok;
  int pctDone, curPctDone;      // Percentage done.
  
  // Read the input raster file.
  ok = raster_band_read_long(input_raster, band, &data, &rasterX, &rasterY);
  if (ok != 0){
    fprintf(stderr, "ERROR. Unable to read input raster file.\n");
    return 1;
  }
  
  // Convert to binary if needed.
  if (binary){
    dataPtr = data;
    for (i = 0; i < (rasterX * rasterY); i++){
      if (*dataPtr >= binaryThreshold) *dataPtr = 1;
      else *dataPtr = 0;
      dataPtr++;
    }
  }
  
  // Get the georeference of the input raster.
  dataset = GDALOpen(input_raster, GA_ReadOnly);
  GDALGetGeoTransform(dataset, georeference);
  
  // Create the output lacunarity data array.
  outRasterX = rasterX - mwin + 1;
  outRasterY = rasterY - mwin + 1;
  lacunarity = (double*)calloc(outRasterX * outRasterY, sizeof(double));
  if (lacunarity == NULL){
    fprintf(stderr,"ERROR. Not enough memory for creating lacunarity raster.\n");
    free(data);
    return 1;
  }
  
  // Compute the lacunarity value for each point in the lacunarity array.
  pctDone = 0;
  lacunarityPtr = lacunarity;
  for (j = 0; j < outRasterY; j++){
    for (i = 0; i < outRasterX; i++){
      *lacunarityPtr = lacunarity_in_window(
        data, rasterX, rasterY, f3d, gbox, i, j, mwin, mwin
      );
      lacunarityPtr++;
      
      curPctDone = floorl(10*(j*outRasterX + i) / (outRasterX*outRasterY));
      if (curPctDone > pctDone){
        pctDone = curPctDone;
        fprintf(stdout, "%i%% done.\n", pctDone*10);
      }
    }
  }
  fprintf(stdout, "100%% done.\n");
  fprintf(stdout, "Writing lacunarity image to file...\n");
  
  // Write the output lacunarity raster to the output file.
  // We need to provide the georeference.
  // The output raster image is smaller by mwin-1 pixels.
  georeference[0] += (mwin-1)*georeference[1];    // Shift the top left x coordinate.
  georeference[3] += (mwin-1)*georeference[5];    // Shift the top left y coordinate.
  ok = raster_band_write_double(output_file, format, 1, georeference, lacunarity, outRasterX, outRasterY);
  if (ok != 0){
    fprintf(stderr, "ERROR. Unable to write output raster file.\n");
    free(data);
    return 1;
  }
  
  free(data);
  return 0;
}







double lacunarity_in_window (
  long *data, int rasterX, int rasterY, 
  int f3d, int gbox, 
  int mwinX, int mwinY, int mwinW, int mwinH)
{

  double lacunarity;                  // The resulting lacunarity.
  long *imgPtr;                       // This will be our pointer for looping through our data.
  int i, j, k;                        // Counter variables.
  int c;                              // Temporary variables.
  int nGlidingStepsX, nGlidingStepsY; // Number of gliding box steps in each direction.
  int nGlidingBoxes;                  // Number of gliding boxes.
  long maxValue;                      // Maximum value in the whole moving window.
  long nLevels;                       // The number of levels.
  long *intensitySum;                 // Table of the sum of intensity values.
  long *intensitySumPtr;              // Pointer to the above table.
  long intensitySumValue;             // The intensity sum value.
  long maxIntensity;                  // Maximum value for intensity.
  int gbx, gby;                       // Coordinates for the gliding box (x and y).
  double *probDens;                   // Pointer to the probability density table.
  double *probDensPtr;
  double M, M2;                       // The distribution moments.
  
  // Place a pointer at the start of the data.
  imgPtr = data;
  imgPtr += (rasterX * mwinY) + mwinX;
  
  // Compute the number of levels.
  // For computing k, we need to find the maximum value in the
  // moving window. We have to loop through all pixels in the moving window.
  maxValue = 0;
  for (j = 0; j < mwinH; j++){
    for (i = 0; i < mwinW; i++){
      if (maxValue < *imgPtr) maxValue = *imgPtr;
      imgPtr++;
    }
    
    // Go to the end of the line:      imgPtr += (rasterX - mwinX - mwinW);
    // Go to the start of the next line:  imgPtr += mwinX;
    imgPtr += (rasterX - mwinW);
  }
  
  // If the maximum value is 0, we do not need to compute the
  // lacunarity for this gliding box. It is equal to 0.
  if (maxValue <= 0) return 0.0f;
  
  // Compute the number of levels.
  if (f3d)
    nLevels = maxValue;
  else
    nLevels = lrint(ceil((double)maxValue / (double)gbox));
  
  // Loops for gliding our box through our local window.
  nGlidingStepsX = mwinW - gbox + 1;
  nGlidingStepsY = mwinH - gbox + 1;

  // The real number of gliding steps is nGlidingStepsX * nGlidingStepsY.
  // We will create the table containing the sum of all intensity values
  // for each cube box. This table corresponds to the graphic 5e, p. 512,
  // in Myint (2005).
  // We need to allocate the necessary memory for this table. The size of 
  // the table is (number of levels) by (number of gliding boxes).
  nGlidingBoxes = nGlidingStepsX * nGlidingStepsY;
  intensitySum = calloc((nLevels * nGlidingBoxes), sizeof(int));
  if (intensitySum == NULL){
    fprintf(stderr, "ERROR. Not enough memory for summing up the intensity values.\n");
    return 0;
  }
  
  // Compute the intensity sum table.
  // Loop in y direction.
  for (j = 0; j < nGlidingStepsY; j++){
    // Loop in x direction; loops faster.
    for (i = 0; i < nGlidingStepsX; i++){
      
      // Place the image pointer at the start of the data.
      imgPtr = data;
      imgPtr += ((mwinY+j)*rasterX) + (mwinX+i);
      
      // Place the pointer at the start of the intensity sum table.
      intensitySumPtr = intensitySum + (((j * nGlidingStepsX) + i)  * nLevels);
      
      // Loop through all cells of the gliding box.
      if (f3d){
        // Here comes the loop for the true 3D gliding box. In this case, the number of levels
        // is the same as the maximum value, and for each step, we sum up all the values 
        // as in 2D case.
        for (gby = 0; gby < gbox; gby++){
          for (gbx = 0; gbx < gbox; gbx++){
            c = *imgPtr;
            for (k = 0; k < nLevels; k++){
              *intensitySumPtr += MIN((c-k), gbox);
              intensitySumPtr++;
            }
            imgPtr++;
            intensitySumPtr -= nLevels;
          }
          // Go to the next line in the image.
          imgPtr += (rasterX - gbox);
        }
      }else{
        // The loop for the layered gliding box as described by Myint and Lam (2005).
        for (gby = 0; gby < gbox; gby++){
          for (gbx = 0; gbx < gbox; gbx++){
            c = *imgPtr;            // The value of the pixel.
            for (k = 0; k < nLevels; k++){    // Loop through all levels.
              if (c >= gbox){
                *intensitySumPtr += gbox;
                c -= gbox;
              }else{
                *intensitySumPtr += c;
                c = 0;
              }
              intensitySumPtr++;
            }
            imgPtr++;
            intensitySumPtr -= nLevels;
          }
          // Go to the next line in the image.
          imgPtr += (rasterX - gbox);
        }
      }
    }  // for (unsigned short i = 0; i < nGlidingSteps; i++)
  }  // for (unsigned short j = 0; j < nGlidingSteps; j++)
  
  
//  // --- DEBUG ---
//  FILE *pFile = fopen("/Temp/lacunarity_intensity.txt", "w");
//  intensitySumPtr = intensitySum;
//  for (i = 0; i < (nLevels * nGlidingBoxes); i++)
//  {
//    int intensity = *intensitySumPtr;
//    fprintf(pFile, "%i\n", intensity);
//    intensitySumPtr++;
//  }
//  fclose(pFile);  
//  // --- END DEBUG ---
  
  
  // Compute the probability density for the intensity sum table.
  
  // Find first the maximum intensity value.
  maxIntensity = 0;
  intensitySumPtr = intensitySum;
  for (i = 0; i < (nGlidingBoxes * nLevels); i++){
    intensitySumValue = *intensitySumPtr;
    if (maxIntensity < intensitySumValue){
      maxIntensity = intensitySumValue;
    }
    intensitySumPtr++;
  }
  // Allocate a table for the probability density values.
  probDens = calloc((maxIntensity+1), sizeof(double));
  if (probDens == NULL){
    fprintf(stderr, "ERROR. Not enough memory to compute probability density values.\n");
    return 0;
  }
  
  intensitySumPtr = intensitySum;
  
  // Fill the probability density table.
  for (i = 0; i < (nGlidingBoxes * nLevels); i++){
    probDensPtr = probDens + *intensitySumPtr;
    double pd = *probDensPtr + 1.0;
    *probDensPtr = pd;
    intensitySumPtr++;
  }  
  // Normalisation of the probability density table.
  probDensPtr = probDens;
  for (i = 0; i < (maxIntensity + 1); i++){
    *probDensPtr = *probDensPtr / (nGlidingBoxes * nLevels);
    probDensPtr++;
  }
  
//  // --- DEBUG ---
//  FILE *pFile = fopen("/Temp/lacunarity.txt", "w");
//  probDensPtr = probDens;
//  for (i = 0; i < (maxIntensity+1); i++)
//  {
//    double probability = *probDensPtr;
//    fprintf(pFile, "%f\n", probability);
//    probDensPtr++;
//  }
//  fclose(pFile);  
//  // --- END DEBUG ---  
  
  // Compute the distribution moments.
  M = 0;
  M2 = 0;
  probDensPtr = probDens;
  for (i = 0; i < (maxIntensity+1); i++){
    //double probability, MTemp;
    //probability = *probDensPtr;
    //MTemp = i * probability;
    M += i * (*probDensPtr);
    //M += MTemp;
    M2 += (i * i) * (*probDensPtr);
    probDensPtr++;
  }
  
  // Compute the lacunarity index.
  lacunarity = (M2 / (M * M));
  
  //  // --- DEBUG ---
  //  FILE *pFile = fopen("/lacunarity.txt", "w");
  //  fprintf(pFile, "M: %f\n", M);
  //  fprintf(pFile, "M2: %f\n", M2);
  //  fprintf(pFile, "Lacunarity: %f\n", lacunarity);
  //  fclose(pFile);  
  //  // --- END DEBUG ---  
  free(intensitySum);
  free(probDens);
  
  return lacunarity;
}
