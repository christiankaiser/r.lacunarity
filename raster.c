

#include "raster.h"


#include <stdlib.h>
#include <stdio.h>



int raster_band_read_double (char *raster, int band, double **data, int *rasterX, int *rasterY)
{
  GDALDatasetH idataset;
  GDALRasterBandH iband;        // The input raster band.
  
  
  // Open the input raster file.
  idataset = GDALOpen(raster, GA_ReadOnly);
  if (idataset == NULL)
  {
    fprintf(stderr, "Error. Unable to open raster '%s'\n", raster);
    return 1;
  }
  
  // Get the size of the input raster.
  *rasterX = GDALGetRasterXSize(idataset);
  *rasterY = GDALGetRasterYSize(idataset);
  
  // Get the input raster band.
  iband = GDALGetRasterBand(idataset, band);
  
  // Fetch the input raster band content.
  *data = (double*) malloc(*rasterX * *rasterY * sizeof(double));
  if (*data == NULL)
  {
    GDALClose(idataset);
    fprintf(stderr, "Error. Not enough memory to read raster '%s'.\n", raster);
    return 1;
  }
  GDALRasterIO(iband, GF_Read, 0, 0, *rasterX, *rasterY, *data, *rasterX, *rasterY, GDT_Float64, 0, 0);
  
  GDALClose(idataset);
  
  return 0;
}




int raster_band_read_long (char *raster, int band, long **data, int *rasterX, int *rasterY)
{
  GDALDatasetH idataset;
  GDALRasterBandH iband;        // The input raster band.
  
  
  // Open the input raster file.
  idataset = GDALOpen(raster, GA_ReadOnly);
  if (idataset == NULL)
  {
    fprintf(stderr, "Error. Unable to open raster '%s'\n", raster);
    return 1;
  }
  
  // Get the size of the input raster.
  *rasterX = GDALGetRasterXSize(idataset);
  *rasterY = GDALGetRasterYSize(idataset);
  
  // Get the input raster band.
  iband = GDALGetRasterBand(idataset, band);
  
  // Fetch the input raster band content.
  *data = (long*) malloc(*rasterX * *rasterY * sizeof(long));
  if (*data == NULL)
  {
    GDALClose(idataset);
    fprintf(stderr, "Error. Not enough memory to read raster '%s'.\n", raster);
    return 1;
  }
  GDALRasterIO(iband, GF_Read, 0, 0, *rasterX, *rasterY, *data, *rasterX, *rasterY, GDT_Int32, 0, 0);
  
  GDALClose(idataset);
  
  return 0;
}







int raster_band_write_double(char *raster, char *format, int band, double *adfGeoTransform, 
               double *data, int rasterX, int rasterY)
{
  
  FILE *fp;
  int rasterExists;
  GDALDriverH hDriver;
  GDALDatasetH hDataset;
  GDALRasterBandH hBand;
  
  
  
  // Check first whether the file already exists. If so, we will write to the existing file.
  rasterExists = 0;
  fp = fopen(raster, "r");
  if (fp) 
  {
    rasterExists = 1;
    fclose(fp);
  }
  
  
  
  if (!rasterExists)
  {
    hDriver = GDALGetDriverByName(format);
    if (hDriver == NULL)
    {
      fprintf(stderr, "ERROR. Unable to create output raster file.\n");
      fprintf(stderr, "No driver found for raster format %s\n\n", format);
      return 1;
    }
    
    hDataset = GDALCreate(hDriver, raster, rasterX, rasterY, 1, GDT_Float64, NULL);
    
    // Set the dataset information.
    GDALSetGeoTransform(hDataset, adfGeoTransform);
    
    // Band must always be 1 when creating a new raster.
    band = 1;
  }
  else
  {
    hDataset = GDALOpen(raster, GA_Update);
  }
  if (hDataset == NULL)
  {
    fprintf(stderr, "ERROR. Unable to open output dataset.\n\n");
    return 1;
  }
  
  
  
  hBand = GDALGetRasterBand(hDataset, band);
  if (hBand == NULL)
  {
    fprintf(stderr, "ERROR. Cannot write to band %i\n\n", band);
    GDALClose(hDataset);
    return 1;
  }
  
  
  // Write the raster band.
  GDALRasterIO(hBand, GF_Write, 0, 0, rasterX, rasterY, data, rasterX, rasterY, GDT_Float64, 0, 0);
  
  
  // Close the dataset.
  GDALClose(hDataset);
  
  
  return 0;
}









void pixel_coord_to_geo(double *padfTransform, double pixelX, double pixelY, double *geoX, double *geoY)
{
  *geoX = padfTransform[0] + pixelX*padfTransform[1] + pixelY*padfTransform[2];
  *geoY = padfTransform[3] + pixelX*padfTransform[4] + pixelY*padfTransform[5];
}




void geo_coord_to_pixel(double *padfTransform, double geoX, double geoY, double *pixelX, double *pixelY)
{
  *pixelX = (geoX*padfTransform[2] - padfTransform[2]*padfTransform[3] - geoX*padfTransform[5] +
         padfTransform[0]*padfTransform[5]) / (padfTransform[2]*padfTransform[4] -
                           padfTransform[1]*padfTransform[5]);
  
  *pixelY = (geoY*padfTransform[1] - padfTransform[3]*padfTransform[1] - geoX*padfTransform[4] +
         padfTransform[0]*padfTransform[4]) / (padfTransform[1]*padfTransform[5] -
                           padfTransform[2]*padfTransform[4]);
                           
}






void raster_bbox (GDALDatasetH raster, double *bbox)
{
  int rasterSizeX, rasterSizeY;
  double x, y, minx, maxx, miny, maxy;
  double padfTransform[6];
  
  
  GDALGetGeoTransform(raster, padfTransform);
  rasterSizeX = GDALGetRasterXSize(raster);
  rasterSizeY = GDALGetRasterYSize(raster);
  
  
  pixel_coord_to_geo(padfTransform, 0.0, 0.0, &x, &y);
  minx = x;
  maxx = x;
  miny = y;
  maxy = y;
  
  pixel_coord_to_geo(padfTransform, 0.0, (double)rasterSizeY, &x, &y);
  if (x < minx) minx = x;
  if (x > maxx) maxx = x;
  if (y < miny) miny = y;
  if (y > maxy) maxy = y;
  
  pixel_coord_to_geo(padfTransform, (double)rasterSizeX, (double)rasterSizeY, &x, &y);
  if (x < minx) minx = x;
  if (x > maxx) maxx = x;
  if (y < miny) miny = y;
  if (y > maxy) maxy = y;
  
  pixel_coord_to_geo(padfTransform, (double)rasterSizeX, 0.0, &x, &y);
  if (x < minx) minx = x;
  if (x > maxx) maxx = x;
  if (y < miny) miny = y;
  if (y > maxy) maxy = y;
  
  bbox[0] = minx;
  bbox[1] = miny;
  bbox[2] = maxx;
  bbox[3] = maxy;
  
}










