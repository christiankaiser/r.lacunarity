
#include <GDAL/gdal.h>



/**
 * Reads a raster band as a double array.
 * Returns 0 in case of success, a non-zero value in case of an error.
 */
int raster_band_read_double (char *raster, int band, double **data, int *rasterX, int *rasterY);


/**
 * Reads a raster band as a long array.
 * Returns 0 in case of success, a non-zero value in case of an error.
 */
int raster_band_read_long (char *raster, int band, long **data, int *rasterX, int *rasterY);



/**
 * Writes a double data array as a raster band into an output file.
 * Return 0 in case of success, and a non-zero value in case of an error.
 */
int raster_band_write_double(char *raster, char *format, int band, double *adfGeoTransform, 
               double *data, int rasterX, int rasterY);




/**
 * Converts pixel coordinates to geographic coordinates using the values of the affine transform.
 * The affine transform needed for this function is returned by the GDALGetGeoTransform() function.
 * Note that for an north-up image, the upper left raster corner is at pixel coordinates 0.0/0.0, 
 * and the center of the first pixel is at coordinates 0.5/0.5.
 * @param  padfTransform  the affine transformation coefficients as returns by GDALGetGeoTransform()
 * @param  pixelX      the pixel coordinate in x
 * @param  pixelY      the pixel coordinate in y
 * @param  geoX      a pointer to the geographic coordinate in x
 * @param  geoY      a pointer to the geographic coordinate in y
 */
void pixel_coord_to_geo(double *padfTransform, double pixelX, double pixelY, double *geoX, double *geoY);



/**
 * Converts geographic coordinates into pixel coordinates using the values of the affine transform.
 * The affine transform needed for this function is returned by the GDALGetGeoTransform() function.
 * This function is the inverse of the pixel_coord_to_geo() function. See the documentation of this
 * function for more details.
 */
void geo_coord_to_pixel(double *padfTransform, double geoX, double geoY, double *pixelX, double *pixelY);





/**
 * Returns the bounding box in geographic coordinates of the provided raster dataset.
 * The bounding box is written to the bbox array with the following values:
 *  bbox[0] = minX
 *  bbox[1] = minY
 *  bbox[2] = maxX
 *  bbox[3] = maxY
 */
void raster_bbox (GDALDatasetH raster, double *bbox);





