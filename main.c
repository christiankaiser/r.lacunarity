/**
 r.lacunarity
 Computes the lacunarity of a raster file
 Author:  Christian Kaiser, christian.kaiser@unil.ch
*/


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "lacunarity.h"
#include "gdal.h"


static char *usage[] = {
"\nr.lacunarity -- computes the lacunarity for a GDAL-compatible raster file.\n\n",
"SYNOPSIS\n",
"   r.lacunarity \n",
"      [--help]\n",
"      [--spatial] [--3d]\n",
"      --input input_raster [--band input_band] [--binary]\n",
"      [--binaryThreshold 1] [--mwin 5]\n",
"      [--gbox 3] [--gboxMin 3] [--gboxMax 30] [--gboxStep 1]\n",
"      [--output output_raster_path] [--format format]\n\n",
"DESCRIPTION\n",
"   The following options are available:\n\n",
"   -h\n",
"   --help\n",
"      Shows this usage note.\n\n",
"   -s\n",
"   --spatial\n",
"      Produces a spatial image of lacunarity using a moving window technique.\n",
"      If this flag is selected, an output image raster (--output) must be \n",
"      provided. This flag is not compatible with the gboxMin, gboxMax and\n",
"      gboxStep options. You may want to provide the moving window size using\n",
"      the mwin option, and an output raster format using the format option.\n\n",
"   --3d\n",
"      For non binary images, considers the gliding box in three dimensions.\n",
"      This flag is therefore not compatible with the binary flag. Preference\n",
"      will be given to the binary flag.\n",
"      If selected, the image is considered as being 3D instead of the layered\n",
"      analysis approach.\n\n",
"   -i input_raster\n",
"   --input input_raster\n",
"      The raster for which we should compute the lacunarity.\n\n",
"   -b band\n",
"   --band band\n",
"      The raster band for which we should compute the lacunarity. Default is 1.\n\n",
"   --binary\n",
"      The input raster should be treated as a binary image instead of\n",
"      grayscale.\n\n",
"   --binaryThreshold\n",
"      The pixel value used as a threshold for binary images. Pixels smaller\n",
"      than this value are converted to 0, pixels greater or equal than the\n",
"      the threshold are converted to 1. This option is ignored if the binary\n",
"      flag is not set. Default value is 1.\n\n",
"   -m moving_window_size\n",
"   --mwin moving_window_size\n",
"      The size of the moving window. If you specify a value for this, you must\n",
"      also choose the spatial flag, otherwise this value will be ignored.\n",
"      The default value for this option is 5.\n\n",
"   -g gliding_box_size\n",
"   --gbox gliding_box_size\n",
"      The size of the gliding box used for estimate the lacunarity.\n\n",
"   --gboxMin gliding_box_minimum_size\n",
"      The minimum gliding box size if you want to compute the lacunarity for\n",
"      more than one gliding box size. This option is not compatible with the\n",
"      gbox option and with the spatial option as it is not possible to compute\n",
"      the spatial lacunarity for several gliding boxes.\n\n",
"   --gboxMax gliding_box_maximum_size\n",
"      The maximum gliding box size if you want to compute the lacunarity for\n",
"      more than one gliding box size. This option is not compatible with the\n",
"      gbox option and with the spatial option.\n\n",
"   --gboxStep gliding_box_step_size\n",
"      If you give a value for the gboxMin and gboxMax options, you can specify\n",
"      a step size for the gliding box size. Default is 1.\n\n",
"   -o output_raster_path\n",
"   --output output_raster_path\n",
"      The path to the output raster file. You need to select the spatial flag\n",
"      in order to make something useful.\n\n",
"   -f format\n",
"      Format for the output raster file. Default is HFA.\n",
"      The following formats are supported:\n",
"      GTiff    : GeoTIFF\n",
"      HFA      : Erdas Imagine (.img)\n",
"      AAIGrid  : Arc/Info ASCII Grid\n",
"      PNG      : Portable Network Graphics\n",
"      JPEG     : JPEG image\n",
"      GIF      : Graphics Interchange Format\n",
"      PCIDSK   : PCIDSK Database File\n",
"      PCRaster : PCRaster Raster File\n",
"      GMT      : GMT NetCDF Grid Format\n",
"      JPEG2000 : JPEG-2000\n",
"      RST      : Idrisi Raster A.1\n",
"      ENVI     : ENVI .hdr Labelled\n\n",
"REFERENCES\n",
"   Mandelbrot, B. (1983). The fractal geometry of nature. New York: Freeman.\n",
"   Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of random\n",
"      and deterministic fractal sets. Physical Review A, 44(6), 3552-3558.\n",
"   Plotnick, R., Gardner, R., Hargrove, W., Prestegaard, K. and Perlmutter, M.\n",
"      Lacunarity analysis: a general technique for the analysis of spatial\n",
"      patterns. Physical Review E, 53(5), 5461-5468.\n",
"   Myint, S. and Lam, N. (2005). A study of lacunarity-based texture analysis\n",
"      approaches to improve urban image classification. Computers, Environment\n",
"      and Urban Systems, 29(5), 501-523.\n",
"\n",
"BUGS\n",
"   Please send any comments or bug reports to christian@361degres.ch.\n\n",
"VERSION\n",
"   1.0.1 (15.7.2009)\n\n",
"AUTHOR\n"
"   Christian Kaiser, University of Lausanne <christian@361degres.ch>\n\n",
"ACKNOWLEDGEMENTS\n",
"   This work has been supported by the Swiss National Foundation, through the projects \"Urbanization\n",
"   Regime and Environmental Impact: Analysis and Modelling of Urban Patters, Clustering and\n",
"   Methamorphoses\" (no 100012-113506, see www.clusterville.org for more information) and \n",
"   \"GeoKernels\": Kernel-Based Methods for Geo- and Environmental Sciences, phase 2 \n",
"   (no 200020-121835, see www.geokernels.org for more information).\n\n",
NULL};





int main (int argc, const char *argv[]){
  int c;
  int index;
  
  int spatial;              // Should we compute the spatial lacunarity.
  char *input_raster;        // Path to the input raster file.
  int band;            // Input raster band.
  int binary;            // Is input raster band binary?
  long binaryThreshold;      // The binary threshold.
  int f3d;            // 3D flag.
  int mwin;            // The size of the moving window for spatial lacunarity.
  int gbox;            // The size of the gliding box.
  int gbox_use_min_max;      // Should we use min/max values for gliding box size?
  int gbox_min;          // The minimum size of the gliding box.
  int gbox_max;          // The maximum size of the gliding box.
  int gbox_step;          // The gliding box step size.
  char *output_file;        // Path to the output image file.
  char *format;          // Output image file format.
  char defaultFormat[] = "HFA";  // Default output image file format.
  
  int ok;
  
  extern int optind;
  extern int optopt;
  extern char * optarg;
  
  // Initialize variables with default values.
  spatial = 0;
  input_raster = NULL;
  band = 1;
  binary = 0;
  binaryThreshold = 1;
  f3d = 0;
  mwin = 5;
  gbox = 3;
  gbox_use_min_max = 0;
  gbox_min = 3;
  gbox_max = 30;
  gbox_step = 1;
  output_file = NULL;
  format = defaultFormat;
  
  // Process command line
  while (1){
    static struct option long_options[] =
    {
      {"help",              no_argument,        0,  'h'},
      {"spatial",           no_argument,        0,  's'},
      {"input",             required_argument,  0,  'i'},
      {"band",              required_argument,  0,  'b'},
      {"binary",            no_argument,        0,  'n'},
      {"binaryThreshold",   required_argument,  0,  'd'},
      {"3d",                no_argument,        0,  '3'},
      {"mwin",              required_argument,  0,  'm'},
      {"gbox",              required_argument,  0,  'g'},
      {"gboxMin",           required_argument,  0,  'p'},
      {"gboxMax",           required_argument,  0,  'q'},
      {"gboxStep",          required_argument,  0,  't'},
      {"output",            required_argument,  0,  'o'},
      {"format",            required_argument,  0,  'f'},
      {0, 0, 0, 0}
    };
    
    c = getopt_long(argc, (char**)argv, "hsi:b:nd:3m:g:p:q:t:o:f:", long_options, NULL);
    
    // Detect the end of the options.
    if (c == -1) break;
    
    switch (c) {
      case 'h':
        index = 0;
        while (usage[index] != NULL) {
          printf("%s", usage[index]);
          index++;
        }
        return 0;
      
      case 's':
        spatial = 1;
        break;
        
      case 'i':
        input_raster = optarg;
        break;
        
      case 'b':
        band = atoi(optarg);
        break;
      
      case 'n':
        binary = 1;
        break;
      
      case 'd':
        binaryThreshold = atol(optarg);
        break;
      
      case '3':
        f3d = 1;
        break;
        
      case 'm':
        mwin = atoi(optarg);
        break;
      
      case 'g':
        gbox = atoi(optarg);
        break;
      
      case 'p':
        gbox_min = atoi(optarg);
        gbox_use_min_max = 1;
        break;
      
      case 'q':
        gbox_max = atoi(optarg);
        gbox_use_min_max = 1;
        break;
      
      case 't':
        gbox_step = atoi(optarg);
        break;
        
      case 'o':
        output_file = optarg;
        break;
        
      case 'f':
        format = optarg;
        break;
        
      case '?':
        return 1;
        
      default:
        abort();
    }
  }
  
  
  if (input_raster == NULL){
    fprintf(stderr, "Error. You must provide at least an input raster file.\n");
    fprintf(stderr, "Use r.lacunarity -h to get help on the input parameters.\n");
    return 1;
  }
  
  GDALAllRegister();
  
  if (spatial == 1){
    ok = spatial_lacunarity(input_raster, band, binary, binaryThreshold, f3d, gbox, mwin, output_file, format);
  }else{
    if (gbox_use_min_max == 0){
      gbox_min = gbox;
      gbox_max = gbox;
      gbox_step = 1;
    }
    ok = lacunarity(input_raster, band, binary, binaryThreshold, f3d, gbox_min, gbox_max, gbox_step);
  }
  
  fprintf(stdout, "r.lacunarity done.\n");
  return ok;
}







