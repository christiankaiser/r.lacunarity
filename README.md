# r.lacunarity

r.lacunarity is a small C program for computing the lacunarity for a GDAL-compatible raster file.

If you don't know what is the lacunarity, have a look into Allain and Cloitre (1991) and/or Plotnick, Gardner, Hargrove, Prestegaard, and Perlmutter (1996); references are below.


## Installation

The program needs the GDAL library installed.

You need to compile the program. There is a Makefile included where you might need to update the path to the GDAL library and headers according to your environment. Once the Makefile is adapted just run `make` and there should be a `r.lacunarity` binary at the end.


## Usage

You can get the information on the input arguments by running `r.lacunarity -h` or `r.lacunarity --help`.

Please note that ther are many bugs in the program and an update would be required. However, the code is provided as is and you can use it if you find it useful.


## Examples

For calculating a global lacunarity with a sliding box size of 5x5:

```bash
r.lacunarity --input testdata/vis_3.tif --gbox 5 --binary
```

## References

- Mandelbrot, B. (1983). The fractal geometry of nature. New York: Freeman.
- Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44(6), 3552-3558.
- Plotnick, R., Gardner, R., Hargrove, W., Prestegaard, K. and Perlmutter, M. (1996). Lacunarity analysis: a general technique for the analysis of spatial patterns. Physical Review E, 53(5), 5461-5468.
- Myint, S. and Lam, N. (2005). A study of lacunarity-based texture analysis approaches to improve urban image classification. Computers, Environment and Urban Systems, 29(5), 501-523.


## Contributions

We are happy to get your comments and contributions. Preferably make a fork of the repo and file your pull requests.


## Author

Christian Kaiser, University of Lausanne <christan.kaiser@unil.ch>


## Acknowledgements

This work has been supported by the Swiss National Foundation, through the projects "Urbanization Regime and Environmental Impact: Analysis and Modelling of Urban Patters, Clustering and Methamorphoses" (no 100012-113506) and  "GeoKernels": Kernel-Based Methods for Geo- and Environmental Sciences, phase 2 (no 200020-121835).