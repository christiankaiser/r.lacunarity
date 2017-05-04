IDIR = /Library/Frameworks/GDAL.framework/unix/include
LDIR = /Library/Frameworks/GDAL.framework/unix/lib
CFLAGS = -O3 -I$(IDIR)
CC = gcc
LIBOPTS =
LIBS = -L$(LDIR) -lgdal -lm

default: all


r_lacunarity:main.o lacunarity.o raster.o
	$(CC) $(CFLAGS) $(LIBOPTS) -o r.lacunarity main.o lacunarity.o raster.o $(LIBS)

lacunarity.o:lacunarity.c Makefile
	$(CC) $(CFLAGS) -c lacunarity.c

raster.o:raster.c Makefile
	$(CC) $(CFLAGS) -c raster.c

main.o:main.c Makefile
	$(CC) $(CFLAGS) -c main.c

all: r_lacunarity

clean:
	rm main.o lacunarity.o raster.o r.lacunarity