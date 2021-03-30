## main.cpp
main.cpp is devoted to doing what example.m did in the Matlab code, image IO. Currently it only supports tifs.

Resources on importing tifs:
https://research.cs.wisc.edu/graphics/Courses/638-f1999/libtiff_tutorial.htm
https://www.cs.rochester.edu/users/faculty/nelson/courses/vision/resources/tiff/libtiff.html#TIFFRGBAImage
http://maptools-org.996276.n3.nabble.com/TIFFWriteScanline-works-on-Windows-and-Linux-but-fails-on-Mac-td969.html

## main_nch.cpp
This just feeds an n-channel image stack to the algorithm, and uses QTiffIO to import and export images.

## qtiffio.cpp
This makes reading and writing tiffs faster and cleaner. It uses the same techniques used in main.cpp and is **not** generalized, so it is only designed to work with the images that we're currently working with (March 2021, B&W 16bit).

The other files are used in both the MATLAB and C++ code, so I'll make a README somewhere that describes their purpose.
