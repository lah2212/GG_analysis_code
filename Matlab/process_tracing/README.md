This folder contains programs specifically for generating images of the
de-noising process in the algorithm. Specifically the Non-Isotropic filter, the
double threshold, and the binary de-noising after that.

In Matlab (after compiling the cpp file with mex), run ex_nchannel_process. 

Usage:
```
ex_nchannel_process("<bundle name>", "<filetype>")
```

So to use this program, make a folder named `Pics` in the same directory as 
ex_nchannel_process.m and add folders of images you want to process in `Pics`.
`<bundle name>` should be replaced by the name of the folder the images are in.
`<bundle name>` should also only have images that you want included in the 
n-channel processing. You can have a folder within `<bundle name>` for
anything else. The algorithm will make a folder named `results` within 
`<bundle name>` with all the resulting images. The default `<bundle name>` (if
no arguments are given) is `bundle1`.

The `<filetype>` argument takes a string of the images filetypes. i.e. "tif", 
"jpg", "png", etc. By default (if no second argument is given) the filetype is
set to "png".
