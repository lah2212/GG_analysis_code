% 2/10/21
% This program is designed to be rudimentary alignment checking software

im1_tif = Tiff("../Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif", 'r');
im2_tif = Tiff("../Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_98.tif", 'r');
im3_tif = Tiff("../Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_99.tif", 'r');

im1 = double(read(im1_tif));
im2 = double(read(im2_tif));
im3 = double(read(im3_tif));

im = cat(3, im1, im2, im3);
m = max(max(max(im)));
im = im * 255/m;

imwrite(uint8(im), "Pics/image_alignment_test.png");
