% 2/10/21
%Example file for calling grain boundary detection algorithm using n-channel images

im1_tif = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif", 'r');
im2_tif = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_98.tif", 'r');
im3_tif = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_99.tif", 'r');
im4_tif = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_100.tif", 'r');

a = double(read(im1_tif));
b = double(read(im2_tif));
c = double(read(im3_tif));
d = double(read(im4_tif));

% a = double(imread('Pics/bundle1/a.png'));
% b = double(imread('Pics/bundle1/b.png'));
% c = double(imread('Pics/bundle1/c.png'));
% d = double(imread('Pics/bundle1/d.png'));
displayImages = false;
fuseImages = true;

%Scale TIF image to be within the range 0-255
m = max(max(max(cat(3, a, b, c, d))));
a = a * (255/m);
b = b * (255/m);
c = c * (255/m);
d = d * (255/m);
%imageData = imageData*(255/m);
%Call algorithm to retrieve skeleton
tic
skeleton = mainMatlab(a, b, c, d);
toc

imwrite(skeleton, 'Pics/skel_out.png');

if (fuseImages)
    J = imresize(skeleton, size(a, 1)/size(skeleton,1));
    C = imfuse(a, uint8(J));
    imwrite(rgb2gray(C), 'Pics/skel_out_overlay.png');
end

if (displayImages && fuseImages)
    %colormap gray;

    %Display the original grain image
    figure;
    imagesc(imageData);
    colormap gray;
    colorbar;
    title('Grain Image')

    %Display the outputted skeleton
    figure;
    imagesc(skeleton);
    colormap gray;
    colorbar;
    title('Skeleton')

    %Resize the skeleton to the size of the TIF image
    figure;
    J = imresize(skeleton,2048/100);
    imagesc(J)
    colormap gray;
    colorbar;
    title('ResizedSkeleton')

    %Overlay the skeleton onto the TIF image
    figure;
    C = imfuse(imageData,J);
    imagesc(rgb2gray(C));
    colormap gray;
    colorbar;
    title('Grains and Resized Skeleton')
end

clear all;

