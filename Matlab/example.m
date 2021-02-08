
%8/4/20

%Example file for calling grain boundary detection algorithm

%{
To run the mainMatlab function, you must compile it first. To do so, run
                mex main Matlab.cpp.
If an error occurs you probably need to download a Mex
Compiler. The easiest way to do this is to go into Add-Ons 
and download MinGW.
%}


%Testing algorithm on 'CCDPRE~1.TIF' image
%Get TIF image
%t = Tiff('Pics/Pt_94kx_Conical 5sec_1fs_70umObj_5frames_02_1.ser_96.png','r');
t = Tiff('Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif','r');
displayImages = false;
fuseImages = false;
% Read in the image and convert it to double
imageData = double(read(t));

% image = imread('Pics/Pt_94kx_Conical 5sec_1fs_70umObj_5frames_02_1.ser_100.png');
% if size(image, 3) > 1
%     imageData = double(rgb2gray(image));
% else
%     imageData = double(image);
% end

% %Display the original grain image
% figure;
% imagesc(imageData);
% colormap gray;
% colorbar;
% title('Grain Image')

%Scale TIF image to be within the range 0-255
m = max(max(imageData));
imageData = imageData*(255/m);
%Call algorithm to retrieve skeleton
tic
skeleton = mainMatlab(imageData);
toc

imwrite(skeleton, 'Pics/skel_out.png');

if (fuseImages)
    J = imresize(skeleton, size(imageData, 1)/size(skeleton,1));
    C = imfuse(imageData, uint8(J));
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
