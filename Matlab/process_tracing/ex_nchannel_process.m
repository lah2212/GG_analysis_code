% 2/10/21

%Example file for calling grain boundary detection algorithm using n-channel images
function ex_nchannel_process(varargin)
    if (nargin)
      folder = varargin{1};
    else
      folder = "bundle1";
    end

    a = double(imread('Pics/' + folder + '/a.png'));
    b = double(imread('Pics/' + folder + '/b.png'));
    c = double(imread('Pics/' + folder + '/c.png'));
    d = double(imread('Pics/' + folder + '/d.png'));

    displayImages = false;
    fuseImages = false;

    %Scale TIF image to be within the range 0-255
    m = max(max(max(cat(3, a, b, c, d))));
    a = a * (255/m);
    b = b * (255/m);
    c = c * (255/m);
    d = d * (255/m);
    %imageData = imageData*(255/m);
    %Call algorithm to retrieve skeleton
    tic
    [edges, threshold, denoise] = mainMatlabProcess(a, b, c, d);
%    [edges, threshold, denoise] = mainMatlabProcess(a);
    toc

    if ~exist('Pics/' + folder + '/results', 'dir')
      mkdir('Pics/' + folder + '/results');
    end

    imwrite(edges, 'Pics/' + folder + '/results/edges.png');
    imwrite(threshold, 'Pics/' + folder + '/results/threshold.png');
    imwrite(denoise, 'Pics/' + folder + '/results/denoise.png');

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
end
