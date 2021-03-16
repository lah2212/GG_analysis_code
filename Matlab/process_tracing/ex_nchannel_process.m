% 2/10/21

%Example file for calling grain boundary detection algorithm using n-channel images
function ex_nchannel_process(varargin)
    folder = "bundle1";
    filetype = "png";
    scale = 0;
    if (nargin == 1)
      folder = varargin{1};
    elseif (nargin == 2)
      folder = varargin{1};
      filetype = varargin{2};
    elseif (nargin >= 2)
      folder = varargin{1};
      filetype = varargin{2};
      scale = double(varargin{3});
    end

    displayImages = false;
    fuseImages = false;

    cd("Pics/" + folder);
    folderContents = dir('*.' + filetype);

    images = [];
    for i = 1:length(folderContents)
      if (filetype == "tif")
        a_tiff = Tiff(folderContents(i).name, 'r');
        a = read(a_tiff);
      else
        a = imread(folderContents(i).name);
      end

      if (size(a,3) == 3)
        a_bw = rgb2gray(a);
        images = cat(3, images, double(a_bw));
      else
        images = cat(3, images, double(a));
      end
    end

    m = max(max(max(images)));
    images = images * (255/m);
    cd('../..');
    args = num2cell(images, [1 2]);

    tic
    [edges, threshold, denoise] = mainMatlabProcess(args{:}, [scale]);
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
