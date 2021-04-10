% 2/10/21
%Example file for calling grain boundary detection algorithm using n-channel images
clear;

fnames = [ "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_96.tif",
      "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif",
      "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_98.tif",
      "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_99.tif",
      "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_100.tif" ];
pregenfname = "Pics/pre-gen-tensor.tif";
NUM_IMAGES = length(fnames);
imgs = [];

for i = 1:NUM_IMAGES
  tif = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif", 'r');
  img = double(read(tif));
  range = max(max(img)) - min(min(img));
  img_n = (255 / range) * (img - min(min(img)));
%  imgs = cat(3, imgs, img_n);
  imgs = cat(3, imgs, imresize(img_n, [700, 700]));
end

displayImages = false;
fuseImages = true;
pregenTensor = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M, N] = size(imgs(:,:,1));
gauss = size(M,N);
scale = .05;
k = .0075;          %parameter edge detector

if (~pregenTensor)

  %Calculate Gaussian
  for i=1:M
      for j=1:N
          gauss(i,j) = sqrt(2*pi)*scale*exp(-2*pi^2*scale^2*(i^2+j^2));
      end
  end
  GAUSS = fftshift(fft2(gauss));

  %%perform convolution

  q = zeros(M+2,N+2);
  for i = 1:NUM_IMAGES
    %calculate and shift FFT
    F = fftshift(fft2(imgs(:,:,i)));

    %performs multiplication
    R = F.*GAUSS;

    %calculate and shift inverse of result
    conv = ifft2(R);
    conv = real(conv);

    for j=1:M
        for k=1:N
            conv(j,k) = conv(j,k)*(-1)^(j+k);
        end
    end

    %calculate q and q*
    q1 = zeros(M+2,N+2);
    q1(2:M+1,2:N+1) = conv;
    gradq1 = zeros(M+2,N+2,2);

    % gradq1(:,:,1) is the y (row) grad and gradq1(:,:,2) is the x (column) grad 
    % and magq1 is the magnitude of those gradients
    gradq1(2:M+1,2:N+1,1) = q1(3:M+2,2:N+1)-q1(1:M,2:N+1);
    gradq1(2:M+1,2:N+1,2) = q1(2:M+1,3:N+2)-q1(2:M+1,1:N);
    magq1(:,:) = sqrt(gradq1(:,:,1).^2+gradq1(:,:,2).^2);

    q = q + (1+1./(1+k*magq1.^2));
  end
  q = q/NUM_IMAGES;
  imgs_s = [];
  for i = 1:NUM_IMAGES
    fprintf("Smoothing image %d\n", i);
    [u, q] = AdaptiveSmoothingUpwind(imgs(:,:,i), 35, 20, q);
    imgs_s = cat(3, imgs_s , u);
  end

  avgfilter = [[1, 1, 1],
              [1, 1, 1],
              [1, 1, 1]] * (1/9);
  lapfilter = [[1, 4, 1],
              [4, -20, 4],
              [1, 4, 1]] * (1/6);
  tensor = zeros(M + 2, N + 2);

  for i = 1:NUM_IMAGES
    img_b = zeros(M+2, N+2);
    img_b(2:M+1, 2:N+1) = imgs_s(:,:,i);
    gradq1 = zeros(M+2,N+2,2);

    % gradq1(:,:,1) is the y (row) grad and gradq1(:,:,2) is the x (column) grad 
    % and tensor is the magnitude of those gradients
    gradq1(2:M+1,2:N+1,1) = img_b(3:M+2, 2:N+1) - img_b(1:M, 2:N+1);
    gradq1(2:M+1,2:N+1,2) = img_b(2:M+1, 3:N+2) - img_b(2:M+1,1:N);
    tensor(:,:) = tensor + (sqrt(gradq1(:,:,1).^2+gradq1(:,:,2).^2)).^2;
  end

else
  pregen_tif = Tiff(pregenfname);
  tensor = double(read(pregen_tif));
end

imwrite(tensor, "Pics/1_first_order.tif");
tensor = exp(-0.1667 * sqrt(tensor));
imwrite(tensor, "Pics/2_tensor.tif");
% tensor = conv2(tensor, avgfilter);
% tensor = conv2(tensor, lapfilter);
imwrite(tensor, "Pics/3_lapfilter.tif");
tensor(tensor < 0) = 0;
tensor = exp(-10 * tensor);
imwrite(tensor, "Pics/4_edges.tif");

edges = tensor(4:M+1, 4:N+1);

% Call algorithm to retrieve skeleton
tic
%args = num2cell(imgs, [1 2]);
%skeleton = mainMatlab(args{:});
skeleton = edgeProcessing(edges);
toc

imwrite(skeleton, 'Pics/skel_out.png');

if (fuseImages)
    J = imresize(skeleton, size(imgs_s(:,:,1), 1)/size(skeleton,1));
    C = imfuse(imgs_s(:,:,1), uint8(J));
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

