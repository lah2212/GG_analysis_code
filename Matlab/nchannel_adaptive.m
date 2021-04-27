% 2/10/21
%Example file for calling grain boundary detection algorithm using n-channel images
function nchannel_adaptive(varargin)
  if (nargin > 1)
    dir_name = "Pics/" + varargin{1} + "/";
    k = varargin{2};
  elseif (nargin > 0)
    dir_name = "Pics/" + varargin{1} + "/";
  else
    dir_name = "Pics/results/";
    k = 0.0023;
  end

  if (~isdir(dir_name))
    mkdir(dir_name);
  end

  tic
  fnames = [ "Pics/Pt170_STEM_110kX_C2(100)_CL205mm_03_noscale.tif"];
%  fnames = [ "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_99.tif",
%        "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_96.tif",
%        "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif",
%        "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_98.tif",
%        "Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_100.tif" ];
  NUM_IMAGES = length(fnames);
  %SCALE_IMG = 1;
  SCALE_IMG = 700/2048;
  imgs = [];

  for i = 1:NUM_IMAGES
    tif = Tiff(fnames(i), 'r');
    img = double(read(tif));
    range = max(max(img)) - min(min(img)); img_n = (255 / range) * (img - min(min(img)));
  %  imgs = cat(3, imgs, img_n);
    imgs = cat(3, imgs, imresize(img_n, SCALE_IMG));
  end

  displayImages = false;
  fuseImages = true;
  pgen_fname = "pgen-700x700";
  pgen_on = false;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Calculate q
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [M, N] = size(imgs(:,:,1));
  gauss = size(M,N);
%  scale = 0.09;
  scale = 0.05;
%  k = .0075;          %parameter edge detector
%  k = .00023;          %parameter edge detector

  % This is just so I can test the edge detection AFTER the smoothing and first order detector
  if (~pgen_on)

    %Calculate Gaussian
    for i=1:M
        for j=1:N
            gauss(i,j) = sqrt(2*pi)*scale*exp(-2*pi^2*scale^2*(i^2+j^2));
        end
    end
    GAUSS = fftshift(fft2(gauss));

    %%perform convolution

    grad_sum = zeros(M+2, N+2);
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
          for l=1:N
              conv(j,l) = conv(j,l)*(-1)^(j+l);
          end
      end

%       lapofgaus = make_lapofgaus(1.4);
%       magq1 = conv2(conv, lapofgaus);
      magq1 = find_grad(conv);

      % grad_sum(:,:,1) = grad_sum(:,:,1) + gradq1(:,:,1);
      % grad_sum(:,:,2) = grad_sum(:,:,2) + gradq1(:,:,2);
      
      % NORMAL SUM %
      grad_sum = grad_sum + magq1;
      % MAX VALUE %
%       grad_sum(grad_sum < magq1) = magq1(grad_sum < magq1);

      % MIN VALUE %
%      if (i == 1)
%        grad_sum = magq1;
%      else
%        grad_sum(grad_sum > magq1) = magq1(grad_sum > magq1);
%      end

  %    q = q + (1+1./(1+k*magq1.^2));
      imwrite_norm(grad_sum, dir_name + "grad_sum_" + num2str(i));
    end
    % q = q + (1 + 1./(1 + k * (grad_sum(:,:,1).^2 + grad_sum(:,:,2).^2)));
    grad_sum = grad_sum ./ NUM_IMAGES;
    q = 1 + 1./(1 + k * grad_sum.^2);
    %q = q/NUM_IMAGES;

%     q_rs = zeros(M + 2, N + 2);
%     b = (size(q, 1) - (M + 2))/2;
%     q_rs = q(b + 1:M + b + 2, b + 1:N + b + 2);

    imwrite_norm(grad_sum, dir_name + "grad_sum_total");
    imwrite_norm(q, dir_name + "q");

    imgs_s = [];
    for i = 1:NUM_IMAGES
      fprintf("Smoothing image %d\n", i);
      [u, q] = AdaptiveSmoothingUpwind(imgs(:,:,i), 35, 20, q);
      imgs_s = cat(3, imgs_s , u);
      imwrite_norm(u, dir_name + "smoothed_img_" + num2str(i));
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
    pgen = load(pgen_fname);
    tensor = pgen.tensor;
  end

  imwrite(uint8(tensor * 255/max(max(tensor))), dir_name + "1_first_order.png");
  tensor = exp(-0.1667 * sqrt(tensor));
  imwrite(uint8(tensor * 255/max(max(tensor))), dir_name + "2_tensor.png");
  % tensor = conv2(tensor, avgfilter);
  % imwrite(uint8(tensor * 255/max(max(tensor))), "Pics/3_avgfilter.png");
  % tensor = conv2(tensor, lapfilter);
  % imwrite(uint8(tensor * 255/max(max(tensor))), "Pics/3_lapfilter.png");
  lapofgaus = make_lapofgaus(1.4);
  tensor = conv2(tensor, lapofgaus);
  imwrite_norm(tensor, dir_name + "3_lapfilter.png");
  tensor(tensor < 0) = 0;
  tensor = exp(-10 .* tensor);
  imwrite_norm(tensor, dir_name +"4_edges.png");

  border = (size(tensor, 1) - M) / 2;
  % edges = 1 - tensor(1 + border:M + border, 1 + border:N + border);
  edges = tensor(2 + border:M + border - 1, 2 + border:N + border - 1);
  %edges = (edges - min(min(edges))) * 1/(max(max(edges)) - min(min(edges)));

  toc

  % Call algorithm to retrieve skeleton
  tic
  skeleton = edgeProcessing(edges);
  toc

  imwrite(skeleton, dir_name + 'skel_out.png');

  if (fuseImages & ~pgen_on)
      J = imresize(skeleton, size(imgs(:,:,1), 1)/size(skeleton,1));
      C = imfuse(imgs(:,:,1), uint8(J));
      imwrite(rgb2gray(C), dir_name + 'skel_out_overlay.png');
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

function imwrite_norm(img, fname)
  img_norm = (img - min(min(img))) * (255/(max(max(img)) - min(min(img))));
  if (contains(fname, [".p", ".j", ".t"]))
    imwrite(uint8(img_norm), fname);
  else
    imwrite(uint8(img_norm), fname + ".png");
  end
end

function gradmag = find_grad(img)
  [M, N] = size(img);
  % make image with one pixel border
  img_b = zeros(M+2,N+2);
  img_b(2:M+1,2:N+1) = img;
  grad = zeros(M+2,N+2,2);

  % gradq1(:,:,1) is the y (row) grad and gradq1(:,:,2) is the x (column) grad 
  % and magq1 is the magnitude of those gradients
  grad(2:M+1,2:N+1,1) = img_b(3:M+2,2:N+1) - img_b(1:M,2:N+1);
  grad(2:M+1,2:N+1,2) = img_b(2:M+1,3:N+2)- img_b(2:M+1,1:N);
  gradmag(:,:) = sqrt(grad(:,:,1).^2 + grad(:,:,2).^2);
end
