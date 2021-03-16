
%tifdata = Tiff("Pics/Pt_94kx_Conical 5sec_1fs_20umObj_5frames_02_1.ser_97.tif", "r");
tifdata = Tiff("../../../C++/Pics/fullsize_tensor/tensor3.tif", "r");
im = double(read(tifdata));
im_range = max(max(im)) - min(min(im));
im = (im - min(min(im))) * (255/im_range);

kernal_gaus = make_gaus(15, 4);
kernal_lapofgaus = make_lapofgaus(9, 1.4);
% print_kernal(kernal_lapofgaus);
print_kernal_cpp(kernal_lapofgaus);

% im_conv = conv2(im, kernal_lapofgaus);
% % im_filt = imdiffusefilt(im);
% 
% im_conv_range = max(max(im_conv)) - min(min(im_conv));
% im_conv = (im_conv - min(min(im_conv))) * ((2^16-1)/im_conv_range);
% 
% imshow(uint16(im_conv));
% %imwrite(uint16(im_conv), "gaussianblur_sigma-" + sigma + ".png")
% imwrite(uint16(im_conv), "LoG-" + 1.4 + ".tif")

function kernal = make_lapofgaus(kernal_width, sigma)
%  kernal_width = 15;
%  sigma = 4;
%  K = 1;

  kernal_center = [ceil(kernal_width/2), ceil(kernal_width/2)];
  kernal = zeros(kernal_width);
  %rs = zeros(kernal_width);

  for i = 1:kernal_width
    for j = 1:kernal_width
%      r = sqrt((i - kernal_center(1))^2 + (j - kernal_center(2))^2);
      r2 = abs((i - kernal_center(1))^2 + (j - kernal_center(2))^2);
  %    rs(i, j) = r;
      K = - 1/(pi * sigma^4);
      kernal(i, j) = K * (1 - (r2)/(2 * sigma^2)) * exp(-(r2)/(2*sigma^2));
    end
  end
end

function kernal = make_gaus(kernal_width, sigma)
%  kernal_width = 15;
%  sigma = 4;
%  K = 1;

  kernal_center = [ceil(kernal_width/2), ceil(kernal_width/2)];
  kernal = zeros(kernal_width);
  %rs = zeros(kernal_width);

  for i = 1:kernal_width
    for j = 1:kernal_width
      r = sqrt((i - kernal_center(1))^2 + (j - kernal_center(2))^2);
  %    rs(i, j) = r;
      K = 1/(2 * pi * sigma^2);
      kernal(i, j) = K * exp(-(r^2)/(2*sigma^2));
    end
  end
end

function print_kernal(kernal)
  kernal_width = size(kernal, 1);

  for i = 1:kernal_width
    for j = 1:kernal_width
      fprintf("%1.2f ", kernal(i,j));
    end
    fprintf("\n");
  end
end

function print_kernal_cpp(kernal)
  kernal_width = size(kernal, 1);

  fprintf("{");
  for i = 1:kernal_width
    for j = 1:kernal_width
      fprintf("%1.4f", kernal(i,j));

      if (i == kernal_width)
        if (j == kernal_width)
          continue;
        end
      end
      fprintf(", ");
    end
    if (i ~= kernal_width)
      fprintf("\n");
    end
  end
  fprintf("}\n");
end
