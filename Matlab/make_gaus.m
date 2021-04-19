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
      %K = 1
      kernal(i, j) = K * exp(-(r^2)/(2*sigma^2));
    end
  end
end
