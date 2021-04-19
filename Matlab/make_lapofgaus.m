function kernal = make_lapofgaus(sigma, kernal_width)
  if (nargin < 2)
    kernal_width = ceil(6 * sigma);
    if (mod(kernal_width, 2) == 0)
      kernal_width = kernal_width + 1;
    end
  end
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
      %kernal(i, j) = (r2 - 2 * sigma^2)/(sigma^4) * exp(-(r2)/(2*sigma^2));
      kernal(i, j) = K * (1 - (r2)/(2 * sigma^2)) * exp(-(r2)/(2*sigma^2));
    end
  end
end
