% This function is supposed to find the translational best match of two images
% by minimizing chi2 values (dist squared). 
% Based on 1998 Carpenter, Rickman, Barmak Section IV.B.5

function [offset_x, offset_y, chi2_min] = chi2goodness(A, B)
  [B_x, B_y] = find(B == 1);
  [A_x, A_y] = find(A == 1);

  %{ ROUGHING PASS %}

  offset_range = [50, -50, 50, -50]; % +x, -x, +y, -y
  offset_n = 10;  % The number of points which I'll be testing on roughing pass

  offset_gap_x = (offset_range(1) - offset_range(2)) / (offset_n - 1);
  offset_gap_y = (offset_range(3) - offset_range(4)) / (offset_n - 1);
  
  offset_map = zeros(offset_n, offset_n, 2);
  for i = 1:offset_n
    offset_map(i, :, 1) = round(offset_range(2) + offset_gap_x * (i - 1));
    offset_map(:, i, 2) = round(offset_range(4) + offset_gap_y * (i - 1));
  end
  
  chi2s = zeros(offset_n);
  for x = 1:offset_n
    for y = 1:offset_n
      offset_x = offset_map(x, y, 1);
      offset_y = offset_map(x, y, 2);
      A_off_x = A_x + offset_x;
      A_off_y = A_y + offset_y;
      
      AId = knnsearch([A_off_x A_off_y], [B_x, B_y]);
      chi2s(x, y) = sum((A_off_x(AId) - B_x) .^ 2 + (A_off_y(AId) - B_y) .^ 2);
    end
  end

  %{ NEIGHBORHOOD TESTS %}

  chi2_min = min(min(chi2s));
  [x, y] = find(chi2_min == chi2s);
  offset_x = offset_map(x, y, 1);
  offset_y = offset_map(x, y, 2);

  while true
    chi2s = zeros(3);
    chi2s(2, 2) = chi2_min;

    for x = -1:1
      for y = -1:1
        if (x == 0 & y == 0)
          continue
        end

        A_off_x = A_x + offset_x + x;
        A_off_y = A_y + offset_y + y;
        
        AId = knnsearch([A_off_x A_off_y], [B_x, B_y]);
        chi2s(x + 2, y + 2) = sum((A_off_x(AId) - B_x) .^ 2 + (A_off_y(AId) - B_y) .^ 2);
      end
    end

    chi2_min = min(min(chi2s));
    [x, y] = find(chi2_min == chi2s);
    if (x == 2 & y == 2) 
      return 
    else
      offset_x = offset_x + (x - 2);
      offset_y = offset_y + (y - 2);
    end
  end
end
