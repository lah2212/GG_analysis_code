function print_kernal(kernal)
  kernal_width = size(kernal, 1);

  for i = 1:kernal_width
    for j = 1:kernal_width
      fprintf("%1.2f ", kernal(i,j));
    end
    fprintf("\n");
  end
end
