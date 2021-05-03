% This is a function for converting smoothed images into a pregen tensor

fname = [
  "Pics/results/smoothed_img_1.png",
  "Pics/results/smoothed_img_2.png",
  "Pics/results/smoothed_img_3.png",
  "Pics/results/smoothed_img_4.png",
  "Pics/results/smoothed_img_5.png",
];
NUM_IMAGES = length(fname);
img = imread(fname(1));
[M, N] = size(img);
tensor = zeros(M+2, N+2);

for i = 1:NUM_IMAGES
  img = imread(fname(i));
  %tensor(:,:) = tensor + find_grad(imgs_s(:,:,i)).^2;
  img_b = zeros(M+2, N+2);
  img_b(2:M+1, 2:N+1) = img;
  gradq1 = zeros(M+2,N+2,2);

  % gradq1(:,:,1) is the y (row) grad and gradq1(:,:,2) is the x (column) grad 
  % and tensor is the magnitude of those gradients
  gradq1(2:M+1,2:N+1,1) = img_b(3:M+2, 2:N+1) - img_b(1:M, 2:N+1);
  gradq1(2:M+1,2:N+1,2) = img_b(2:M+1, 3:N+2) - img_b(2:M+1,1:N);
  tensor(:,:) = tensor + (sqrt(gradq1(:,:,1).^2+gradq1(:,:,2).^2)).^2;
end

save("pgen-fullsize");
