imag = imread('C:\Users\18143\Desktop\lake.bmp');

gray_image1 =  rgb2gray(imag);
gray_image = im2double(gray_image1);

[m,n] = size(gray_image);
s = svd(gray_image);
[U,S,V] = svd(gray_image);
[m1,~] = size(s);

B = zeros(m,n);
singularWeight = 0;
sum_sigular = s'*ones(m1,1);
i = 1;
while singularWeight < 0.9
    B = B + s(i)*U(:,i)*V(:,i)';
    imshow(B)
    singularWeight = singularWeight + s(i)/sum_sigular;
    fprintf("前%d个奇异值所占的比重为:%f\n",i,singularWeight);
    i = i+1;
end
