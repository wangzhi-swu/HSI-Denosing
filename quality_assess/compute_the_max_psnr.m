clear;
close all;
% load('result_31-49_abs.mat');
load('result_220k-300k.mat')
% load('result_305k-360k.mat')
% load('result_405k-500k.mat')
psnr = zeros(1,size(result,2)/9);
for i = 1:size(result,2)/9
    for j = 1:9
        k = (i-1)*9+j;
        psnr(i)=psnr(i) + result(k).result(1);
    end
%         disp(k)
end
[max_psnr, index] = max(psnr)

mean_psnr = mean(psnr)
