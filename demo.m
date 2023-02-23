clear all;
%% Load HSI data
addpath Dataset\
load wdc_demo.mat; % Washionton DC dataset in case3
% Noisy_Img = Noisy_Img(150:200,150:200,:);
% Img = Img(150:200,150:200,:);% for quick display
[m,n,p]  = size(Noisy_Img);
%% Call the main function
blocksize = 20;  
stepsize  = 4;
a = 20;
mu = 0.008; %
lambda = 0.006; % mu and lambda should be changed in different datasets
fprintf("=== a = %d ; mu = %f ; lambda = %f ===\n",a,mu,lambda);
[ output_image ] = NFF_HSIdenoise(Noisy_Img,blocksize,stepsize,a,mu,lambda);
% save nff_dc3 Img Noisy_Img output_image;
%% quality assess
addpath(genpath('quality_assess'));
npsnr_v = zeros(1,p);nssim_v = zeros(1,p);nfism_v = zeros(1,p);nsam_v = zeros(1,p);nergas_v = zeros(1,p);
mpsnr_v = zeros(1,p);mssim_v = zeros(1,p);mfism_v = zeros(1,p);msam_v = zeros(1,p);mergas_v = zeros(1,p);
for i=1:1:p
    J=255*Img(:,:,i);
    K=255*Noisy_Img(:,:,i);
    I=255*output_image(:,:,i);
    [npsnr_v(i),nssim_v(i),nfsim_v(i),nergas_v(i),nmsam_v(i)] = MSIQA(J,K);
    [mpsnr_v(i),mssim_v(i),mfsim_v(i),mergas_v(i),mmsam_v(i)] = MSIQA(J,I);
end
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v));





