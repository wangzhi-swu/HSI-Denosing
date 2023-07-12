clear all;
%% Load HSI data
% sig = {'1','2','3','4','5'};
% for noiK = 3:3
% noiName = strcat('new_case\case_g', sig{noiK});
% load(noiName, '-mat');
[m,n,p]  = size(Noisy_Img);          
Noisy_Img = myNormalization(Noisy_Img,2);
Img = myNormalization(Img,2);
%% Call the main function
blocksize = 20;  
stepsize  = 4;
a = 2;
lambda = 0.006; % (need be tuned)
fprintf('=== a = %d ;  lambda = %f \n',a,lambda);
[ output_image ] = NFF_HSIdenoise(Noisy_Img,blocksize,stepsize,a,mu,lambda);
% save nff_dc3 Img Noisy_Img output_image;
output_image(output_image>1)=1;
output_image(output_image<0)=0;
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
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(npsnr_v),mean(nssim_v),mean(nfsim_v),mean(nergas_v));
fprintf('psnr = %.4f , ssim = %.4f, fsim = %.4f , ergas = %.4f \n',mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v));
results(:,noiK) = [mean(mpsnr_v),mean(mssim_v),mean(mfsim_v),mean(mergas_v)]';
% end


