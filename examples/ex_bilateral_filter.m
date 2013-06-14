% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

clearvars all;
close all;

initialData = imread('cat.png');

originalGrayImg = double(rgb2gray(initialData));

subplot(221);
imagesc(originalGrayImg);
title('original image');
colormap gray;

subplot(222);
% noise functions expects integers
noisyGrayImg = imnoise(rgb2gray(initialData),'gaussian', 0, 0.01);
noisyGrayImg = double(noisyGrayImg);
imagesc(noisyGrayImg);
title('noisy image');

% parameters
sigmaSpatial = 10;
sigmaRange = 100;

samplingSpatial = sigmaSpatial;
samplingRange = sigmaRange;

% brute force Bilateral filter
subplot(223);
tStart = cputime;
filteredImg = regBilateralFilter(noisyGrayImg, 0, sigmaSpatial, sigmaRange, samplingSpatial); 
tEnd = cputime;
imagesc(filteredImg);
title(['brute force BF. Time: ' int2str(tEnd - tStart) 'sec']);

% fast bilateral filter - terms and notations are from the paper
subplot(224);
tStart = cputime;
filteredImg = fastBilateralFilter(noisyGrayImg, sigmaSpatial, sigmaRange, samplingSpatial, samplingRange, 0, 15);
tEnd = cputime;
imagesc(filteredImg);
title(['fast BF. Time: ' int2str(tEnd - tStart) 'sec']);
