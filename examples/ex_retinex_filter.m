% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

clearvars all;
close all;

initialData = imread('yourpicture.jpg');

[h,s,v] = rgb2hsv(initialData); % Elad used HSV and applied retinex to V

originalGrayImg = double(rgb2gray(initialData));
originalRGBImage = double(initialData);

subplot(221);
imagesc(initialData);
title('original image');
%colormap gray;

% parameters
sigmaSpatial = 30;
sigmaRange = 0.2;

samplingSpatial = sigmaSpatial;
samplingRange = sigmaRange;

gamma = 1.5;

subplot(224);
v = retinexFilter(v, sigmaSpatial, sigmaRange, samplingSpatial, samplingRange, gamma, 0);

imhsv(:,:,1) = h;
imhsv(:,:,2) = s;
imhsv(:,:,3) = v;
outRGBImg = hsv2rgb(imhsv);

imagesc(outRGBImg);
title('after retinex');