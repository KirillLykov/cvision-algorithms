% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

clearvars all;
close all;

image = double(rgb2gray(imread('sheef.jpg')));
depth = shapeFromShading(image, 1000, 1/32);

figure;
mesh(depth);
title('Surface reconstruction');
