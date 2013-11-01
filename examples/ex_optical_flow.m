% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

clearvars all;
close all;

img1 = imread('/Users/kirill/Desktop/smith1.png');
img1 = double(rgb2gray(img1));

img2 = imread('/Users/kirill/Desktop/smith2.png');
img2 = double(rgb2gray(img2));

images = zeros(size(img2, 1), size(img2, 2), 2);
images(:, :, 1) = img1;
images(:, :, 2) = img2;

[Vx, Vy] = opticalFlow(images, 1, 10);

% Show flow on the picture without showing all vectors (othervise it will
% be all green)

rSize = 5;
for i=1:size(Vx, 1)
    for j=1:size(Vx, 2)
        if (floor(i/rSize) ~= i/rSize || floor(j/rSize) ~= j/rSize)
            Vx(i, j) = 0;
            Vy(i, j) = 0;
        end
    end
end

figure(); 
imshow(uint8(img2)); 
hold on;

% To avoid shoing too many 0 vectos
validIndex = (abs(Vx) ~= 0); 

[gridX, gridY] = meshgrid(1:size(img1, 2), 1:size(img1, 1));
% plot only valid flow vectors
quiver(gridX(validIndex), gridY(validIndex), ...
    Vx(validIndex), Vy(validIndex), 3, 'color', 'g', 'linewidth', 1);


