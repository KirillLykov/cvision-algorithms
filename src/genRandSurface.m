% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function [x, y, z, backwardMapping] = genRandSurface(dx, dy, coord, countPointInCurve, distrortionModule, doPlot)
    %if distrortionModule==0 than it is level set
    %as domain use 2 x bounding box
    xmin = 2 * floor(min(coord(1,:)));
    xmax = 2 * ceil(max(coord(1,:)));
    ymin = 2 * floor(min(coord(2,:)));
    ymax = 2 * ceil(max(coord(2,:)));

    %construct 0,1 map. for every voxel it indicates whether there is a curve
    %point inside or not.
    dim = [fix((xmax - xmin)/dx) + 1, fix((ymax - ymin)/dy) + 1];
    %[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
    curveMask = zeros(dim(1), dim(2));
    backwardMapping = zeros(1, countPointInCurve);
    for i = 1:countPointInCurve
        indX = fix((coord(1, i) - xmin)/dx);
        indY = fix((coord(2, i) - ymin)/dy);

        backwardMapping(1, i) = indX;
        backwardMapping(2, i) = indY;
        
        curveMask(indX, indY) = 1;
    end

    signMask = imfill(int32(curveMask), 'holes');

    signMask = -2 * signMask + 1;
    
    %subplot(2,2,2);
    %imagesc(curveMask);

    %subplot(2,2,3);
    %imagesc(signMask);

    [x, y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
    z = double(bwdist(curveMask));
    z = z.*double(signMask);
    
    if distrortionModule ~= 0
        negativeCurveMask = 1 - curveMask;
        noise = distrortionModule * rand(size(z));
        noise = noise.*negativeCurveMask;
        z = z + noise;
    end;
    
    % it must be dx == dy, otherwise it is not convinient to use bwdist
    z = dx*z';
    
    if doPlot == 1
        surf(x, y, z);
    end;
end