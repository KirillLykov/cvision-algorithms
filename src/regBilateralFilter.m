% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function filteredImg = regBilateralFilter(inputImage, isRetinex, sigmaD, sigmaR, stencilSz)
    % Implementation of the Bilateral Filter by C. Tomasi and R. Manduchi
    % "Bilateral filtering for gray and color images"
    % isRetinex flag is used for purpose of retinex algorithm
    func = @(x) bilateralCore(isRetinex, sigmaD, sigmaR, x);
    filteredImg = nlfilter(inputImage, [stencilSz stencilSz], func);
end

function res = bilateralCore(isRetinex, sigmaD, sigmaR, pixelVicinity)
    
    sigmaD2 = sigmaD*sigmaD;
    sigmaR2 = sigmaR*sigmaR;

    dist2 = @(x,y) double((x(1) - y(1))^2 + (x(2) - y(2))^2); 

    sz = size(pixelVicinity);
    center = int32(sz/2);
    f_center = pixelVicinity(center(1), center(2));
    
    totalMass = 0.0;
    filteredVal = 0.0;
    for i=1:sz(1)
        for j = 1:sz(2)
            
            spatialDist = exp(-(f_center - pixelVicinity(i,j))^2/2/sigmaR2);
            if (f_center > pixelVicinity(i,j) && isRetinex == 1)
                spatialDist = 0;
            end;
            w = exp(-dist2(center, [i, j])/2/sigmaD2) * spatialDist;
            
            totalMass = totalMass + w;
            
            filteredVal = filteredVal + pixelVicinity(i, j) * w;
        end;
    end;
    
    res = filteredVal / totalMass;
end