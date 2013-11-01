% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function [X,Y] = genRandCurve( dt, distrortionModule, outSize, scaleFactor, doPlot )

    pointCount = 2 * pi/ dt;
  
    randGrid = linspace(0, 2*pi, pointCount);
    splineGrid = linspace(0, 2*pi, outSize);

    % if distrortionModule=0, than it is "smooth" i.e. circle
    % the bigger is distortion the bigger is curvature
    r = scaleFactor*(1 + distrortionModule * rand(size(randGrid)));
    r(end) = r(1); %closed curve
    rr = spline(randGrid, [0 r 0], splineGrid);

    X = rr.*cos(splineGrid);
    Y = rr.*sin(splineGrid);

    if (doPlot == 1)
        plot(X, Y); 
    end
end

