% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function output = fastBilateralFilter( image, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange, regularMode, stencilSz )
% Implementation of the approach described in the paper by S. Paris and
% Dorsey (A fast approximation of the bilateral filter using a signal processing approach '07)
% I also used their presentation slides: http://cs.brown.edu/courses/cs129/lectures/bf_course_Brown_Oct2012.pdf
    if ~exist( 'regularMode', 'var' )
        regularMode = 0;
    end;
    
    autoStencilSize = 1;
    if exist( 'stencilSz', 'var' )
        autoStencilSize = 0;
    end;

    sizeImg = size(image);

    minImg = min( image(:) );
    maxImg = max( image(:) );

    spanImg = maxImg - minImg;

    derivedSigmaSpatial = sigmaSpatial / samplingSpatial;
    derivedSigmaRange = sigmaRange / samplingRange;

    % to avoid checking points on borders in gausian convolution, create a 0 band 
    % width and height 'padding' around the grid
    paddingXY = floor(2 * derivedSigmaSpatial) + 1;
    paddingZ = floor(2 * derivedSigmaRange) + 1;
    if (autoStencilSize == 0)
        paddingXY = floor((stencilSz - 1)/2);
        paddingZ = floor((stencilSz - 1)/2);
    end;

    % 1. create grid
    downsampledWidth = floor((sizeImg(2) - 1) / samplingSpatial) + 1 + 2 * paddingXY;
    downsampledHeight = floor((sizeImg(1) - 1) / samplingSpatial) + 1 + 2 * paddingXY;
    downsampledDepth = floor(spanImg / samplingRange) + 1 + 2 * paddingZ;

    grid_wi = zeros(downsampledHeight, downsampledWidth, downsampledDepth);
    grid_w = zeros(downsampledHeight, downsampledWidth, downsampledDepth);

    % 2. downsampling
    for i = 1 : sizeImg(1)
        for j = 1 : sizeImg(2)
            valImg = image(i, j);

            if ~isnan( valImg  ) % skip point in the 0 band
                x = round(i / samplingSpatial) + paddingXY + 1;
                y = round(j / samplingSpatial) + paddingXY + 1;
                ksi = round((valImg - minImg ) / samplingRange) + paddingZ + 1;
            end;

            grid_wi(x, y, ksi) = grid_wi(x, y, ksi) + valImg;
            grid_w(x, y, ksi) = grid_w(x, y, ksi) + 1;
        end;
    end;

    % create gaussian kernel
    kernelWidth = 2 * derivedSigmaSpatial + 1;
    kernelHeight = kernelWidth;
    kernelDepth = 2 * derivedSigmaRange + 1;
    if (autoStencilSize == 0)
        kernelWidth = stencilSz;
        kernelHeight = stencilSz;
        kernelDepth = stencilSz;
    end;

    if (regularMode == 0)
        % convolution with normal gaussian kernel
        grid_wi_b = gausConv3D(kernelWidth, kernelHeight, kernelDepth, ...
            derivedSigmaSpatial, derivedSigmaRange, grid_wi);
        grid_w_b = gausConv3D(kernelWidth, kernelHeight, kernelDepth, ...
            derivedSigmaSpatial, derivedSigmaRange, grid_w);
    else
        % using notation from
        % http://www.cs.technion.ac.il/~elad/publications/conferences/2005/24_ScaleSpace_Retinex_Bilateral.pdf
        % formula 14
        % TODO Instead of regular convolution use a custom which takes only
        % pixels s_ij > s_center
        %nlfilter doesn't work for 3D ;-( 
        %func = @(x) fastBilateralCore(kernel, x);
        %grid_wi_b = nlfilter(grid_wi, [stencilSz stencilSz stencilSz], func);
        %grid_w_b = nlfilter(grid_w, [stencilSz stencilSz stencilSz], func);
    end;

    % 3. divide wi_b / w_b
    normalizedBlurredGrid = zeros(downsampledHeight, downsampledWidth, downsampledDepth);
    for i = 1 : downsampledHeight
        for j = 1 : downsampledWidth
            for k = 1 : downsampledDepth
                if grid_w_b(i, j, k) ~= 0
                    normalizedBlurredGrid(i, j, k) = grid_wi_b(i, j, k) / grid_w_b(i, j, k);
                else
                    normalizedBlurredGrid(i, j, k) = 0;
                end
            end;
        end;
    end;
    
    % 4. upsample
    [ jj, ii ] = meshgrid( 0 : sizeImg(2) - 1, 0 : sizeImg(1) - 1 );

    di = (ii / samplingSpatial) + paddingXY + 1;
    dj = (jj / samplingSpatial) + paddingXY + 1;
    dz = (image - minImg) / samplingRange + paddingZ + 1;

    output = interpn(normalizedBlurredGrid, di, dj, dz);
end

function res = gausConv3D(kernelWidth, kernelHeight, kernelDepth, ...
    derivedSigmaSpatial, derivedSigmaRange, grid)
    kernelX = fspecial('Gaussian', [kernelWidth 1], derivedSigmaSpatial);
    kernelY = fspecial('Gaussian', [1 kernelHeight], derivedSigmaSpatial);
    kernelZ = fspecial('Gaussian', [kernelDepth 1], derivedSigmaRange);

    k = zeros([1, 1, kernelDepth]);
    k(:) = kernelZ;

    data2 = grid;
    data2 = convn(data2, kernelX, 'same');
    data2 = convn(data2, kernelY, 'same');
    data2 = convn(data2, k, 'same');
    res = data2;
end
