% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function output = retinexFilter( image, sigmaSpatial, sigmaRange, ... 
     samplingSpatial, samplingRange, gamma, showPlots )
    % should be applied to the v channel of HSV picture
    % following notations by Elad. Based on paper Retinex by two bilateral
    % filters.
    % http://www.cs.technion.ac.il/~elad/talks/2005/Retinex-ScaleSpace-Short.pdf

    % 1. transform to ln
    image( image == 0 ) = image( image == 0 ) + 0.001; % to avoid inf values for log
    illumination = log(image);
    reflection = illumination;

    % 2. find illumination by filtering with envelope mode
    %illumination = fastBilateralFilter(illumination, sigmaSpatial, sigmaRange, samplingSpatial, samplingRange);
    illumination = regBilateralFilter(illumination, 1, sigmaSpatial, sigmaRange, 15);
    if (showPlots == 1)
        subplot(222);
        imagesc(illumination);
    end;

    % 3. find reflection by filtering with regular mode
    % at this point reflection stores original image
    reflection = (reflection - illumination);
    %reflection = fastBilateralFilter(reflection, sigmaSpatial, sigmaRange, samplingSpatial, samplingRange);
    reflection = regBilateralFilter(reflection, 0, sigmaSpatial, sigmaRange, 5);
    subplot(223);
    imagesc(exp(illumination));

    % 4. apply gamma correction to illumination
    illumination = 1/gamma*illumination; %1.0/gamma*(illumination - log(255)) + log(255); % for [1,255]
    if (showPlots == 1)
        subplot(224);
        imagesc(exp(illumination));
    end;

    % 5. S_res = exp(s_res) = exp( illumination_corrected + reflection )
    output = exp( reflection + illumination);

end
