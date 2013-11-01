% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

function [Vj,Vi] = opticalFlow(images, alpha, iterations)
% implementation of Horn-Shunck flow for arrai of input images
%http://dspace.mit.edu/bitstream/handle/1721.1/6337/AIM-572.pdf?sequence=2

[height, width, frames] = size(images);
 
Vj = zeros(height, width);
Vi = zeros(height, width);
 
for k = 1:frames-1
     
    Ej = zeros(height-1, width-1, frames-1);
    Ei = zeros(height-1, width-1, frames-1);
    Et = zeros(height-1, width-1, frames-1);
     
    % page 6, quantization and noise
    for j = 2:width-1
        for i = 2:height-1
            Ej(i, j, k) = (images(i + 1, j + 1, k) - images(i + 1, j, k) + images(i, j + 1, k)...
                - images(i, j, k) + images(i + 1, j + 1, k + 1) - images(i + 1, j, k + 1)...
                + images(i, j + 1,k+1)-images(i,j,k+1))/4;
             
            Ei(i, j, k) = (images(i, j, k) - images(i + 1, j, k) + images(i, j + 1, k)...
                - images(i + 1, j + 1, k) + images(i, j, k + 1)-images(i + 1, j, k + 1)...
                + images(i, j + 1, k + 1) - images(i + 1, j + 1, k + 1))/4;
             
            Et(i,j,k) = (images(i + 1, j, k + 1) - images(i + 1, j, k) + images(i, j, k + 1)...
                - images(i, j, k) + images(i + 1, j + 1, k + 1)- images(i + 1, j + 1, k)...
                + images(i, j + 1, k + 1) - images(i, j + 1, k))/4;
        end
    end
 
    for nn = 1:iterations
        for j = 2:width-1
            for i = 2:height-1
                 
                % page 6, estimating the laplacian of the flow velocities
                Vjbar = (Vj(i - 1, j) + Vj(i, j + 1) + Vj(i + 1,j) + Vj(i, j - 1))/6 +...
                     (Vj(i - 1,j - 1) + Vj(i - 1,j + 1) + Vj(i + 1,j + 1) + Vj(i + 1,j - 1)) / 12;
                 
                Vibar = (Vi(i - 1,j) + Vi(i,j + 1) + Vi(i + 1,j) + Vi(i, j - 1))/6+...
                    (Vi(i - 1,j - 1) + Vi(i - 1,j + 1) + Vi(i + 1, j + 1) + Vi(i + 1, j - 1)) / 12;
 
                %page 12, iterative solution
                coef = (Ej(i,j,k) * Vjbar + Ei(i,j,k) * Vibar + Et(i,j,k)) ...
                    /(alpha^2 + Ej(i,j,k)^2 + Ei(i,j,k)^2);
                %// update u and v 
                Vj(i,j) = Vjbar - Ej(i,j,k) * coef;
                Vi(i,j) = Vibar - Ei(i,j,k) * coef;
            end
        end
    end
     
end

