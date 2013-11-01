% (C) Copyright Kirill Lykov 2013.
%
% Distributed under the FreeBSD Software License (See accompanying file license.txt)

countPointInCurve = 1000;
clearvars all;
close all;
figure
subplot(1,2,1)
[x,y] = genRandCurve(0.4188, 2, countPointInCurve, 1, 1);

dx = 0.2;
dy = 0.2;
subplot(1,2,2);
genRandSurface(dx, dy, [x;y], countPointInCurve, 1, 1);
