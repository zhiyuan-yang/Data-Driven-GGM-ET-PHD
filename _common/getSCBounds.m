function [x,y] = getSCBounds(cx,cy,R,objectBounds)
tmp = R*objectBounds;
x = tmp(1,:) + cx;
y = tmp(2,:) + cy;