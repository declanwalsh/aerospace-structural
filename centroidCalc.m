% function that calculates the centroid of a idealised structure
% B is boom area
% x, y are boom co-ordinates with respect to an arbitrary origin
% centroid is returned relative to this same origin
function [ centroid ] = centroidCalc( B, x, y )

    xC = sum(B.*x)/sum(B);
    yC = sum(B.*y)/sum(B);
    centroid = [xC, yC];

end

