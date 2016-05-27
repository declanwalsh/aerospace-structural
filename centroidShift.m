% co-odinates x and y are input relative to an arbitray origin that the centroid is defined from
% returns co-ordinates relative to the centroid
function [ xB, yB ] = centroidShift( x, y, centroid )

    xB = x - centroid(1);
    yB = y - centroid(2);

end

