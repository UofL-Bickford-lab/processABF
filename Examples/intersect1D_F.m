function [Xints] = intersect1D_F(voltage, yVal)
% Finds all X indices where voltage intersects the horizontal line at yVal
numPts = numel(voltage);
yLine(1:numPts, 1) = yVal;
for i = 1:1:numPts
    voltage(i,2) = i;
    yLine(i,2) = i;
end
[X0, ~] = intersections(voltage(:,2), voltage(:,1), yLine(:,2), yLine(:,1));
Xints = unique(round(X0));
end

