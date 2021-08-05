clear all;
figure;

% parameters
stepSize = 100;     % integration step size
t0 = 5000;          % start time in seconds
duration = 500000;  % integration duration in seconds
dtest = 5000;       % termination distance for close streamlines
dsep = 10000;       % separation distance when seeding streamlines

% read data
f = flow('data/guadalupe.nc');

% plot reflectance map
slice = (t0 - f.DomainMin(3)) / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));

[vert, offset, len] = f.evstreamline(stepSize, t0, duration, dtest, dsep);
numLines = length(offset);
    
% trace and draw streamline for each seed point
for iLine = 1:numLines
    from = offset(iLine) + 1;
    to = offset(iLine)+len(iLine);
    zone = repmat('11 R', size(vert(from:to,1),1), 1);
    [lat,lon] = utm2deg(vert(from:to,1), ...
         vert(from:to,2), ...
         zone);
    geoshow(lat, lon, 'LineWidth', 1.5, 'color',[0.7 0.7 1]);
end

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);
