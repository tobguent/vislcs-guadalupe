clear all;
figure;

% parameters
stepSize = 100;
t0 = 0;
duration = 22000;
numVerts = 5000;
numLines = 250;

% read data
f = flow('data/guadalupe.nc');

% plot reflectance map
slice = (t0 - f.DomainMin(3)) / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));
    
% generate random seed points
rng(14,'twister');
seedsX = rand(numLines, 1) * 250000;
seedsY = rand(numLines, 1) * 400000;
seedsX = seedsX + (400000 - seedsY) * 0.5;
seedsX = seedsX + f.DomainMin(1);
seedsY = seedsY + f.DomainMin(2);
    
% trace and draw streakline for each seed point
for iLine = 1:numLines
    seed(1) = seedsX(iLine);
    seed(2) = seedsY(iLine);
    seed(3) = t0;
    line = f.streakline(stepSize, duration, seed, numVerts);
    zone = repmat('11 R', size(line,1), 1);
    [lat,lon] = utm2deg(line(:,1), line(:,2), zone);
    geoshow(lat, lon, 'LineWidth', 1.5, 'color',[0.7 0.7 1]);
end

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);