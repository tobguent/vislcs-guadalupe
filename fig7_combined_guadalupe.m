clear all;
figure;

% parameters
stepSize = 100;         % integration step size in seconds
durationFTLE = 10000;   % integration duration for FTLE
durationLAVD = 10000;   % integration duration for LAVD
windowSize = 10;        % neighborhood size for LAVD
slice = 46;             % index of time slice

% read data
f = flow('data/guadalupe.nc');
    
% plot reflectance map
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));
    
% compute LAVD
lavd = -f.lavd(stepSize, durationLAVD, windowSize, [700,640], slice);
ilat = imresize(f.Latitude, size(lavd));
ilon = imresize(f.Longitude, size(lavd));
[C,A] = colormapping(lavd, -0.0004, 0.0004, 'blue2red');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);
    
% compute FTLE
ftleB = f.ftle(-stepSize, durationFTLE, [700,640], slice);
ilat = imresize(f.Latitude, size(ftleB));
ilon = imresize(f.Longitude, size(ftleB));
[C,A] = colormapping(ftleB, 0, 0.0003, 'yellow');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);
    
% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);