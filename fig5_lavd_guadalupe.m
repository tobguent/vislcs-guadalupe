clear all;
figure;

% parameters
stepSize = 100;     % integration step size in seconds
duration = 10000;   % integration duration in seconds
windowSize = 10;    % window size for spatial averaging
slice = 15;         % index of time slice

% read data
f = flow('data/guadalupe.nc');

% plot reflectance map
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));

% compute Lagrangian averaged vorticity deviation
lavd = -f.lavd(stepSize, duration, windowSize, [700,640], slice);
ilat = imresize(f.Latitude, size(lavd));
ilon = imresize(f.Longitude, size(lavd));
[C,A] = colormapping(lavd, -0.0005, 0.0005, 'blue2red');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);