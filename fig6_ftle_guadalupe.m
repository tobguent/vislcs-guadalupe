clear all;

% parameters
stepSize = 100;     % integration step size in seconds
duration = 10000;   % integration duration in seconds
slice = 46;         % index of time slice

% read data
f = flow('data/guadalupe.nc');

% compute forward and backward FTLE    
ftleF = f.ftle(stepSize, duration, [700,640], slice);
ftleB = f.ftle(-stepSize, duration, [700,640], slice);
ilat = imresize(f.Latitude, size(ftleF));
ilon = imresize(f.Longitude, size(ftleF));

% plot reflectance map
figure;
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));
    
% apply color map and draw
[C,A] = colormapping(ftleF, 0, 0.0002, 'white2red');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);
   
% plot reflectance map
figure;
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));
    
% apply color map and draw
[C,A] = colormapping(ftleB, 0, 0.0002, 'white2blue');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);