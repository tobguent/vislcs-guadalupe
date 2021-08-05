clear all;
figure;

% parameters
slice = 15;     % index of time slice

% read data
f = flow('data/guadalupe.nc');

% plot reflectance map
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));
    
% compute vorticity and upsample
vort = -f.vorticity(slice);
vort = interp2(vort,5,'cubic'); % upsample
ilat = interp2(f.Latitude,5,'cubic'); % upsample
ilon = interp2(f.Longitude,5,'cubic'); % upsample
[C,A] = colormapping(vort, -0.0005, 0.0005, 'blue2red');
geoshow(ilat, ilon, C, 'FaceAlpha', 'texturemap', 'AlphaData', A);

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);