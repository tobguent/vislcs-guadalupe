clear all;
figure;

% parameters
slice = 30;

% read data
f = flow('data/guadalupe.nc');
[U,V] = f.components();
[X,Y] = f.grid();

% plot reflectance map
worldmap(f.AxisY, f.AxisX);
R = f.Reflectance(1:int32(3*end/5),:,slice);
Rlat = f.LatReflectance(1:int32(3*end/5),:);
Rlon = f.LonReflectance(1:int32(3*end/5),:);
geoshow(Rlat, Rlon, double(R));

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);