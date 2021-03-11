clear all;
figure;

% compute
res = [1000,1000];      % resolution of the output texture
stepSize = 100;         % integration step size in seconds
duration = 80000;       % integration duration in seconds
slice = 30;             % time slice
noise_res = [250, 250]; % resolution of the noise texture

% read data
f = flow('data/guadalupe.nc');

% compute line integral convolution
lic = f.lic(stepSize, duration, noise_res, res, slice);
lic = (lic - 0.5) * 1.2 + 0.7; % increase contrast
ilat = imresize(f.Latitude, size(lic));
ilon = imresize(f.Longitude, size(lic));
worldmap(f.AxisY, f.AxisX);
geoshow(ilat, ilon, double(lic));

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [0 0 0]);