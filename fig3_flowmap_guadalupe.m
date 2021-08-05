clear all;
figure;

% parameters
stepSize = 100;
duration = 10000;   % integration duration in seconds
slice = 46;         % index of time slice
invalidValue = -9900;
res = [400, 400];

% read data
f = flow('data/guadalupe.nc');

% compute and map the flow map values to the unit range [0,1] for color coding
phi = f.flowmap(stepSize, duration, invalidValue, res, slice);
image = zeros([res(1),res(2),3]);
image(:,:,1) = (phi(:,:,1) - f.DomainMin(1)) / (f.DomainMax(1) - f.DomainMin(1));
image(:,:,2) = (phi(:,:,2) - f.DomainMin(2)) / (f.DomainMax(2) - f.DomainMin(2));

% set the invalid pixels to white
invalid = 1 - (phi(:,:,1) < -9000) .* (phi(:,:,2) < -9000);
image(:,:,1) = image(:,:,1) .* invalid + 1 - invalid;
image(:,:,2) = image(:,:,2) .* invalid + 1 - invalid;
image(:,:,3) = image(:,:,3) .* invalid + 1 - invalid;

% plot reflectance map
worldmap(f.AxisY, f.AxisX);
geoshow(f.Latitude, f.Longitude, double(image));

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);