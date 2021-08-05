clear all;
figure;

% parameters
stepSize = 0.1;
duration = 1.5;
slice = 1200;
invalidValue = -9900;
res = [1280, 160];

% read data
f = flow('data/cylinder2d.nc');

% compute vorticity
phi = f.flowmap(stepSize, duration, invalidValue, res, slice);

% map the flow map values to the unit range [0,1] for color coding
image = zeros([res(1),res(2),3]);
image(:,:,1) = (phi(:,:,1) - f.AxisX(1)) / (f.AxisX(2) - f.AxisX(1));
image(:,:,2) = (phi(:,:,2) - f.AxisY(1)) / (f.AxisY(2) - f.AxisY(1));

% set the invalid pixels to white
invalid = 1 - (phi(:,:,1) < -9000) .* (phi(:,:,2) < -9000);
image(:,:,1) = image(:,:,1) .* invalid + 1 - invalid;
image(:,:,2) = image(:,:,2) .* invalid + 1 - invalid;
image(:,:,3) = image(:,:,3) .* invalid + 1 - invalid;

% permute the X and Y axis, since imshow has flipped orientation
image = permute(image, [2, 1, 3]);
    
% drawing
RI = imref2d(size(image));
RI.XWorldLimits = f.AxisX';
RI.YWorldLimits = f.AxisY';
imshow(image, RI);
set(gca,'YDir','normal');
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');