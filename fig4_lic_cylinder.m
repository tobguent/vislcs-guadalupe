clear all;
figure;

% parameters
res = [2560, 320];
stepSize = 0.001;
duration = 0.2;
slice = 1000;
noise_res = [320, 40];

% read data
f = flow('data/cylinder2d.nc');

% steady frame approximation
f.Velocity(1:2:end) = f.Velocity(1:2:end) - 0.9;

% compute line integral convolution
lic = f.lic(stepSize, duration, noise_res, res, slice);
lic = (lic - 0.5) * 4 + 0.7; % increase contrast

% drawing
x = linspace(f.AxisX(1), f.AxisX(2), res(1));
y = linspace(f.AxisY(1), f.AxisY(2), res(2));
imagesc(x, y, lic');
set(gca,'YDir','normal');
colormap('gray');
caxis([0 1])
axis equal;
axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');