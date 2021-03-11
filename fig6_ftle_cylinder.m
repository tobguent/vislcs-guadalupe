clear all;

% compute
stepSize = 0.1;
duration = 1.5;
slice = 1000;

% read data
f = flow('data/cylinder2d.nc');

% compute forward and backward FTLE
ftleF = f.ftle(stepSize, duration, [1280, 160], slice);
ftleB = f.ftle(-stepSize, duration, [1280, 160], slice);

% define color maps
cmapF = [255 255 255; 254 224 210; 252 187 161; 252 146 114; 251 106 74; 239 59 44; 203 24 29; 165 15 21; 103 0 13] / 255;
xmap = linspace(0,1,size(cmapF,1));
cmapF = interp1(xmap, cmapF, linspace(0,1,128));
cmapB = [255 255 255; 222 235 247; 198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181; 8 81 156; 8 48 107] / 255;
xmap = linspace(0,1,size(cmapB,1));
cmapB = interp1(xmap, cmapB, linspace(0,1,128));


% drawing forward
figure;
x = linspace(f.AxisX(1), f.AxisX(2), f.Resolution(1));
y = linspace(f.AxisY(1), f.AxisY(2), f.Resolution(2));
imagesc(x, y, ftleF');
set(gca,'YDir','normal');
colormap(cmapF);
colorbar;
caxis([0 2.7])
axis equal;
axis([f.AxisX(1)-eps f.AxisX(2)+eps f.AxisY(1)-eps f.AxisY(2)+eps]);
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');


% drawing backward
figure;
x = linspace(f.AxisX(1), f.AxisX(2), f.Resolution(1));
y = linspace(f.AxisY(1), f.AxisY(2), f.Resolution(2));
imagesc(x, y, ftleB');
set(gca,'YDir','normal');
colormap(cmapB);
colorbar;
caxis([0 2.7])
axis equal;
axis([f.AxisX(1)-eps f.AxisX(2)+eps f.AxisY(1)-eps f.AxisY(2)+eps]);
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');