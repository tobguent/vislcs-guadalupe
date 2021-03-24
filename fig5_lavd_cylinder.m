clear all;
figure;

% parameters
stepSize = 0.1;
duration = 2;
windowSize = 20;
slice = 1200;

% read data
f = flow('data/cylinder2d.nc');

% compute Lagrangian average vorticity deviation
lavd = f.lavd(stepSize, duration, windowSize, [1280, 160], slice);

% create color map
cmap = [33 102 172; 67 147 195;146 197 222; 209 229 240; 247 247 247; 253 219 199; 244 165 130; 214 96 77; 178 24 43] / 255;
xmap = linspace(0,1,size(cmap,1));
cmap = interp1(xmap, cmap, linspace(0,1,128));

% drawing
x = linspace(f.AxisX(1), f.AxisX(2), f.Resolution(1));
y = linspace(f.AxisY(1), f.AxisY(2), f.Resolution(2));
imagesc(x, y, lavd');
set(gca,'YDir','normal');
colormap(cmap);
colorbar;
caxis([-10 10])
axis equal;
axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');