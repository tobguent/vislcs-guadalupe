clear all;
figure;

% parameters
stepSize = 0.01;        % integration step size
t0 = 7;                 % start time
duration = 5;           % integration duration
numVerts = 100000;      % number of vertices per line
numLines = 2;           % number of lines

% read data
f = flow('data/cylinder2d.nc');
    
% plot reflectance
slice = (t0 + duration - f.DomainMin(3)) / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
R = f.Reflectance(:,:,slice);
x = linspace(f.AxisX(1), f.AxisX(2), size(R,1));
y = linspace(f.AxisY(1), f.AxisY(2), size(R,2));
h = imagesc(x, y, R');
colormap(gray());
caxis([0 1]);
hold all;

% draw lines
for iLine = 1:numLines
    seed = [0, f.DomainMin(2) + (f.DomainMax(2) - f.DomainMin(2)) * (iLine) / (numLines + 1), t0];
    line = f.streakline(stepSize, duration, seed, numVerts);
    plot(line(:,1), line(:,2), 'LineWidth', 1.1);
    axis equal;
    axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
end

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'white');