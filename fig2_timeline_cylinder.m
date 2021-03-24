clear all;
figure;

% parameters
stepSize = 0.001;       % integration step size
t0 = 8;                 % start time
duration = 1.2;         % integration duration
numVerts = 2000;        % number of vertices per line
numLines = 2;           % number of lines

% read data
f = flow('data/cylinder2d.nc');

% plot reflectance
slice = (t0 + duration f.DomainMin(3)) / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
R = f.Reflectance(:,:,slice);
x = linspace(f.AxisX(1), f.AxisX(2), size(R,1));
y = linspace(f.AxisY(1), f.AxisY(2), size(R,2));
h = imagesc(x, y, R');
colormap(gray());
caxis([0 1]);
hold all;

% draw lines
for iLine = 1:numLines
    seedA = [-0.5, f.DomainMin(2) + (f.DomainMax(2) - f.DomainMin(2)) * (iLine) / (numLines + 1), t0];
    seedB = [ 6.0, f.DomainMin(2) + (f.DomainMax(2) - f.DomainMin(2)) * (iLine) / (numLines + 1), t0];
    line = f.timeline(stepSize, duration, seedA, seedB, numVerts);
    plot(line(:,1), line(:,2), 'LineWidth', 1.1);
    axis equal;
    axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
    hold all;
end

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'white');