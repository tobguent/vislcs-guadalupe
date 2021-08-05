clear all;
figure;

% parameters
stepSize = 0.001;   % integration step size
t0 = 10;            % start time
duration = 9;       % integration duration
numVerts = 200;     % number of vertices per line
numLines = 20;      % number of lines

% read data
f = flow('data/cylinder2d.nc');

% plot reflectance
slice = t0 / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
R = f.Reflectance(:,:,slice);
x = linspace(f.AxisX(1), f.AxisX(2), size(R,1));
y = linspace(f.AxisY(1), f.AxisY(2), size(R,2));
h = imagesc(x, y, R');
colormap(gray());
caxis([0 1]);
hold all;

% draw lines
for iLine = 1:numLines
    seed = [-0.5, f.DomainMin(2) + (f.DomainMax(2) - f.DomainMin(2)) * (iLine - 1.0) / (numLines - 1.0), t0];
    line = f.streamline(stepSize, duration, seed, numVerts);
    plot(line(:,1), line(:,2), 'Color',[0.7 0.7 1], 'LineWidth', 1.1);
    axis equal;
    axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
end

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'white');