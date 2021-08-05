clear all;
figure;

% parameters
stepSize = 0.001;   % integration step size
t0 = 10;            % start time
duration = 9;       % integration duration
dtest = 0.025;      % termination distance for close streamlines
dsep = 0.05;        % separation distance when seeding streamlines

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

[vert, offset, len] = f.evstreamline(stepSize, t0, duration, dtest, dsep);
numLines = length(offset);

% draw lines
for iLine = 1:numLines
    from = offset(iLine) + 1;
    to = offset(iLine)+len(iLine);
    plot(vert(from:to,1), ...
         vert(from:to,2), ...
         'Color',[0.7 0.7 1], 'LineWidth', 1.1);
    axis equal;
    axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
end

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'white');