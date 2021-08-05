clear all;
figure;

% parameters
stepSize = 100;
duration = 27000;
numVerts = 5000;
numLines = 30;
t0 = 0;

% read data
f = flow('data/guadalupe.nc');

for iSeed = 1:3
    subplot(1,3,iSeed);
    
    % plot reflectance map
    slice = (t0 + duration - f.DomainMin(3)) / (f.DomainMax(3) - f.DomainMin(3)) * (f.Resolution(3) - 1) + 1;
    worldmap(f.AxisY, f.AxisX);
    R = f.Reflectance(1:int32(3*end/5),:,slice);
    Rlat = f.LatReflectance(1:int32(3*end/5),:);
    Rlon = f.LonReflectance(1:int32(3*end/5),:);
    geoshow(Rlat, Rlon, double(R));

    % define seeding curve
    if iSeed == 1
        minSeedX = 3.398121819096623e+05;
        maxSeedX = 5.212017375388895e+05;
        seedY = 3.133510481207515e+06;
    else
        if iSeed == 2
            minSeedX = 3.563021415123193e+05;
            maxSeedX = 5.376916971415465e+05;
            seedY = 3.103200527365279e+06;
        else
            minSeedX = 3.727921011149763e+05;
            maxSeedX = 5.541816567442035e+05;
            seedY = 3.072890573523041e+06;
        end
    end
    
    % trace and draw streakline for each seed point
    for iLine = 1:numLines
        seed = [minSeedX + (maxSeedX - minSeedX) * (iLine-1.0)/(numLines-1.0), seedY, t0];
        line = f.streakline(stepSize, duration, seed, numVerts);
        zone = repmat('11 R', size(line,1), 1);
        [lat,lon] = utm2deg(line(:,1), line(:,2), zone);
        geoshow(lat, lon, 'LineWidth', 1.5, 'color',[0.7 0.7 1]);
    end

    % draw contour
    geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [1 1 1]);
end