classdef flow < handle
    %FLOW Class that stores a vector field.
    properties
        Velocity;       % interleaved velocity components (space-time m/s), single precision
        Resolution;     % resolution of the regular grid (X x Y x T), integer32
        DomainMin;      % min corner of the domain (space m)
        DomainMax;      % max corner of the domain (space m)
        AxisX;          % x-axis range for plotting
        AxisY;          % y-axis range for plotting
        AxisT;          % t-axis range for plotting
        Longitude;      % longitude coordinates of the flow
        Latitude;       % latitude coordinates of the flow
        Contour;        % polyline of a decoration object (points x 2) 
        Reflectance;    % reflectance map (if available)
        LonReflectance; % longitude coordinates of reflectance map (if available)
        LatReflectance; % latitude coordinates of reflectance map (if available)
    end
    
    methods
        % Constructor: Reads the flow from file
        function obj = flow(path)
            obj.Velocity = ncread(path, 'Velocity');
            obj.Resolution = ncread(path, 'Resolution')';
            obj.DomainMin = ncread(path, 'DomainMin')';
            obj.DomainMax = ncread(path, 'DomainMax')';
            obj.Contour = ncread(path, 'Contour');
            try
                obj.AxisX = ncread(path, 'AxisX');
            catch
                obj.AxisX = [obj.DomainMin(1) obj.DomainMax(1)]';
            end
            try
                obj.AxisY = ncread(path, 'AxisY');
            catch
                obj.AxisY = [obj.DomainMin(2) obj.DomainMax(2)]';
            end
            try
                obj.AxisT = ncread(path, 'AxisT');
            catch
                obj.AxisT = [obj.DomainMin(3) obj.DomainMax(3)]';
            end
            try
                res = ncread(path, 'Reflectance_Resolution');
                obj.Reflectance = reshape(ncread(path, 'Reflectance'), res(1), res(2), res(3));
            catch
                obj.Reflectance = [];
            end
            try
                res = ncread(path, 'Reflectance_Resolution');
                obj.LonReflectance = reshape(ncread(path, 'Reflectance_Longitude'), res(1), res(2));
            catch
                obj.LonReflectance = [];
            end
            try
                res = ncread(path, 'Reflectance_Resolution');
                obj.LatReflectance = reshape(ncread(path, 'Reflectance_Latitude'), res(1), res(2));
            catch
                obj.LatReflectance = [];
            end
            try
                obj.Longitude = double(reshape(ncread(path, 'Longitude'), obj.Resolution(1), obj.Resolution(2)));
            catch
                obj.Longitude = [];
            end
            try
                obj.Latitude = double(reshape(ncread(path, 'Latitude'), obj.Resolution(1), obj.Resolution(2)));
            catch
                obj.Latitude = [];
            end
        end
        
        % Gets the U and the V components as separate scalar fields.
        function [U,V] = components(obj)
            res = obj.Resolution(1) * obj.Resolution(2) * obj.Resolution(3);
            vel = reshape(obj.Velocity, 2, res);
            vel = vel';
            U = reshape(vel(:,1), obj.Resolution);
            V = reshape(vel(:,2), obj.Resolution);
        end
        
        % Gets the grid coordinates as separate scalar fields.
        function [X,Y] = grid(obj)
            if isempty(obj.Longitude)
                x = linspace(obj.AxisX(1), obj.AxisX(2), obj.Resolution(1));
                y = linspace(obj.AxisY(1), obj.AxisY(2), obj.Resolution(2));
                [X,Y] = meshgrid(x,y);
                X = X';
                Y = Y';
            else
                X = obj.Longitude;
                Y = obj.Latitude;
            end
        end
        
        % Computes the vorticity of the field. 'slice' is an optional
        % parameter. If it is not specified, vorticity is computed for the
        % entire data set. The result 'vort' is therefore either a slice
        % (2D array) or a volume (3D array)
        function [vort] = vorticity(obj, slice)
            if nargin > 1
                vort = cvorticity(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, slice - 1);
            else
                vort = cvorticity(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity);
                vort = reshape(vort, obj.Resolution);
            end
        end
                
        % Computes the flow map. Particles that leave the domain receive
        % the value 'invalidValue#. 'slice' is an optional parameter. 
        % If it is not specified, FTLE is computed for the entire data set.
        % The result 'f' is therefore either a slice (2D array) or a 
        % volume (3D array).
        function [f] = flowmap(obj, stepSize, duration, invalidValue, out_res, slice)
            if nargin < 5
                out_res = obj.Resolution;
            else
                if length(out_res) == 2
                    out_res = [out_res(1), out_res(2), 1];
                end
            end
            if nargin > 5
                f = cflowmap(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, invalidValue, int32(out_res), slice - 1);
                f = reshape(f, [2,out_res(1),out_res(2)]);
                f = permute(f, [2, 3, 1]);
            else
                f = cflowmap(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, invalidValue, int32(out_res));
                f = reshape(f, [2,out_res(1),out_res(2),out_res(3)]);
                f = permute(f, [2, 3, 4, 1]);
            end
        end
        
        % Computes finite-time Lyapunov exponent. 'slice' is an optional
        % parameter. If it is not specified, FTLE is computed for the
        % entire data set. The result 'f' is therefore either a slice
        % (2D array) or a volume (3D array)
        function [f] = ftle(obj, stepSize, duration, out_res, slice)
            if nargin < 4
                out_res = obj.Resolution;
            else
                if length(out_res) == 2
                    out_res = [out_res(1), out_res(2), 1];
                end
            end
            if nargin > 4
                f = cftle(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, int32(out_res), slice - 1);
            else
                f = cftle(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, int32(out_res));
                f = reshape(f, out_res);
            end
        end
        
        % Compute instantaneous vorticity deviation. 'slice' is an optional
        % parameter. If it is not specified, IVD is computed for the
        % entire data set. The result 'vort' is therefore either a slice
        % (2D array) or a volume (3D array)
        function [vort] = ivd(obj, windowSize, slice)
            if nargin > 2
                vort = civd(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, windowSize, slice - 1);
            else
                vort = civd(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, windowSize);
                vort = reshape(vort, obj.Resolution);
            end
        end
        
        % Compute Lagrangian averaged vorticity deviation. 'slice' is an
        % optional parameter. If it is not specified, LAVD is computed for
        % the entire data set. The result 'f' is therefore either a slice
        % (2D array) or a volume (3D array)
        function [f] = lavd(obj, stepSize, duration, windowSize, out_res, slice)
            if nargin < 5
                out_res = obj.Resolution;
            else
                if length(out_res) == 2
                    out_res = [out_res(1), out_res(2), 1];
                end
            end
            if nargin > 5
                f = clavd(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, windowSize, int32(out_res), slice - 1);
            else
                f = clavd(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, windowSize, int32(out_res));
                f = reshape(f, out_res);
            end
        end
        
        % Trace streamline for a certain duration. The seed point is given
        % in space-time. The result is an Nx3 matrix with N vertices.
        function [f] = streamline(obj, stepSize, duration, seed, numSteps)
            f = cstreamline(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, seed, int32(numSteps));
            f = reshape(f, 3, [])';
        end
        
        % Trace pathline for a certain duration. The seed point is given
        % in space-time. The result is an Nx3 matrix with N vertices.
        function [f] = pathline(obj, stepSize, duration, seed, numSteps)
            f = cpathline(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, seed, int32(numSteps));
            f = reshape(f, 3, [])';
        end
        
        % Trace streakline for a certain duration without refinement. The 
        % seed point is given in space-time. The result is an Nx3 matrix 
        % with N vertices.
        function [f] = streakline(obj, stepSize, duration, seed, numSteps)
            f = cstreakline(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, seed, int32(numSteps));
            f = reshape(f, 3, [])';
        end
        
        % Trace timeline for a certain duration without refinement. The 
        % seed point is given in space-time. The result is an Nx3 matrix 
        % with N vertices.
        function [f] = timeline(obj, stepSize, duration, seedA, seedB, numSteps)
            f = ctimeline(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, duration, seedA, seedB, int32(numSteps));
            f = reshape(f, 3, [])';
        end
        
        % Places streamlines with even spacing using the Jobard-Lefer
        % algorithm for a single time slice at time t0.
        function [vert,offset,length] = evstreamline(obj, stepSize, t0, duration, dtest, dsep)
            [vert,offset,length] = cevstreamline(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, stepSize, t0, duration, dtest, dsep);
            vert = reshape(vert, 2, [])';
        end
        
        % Compute a line integral convolution (LIC). 'slice' is an
        % optional parameter. If it is not specified, LIC is computed for
        % the entire data set. The result 'f' is therefore either a slice
        % (2D array) or a volume (3D array)
        function [f] = lic(obj, stepSize, duration, noise_res, out_res, slice)
            if nargin < 5
                out_res = obj.Resolution;
            else
                if length(out_res) == 2
                    out_res = [out_res(1), out_res(2), 1];
                end
            end
            if nargin > 5
                f = clic(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, int32(noise_res), stepSize, duration, int32(out_res), slice - 1);
            else
                f = clic(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, int32(noise_res), stepSize, duration, int32(out_res));
                f = reshape(f, out_res);
            end
        end
        
        % Advects a texture in the vector field from a slice. 'slice' is an
        % optional parameter. If it is not specified, the texture is
        % advected for the entire data set. The result 'f' is therefore
        % either a slice (2D array) or a volume (3D array)
        function [f] = texadvect(obj, stepSize, duration, tex, out_res, slice)
            if nargin < 5
                out_res = obj.Resolution;
            else
                if length(out_res) == 2
                    out_res = [out_res(1), out_res(2), 1];
                end
            end
            if nargin > 5
                f = ctexadvect(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, single(tex), int32(size(tex)), stepSize, duration, int32(out_res), slice - 1);
            else
                f = ctexadvect(obj.Resolution, obj.DomainMin, obj.DomainMax, obj.Velocity, single(tex), int32(size(tex)), stepSize, duration, int32(out_res));
                f = reshape(f, out_res);
            end
        end
        
    end
end

