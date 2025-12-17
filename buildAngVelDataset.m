function [angVelBody] = buildAngVelDataset(fs, segmentsAngVel)
% segments is an array of structs with fields:
%   .type, .duration, .params (struct)
dt = 1/fs;

angVelBody = [];

for k = 1:numel(segmentsAngVel)
    N = round(segmentsAngVel(k).duration * fs);
    t = (0:N-1)' * dt;

    w = zeros(N,3);

    p = segmentsAngVel(k).params;

    switch lower(segmentsAngVel(k).type)

        case "static"
            % no motion

        case "const_rot"
            u = p.axisUnit(:).' / norm(p.axisUnit);

            w = repmat(p.omega .* u, N, 1);  % p.omega = [wx wy wz]

        case "sin_rot"
            %Create fixed axis u which the body rotates about
            u = p.axisUnit(:).' / norm(p.axisUnit);
            
            %Magnitude of the sinusoidal rotation
            omegaMag = p.amp * sin(2*pi*p.freq*t + p.phase);
            
            %Apply the magnitude to the axis we are rotating
            w = omegaMag .* u;

        case "chirp_rot"
            u = p.axisUnit(:).' / norm(p.axisUnit);
            
            %Sweep between the 2 frequencies over a time period
            f0 = p.f0; 
            f1 = p.f1; 
            T = segments(k).duration;
            kchirp = (f1 - f0)/T;
            
            %Magnuitude of chirped rotation
            omegaMag = p.amp * sin(2*pi*(f0*t + 0.5*kchirp*t.^2) + p.phase);

            w = omegaMag .* u;

        case "bandlimited_random_rot"
            axis = p.axis;
            x = randn(N,1);
            % simple 1st-order low-pass (cheap + OK for sim)
            alpha = exp(-2*pi*p.cutoff*dt);
            y = zeros(N,1);
            for i=2:N
                y(i) = alpha*y(i-1) + (1-alpha)*x(i);
            end
            w(:,axis) = p.amp * y / std(y);

        case "random_walk_rot"
            axis = p.axis;
            w(:,axis) = cumsum(p.sigma*sqrt(dt)*randn(N,1));

        case "stop_go"
            axis = p.axis;
            % smoothstep from 0 to omega over rampTime, hold, then back
            rampN = round(p.rampTime*fs);
            holdN = max(0, N - 2*rampN);
            s = linspace(0,1,rampN)'; s = s.^2 .* (3 - 2*s); % smoothstep
            prof = [s; ones(holdN,1); flipud(s)] * p.omega;
            prof = prof(1:N);
            w(:,axis) = prof;

        otherwise
            error("Unknown segment type: %s", segments(k).type);
    end

    angVelBody = [angVelBody; w];
end
end
