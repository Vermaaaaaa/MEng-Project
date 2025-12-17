function [accBody] = buildAccDataset(fs, segmentsAcc)
% segments is an array of structs with fields:
%   .type, .duration, .params (struct)
dt = 1/fs;

accBody = [];

for k = 1:numel(segmentsAcc)
    N = round(segmentsAcc(k).duration * fs);
    t = (0:N-1)' * dt;

    a = zeros(N,3);

    p = segmentsAcc(k).params;

    switch lower(segmentsAcc(k).type)

        case "static"
            % no motion
        
        case "const_accel"
            u = p.axisUnit(:).' / norm(p.axisUnit);
            a = repmat(p.accel .* u, N, 1);  % p.accel = [ax ay az]
        
        case "sin_accel"
            axis = p.axis;
            a(:,axis) = p.amp * sin(2*pi*p.freq*t + p.phase);

        otherwise
            error("Unknown segment type: %s", segments(k).type);
    end

    accBody = [accBody; a];
end
end
