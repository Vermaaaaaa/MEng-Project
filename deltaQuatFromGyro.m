function dq = deltaQuatFromGyro(omega, dt)
% deltaQuatFromGyro  Compute incremental quaternion from gyro reading.
%
%   omega : 1x3 [wx wy wz] angular rate in rad/s (body frame)
%   dt    : scalar time step in seconds
%   dq    : 1x4 quaternion [w x y z] representing rotation over dt

    % Rotation angle over the interval
    theta = norm(omega) * dt;

    if theta < 1e-12
        % Very small rotation → approximate as identity
        dq = [1 0 0 0];
        return;
    end

    % Unit rotation axis
    u = omega / norm(omega);

    halfTheta = 0.5 * theta;
    s = sin(halfTheta);
    c = cos(halfTheta);

    % Quaternion [w x y z]
    dq = [c, u(1)*s, u(2)*s, u(3)*s];
end
