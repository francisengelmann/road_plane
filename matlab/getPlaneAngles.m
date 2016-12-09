function [pitch, roll] = getPlaneAngles(pnormal)
    % If the plane is flat on the ground the normal will be [0 1 0]
    % get the relative pitch (z slope) and roll (x slope) of the plane
    groundn = [0 1 0];
    % Get the yz component of the normal
    yz = [0 pnormal(2:3)];
    normyz = norm(yz);
    if normyz > 0
        yz = yz ./ normyz;
        % Get the angle between it and the ground normal [0 1 0]
        % (dot prod same as yz(2)), magnitude of both are 1 so no division
        % Use absolute value so we only get between 0-90 not 0-180 with 179
        % being the same as 1
        pitch = real(acosd(abs(yz * groundn')));
    else
        pitch = 0;
    end
    % Now get the twist to the side
    % Sometimes precision error produces a small imaginary number so get
    % the real only
    roll = real(acosd(abs(pnormal * yz' / norm(pnormal))));
    
    if yz(2) * yz(3) > 0
        % The y and z components have the same sign, meaning the plane
        % is tilted forward
        pitch = pitch * -1;
    end
    if pnormal(1) * yz(2) < 0
        % Plane is tilted to the left
        roll = roll * -1;
    end
end
