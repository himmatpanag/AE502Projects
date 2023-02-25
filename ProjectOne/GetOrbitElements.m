function orbitElements = GetOrbitElements(r0,v0,mu)
    % Inputs:
    %   r0 - position vector [km]
    %   v0 - velocity vector [km/s]
    %   mu - gravitational parameter of the central body [km^3/s^2]
    % Outputs:
    %   a - semi-major axis [km]
    %   e - eccentricity
    %   i - inclination [rad]
    %   Omega - longitude of the ascending node [rad]
    %   omega - argument of periapsis [rad]
    %   M_or_F - mean anomaly [rad] for elliptical orbits, or hyperbolic eccentric anomaly [rad] for hyperbolic orbits

    r = norm(r0); v = norm(v0);

    h = cross(r0,v0);
    p = norm(h)^2/mu;

    eVec = 1/mu * cross(v0,h) - r0./r; 
    e = norm(eVec);

    a = p/(1-e^2);
    inclination = acos(h(3)/norm(h));

    K = [0;0;1];
    n = cross(K,h./norm(h)); 
    
    Omega = acos(n(1)./norm(n));
    
    if n(2) < 0
        Omega = 2*pi - Omega;
    end 

    omega = acos(dot(n,eVec)./norm(n)./norm(eVec));
    if eVec(3) < 0
        omega = 2*pi - omega;
    end 

    if e==0
        omega = 0; 
        f = 0; 
    elseif e < 1 % elliptical orbit 
        f = acos(min(dot(eVec,r0)./norm(eVec)./r,1));
        if dot(r0,v0) < 0
            f = 2*pi - f;
        end 
    else % hyperbolic orbit 
        % Calculate hyperbolic eccentric anomaly from r,v 
        F = acosh((r*norm(e) - dot(r0, v0)) / (norm(e)*r));
        if dot(r0, v0) < 0
            F = -F;
        end
        % Calculate true anomaly of hyberbolic orbit from r,v
        f = 2*atan(sqrt((e+1)/(e-1))*tanh(F/2));
    end

    orbitElements.semiMajorAxis = a; 
    orbitElements.eccentricity = norm(eVec); 
    orbitElements.inclination = inclination;
    orbitElements.argPeriapse = omega;
    orbitElements.longitudeAscendingNode = Omega;
    orbitElements.trueAnomaly = f;
    orbitElements.nodeVector = n;
    orbitElements.eccentricityVec = eVec;
end 