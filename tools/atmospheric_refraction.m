function correction = atmospheric_refraction(true_elevation, baroPressure, temperature)
% INPUT true_elevation in [degree], baroPressure in [Pa], temperatire in
% [celsius degrees]. Return the atmospheric refraction for elevation in
% [degree].
% See https://en.wikipedia.org/wiki/Atmospheric_refraction

factor = (baroPressure / 101000.0) * (283.0 / (273.0 + temperature));
ang = true_elevation + 10.3 / (true_elevation + 5.11); % in degree
R = 1.02 / tan(deg2rad(ang)) * factor; % in arcmin
correction = (R / 60.0); % in degree
end