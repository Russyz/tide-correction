function gamma = refrac_index(phpa,  tc,  rh,  wl)
%  Given:
%    phpa       pressure at the observer (hPa = millibar)
%    tc         ambient temperature at the observer (deg C)
%    rh         relative humidity at the observer (range 0-1)
%    wl         wavelength (micrometers)
%
%  Returned:
%    Refractive index minus 1 at the observer.

% Restrict parameters to safe values.
t = max ( tc, -150 );
t = min ( t, 200 );
p = max ( phpa, 0 );
p = min ( p, 10000 );
r = max ( rh, 0 );
r = min ( r, 1 );
w = max ( wl, 0.1 );
w = min ( w, 1e6 );

% Water vapour pressure at the observer.
if ( p > 0 )
    ps = ( 10^( (0.7859 + 0.03477*t)/(1 + 0.00412*t) ) ) * ...
         ( 1 + p*(4.5e-6 + 6e-10*t*t) );
    pw = r * ps / ( 1 - (1-r)*ps/p );
else
    pw = 0;
end

% Refractive index minus 1 at the observer.
tk = t + 273.15;

    wlsq = w * w;
    gamma = ( ( 77.53484e-6 + ( 4.39108e-7 + 3.666e-9/wlsq ) / wlsq ) * p ...
               - 11.2684e-6*pw ) / tk;

end
