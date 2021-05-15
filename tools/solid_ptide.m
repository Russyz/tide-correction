function dxptide = solid_ptide(YR,MONTH,DATE,hour,min,sec)

% Input YR,MONTH,DATE,hour,min and sec in UTC time, then
% you will get Observtory of SYSU solid pole tide correction.

eop_data = get_eop(YR,MONTH,DATE,hour,min,sec);
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(eop_data(3),eop_data(2));

MJD_TT = eop_data(1) + TT_UTC/86400;
MJD_TDB = Mjday_TDB(MJD_TT);



t = day(datetime(YR,MONTH,DATE),'dayofyear') ;
% where t is the date in years of 365.25 days.
% Equation (7.21) IERS Conventions 2010
x_s = (55.0 + 1.677 * (t - 2000)) * 1e-3; % arcseconds 
y_s = (320.5 + 3.460 * (t - 2000)) * 1e-3; % arcseconds


%Equation (7.25) IERS Conventions 2010
m_1 = (eop_data(4) - x_s);
m_2 = (y_s - eop_data(5));

lat = 22.3464/180*pi;
lon = 113.5543/180*pi;
theta = pi / 2 - lat;
coslon = cos(lon);
sinlon = sin(lon);
coslat = cos(lat);
sinlat = sin(lat);
% Equation (7.26) IERS Conventions 2010
% south displacement
d_theta = -9 * cos(2 * theta) * (m_1 * coslon + m_2 * sinlon); % Unit.mm
% east displacement
d_lambda = 9 * cos(theta) * (m_1 * coslon - m_2 * sinlon); % Unit.mm
% radial displacement
d_radial = -33 * sin(2 * theta) * (m_1 * coslon + m_2 * sinlon); % Unit.mm

%dxptide = [dup,dsouth,deast] * 1e-3; % Unit.m
% transfrom from coordinate eun to itrs
R1 = [coslat*coslon, coslat*sinlon, -sinlat
    -sinlon, coslon, 0
sinlat*coslon, sinlat*sinlon, coslat];

R = [-sinlat*coslon, -sinlat*sinlon, coslat
-sinlon, coslon, 0
coslat*coslon, coslat*sinlon, sinlat];

dxptide = (inv(R) *(R1 *  [d_theta,d_lambda,d_radial]') * 1e-3)';
%dxptide = ([deast,-dsouth,dup]' * 1e-3)';
end