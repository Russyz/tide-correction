function dxoptide = ocean_ptide(YR,MONTH,DATE,hour,min,sec)

% Input YR,MONTH,DATE,hour,min and sec in UTC time, then
% you will get Observtory of SYSU solid pole tide correction.
global const
SAT_Const
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
G=6.67428e-11;
rho_w = 1025;
u_enu_R = [-0.055815    0.032450   -0.021149];
u_enu_I=[-0.001157    0.049369    0.043607];
g_E = 9.7803278;
m_1 = (eop_data(4) - x_s)/const.Arcs;
m_2 = (y_s - eop_data(5))/const.Arcs;

H_p = sqrt(8 * pi / 15) * const.omega_Earth^2 * const.R_Earth^4 / const.GM_Earth;
K = 4 * pi * G * const.R_Earth * rho_w * H_p / (3 * g_E);
gamma_2_R = 0.6870;
gamma_2_I = 0.0036;
lat = 22.3464/180*pi;
lon = 113.5543/180*pi;
theta = pi / 2 - lat;
coslon = cos(lon);
sinlon = sin(lon);
coslat = cos(lat);
sinlat = sin(lat);
denu = (K * ((m_1 * gamma_2_R + m_2 * gamma_2_I) * u_enu_R + (m_1 * gamma_2_R - m_2 * gamma_2_I) * u_enu_I));

neu = [denu(2),denu(1),denu(3)];
R = [-sinlat*coslon, -sinlat*sinlon, coslat
-sinlon, coslon, 0
coslat*coslon, coslat*sinlon, sinlat];
% ITRF <--- N E U
dxoptide = (inv(R) * neu')';
end