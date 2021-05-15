function dxtide = solid_tide_dehan(YR,MONTH,DATE,hour,min,sec)

% Input YR,MONTH,DATE,hour,minutes and seconds in UTC time, 
% then you will get Observtory of SYSU ocean tide correction
% in ITRF x,y,z in meters.

moon_position_gcrs = planetEphemeris(juliandate(YR, MONTH, DATE, hour, min,  sec),'Earth','Moon','430');
sun_position_gcrs = planetEphemeris(juliandate(YR, MONTH, DATE, hour, min,  sec),'Earth','Sun','430');

FHR = hour + min/60 + sec/3600;

% EOP information
addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
EOP_data = get_eop(YR,MONTH,DATE,hour,min,sec);
addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
% GC2IT = IERS.GCRS2ITRS(MJD_UTC,TAI_UTC,dt,xp,yp);
GC2IT = GCRS2ITRS1(YR,MONTH,DATE,hour,min,sec);
addpath(genpath('E:\潮汐修正\test file\dehanttideinel')) %路径

xsta = [-2358691.210,5410611.484,2410087.607];

% GCRS --> ITRS
xmon = GC2IT*moon_position_gcrs'*1e3;
xsun = GC2IT*sun_position_gcrs'*1e3;
dxtide =dehanttideinel(xsta,YR,MONTH,DATE,FHR,xsun,xmon);
end

