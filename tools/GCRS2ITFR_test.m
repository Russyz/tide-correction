addpath(genpath('E:\潮汐修正\坐标变换')) %路径
% set the date [UTC]
MJD  = 53101;

% set the time [UTC]
hour = 7; min = 51;  sec = 28.386009;
% hour = 0; min = 0;  sec = 0;
% compute the seconds of day
sod  = hour*3600 + min*60 + sec;

% EOP information
xp = -0.140682;  yp =  0.333309;  du = -0.439962; dt = 32;

% compute the date+time
fMJD_UTC = MJD +  sod/86400;

% Vallado et al. 2006, AIAA NOTE: using IERS 2003~ [meters]
X_itrs = [-1033.4793830, 7901.2952754, 6380.3565958]';

% dx,dy = 0                                        [meters]
X_gcrs = [5102.5089592, 6123.0114033, 6378.1369247]';

% compute the xformation matrix
GC2IT = IERS.GCRS2ITRS(fMJD_UTC,dt,du,xp,yp);

% perform the transformation
X = GC2IT*X_gcrs
