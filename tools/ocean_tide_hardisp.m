function dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour,min,sec)
% Input YR,MONTH,DATE,hour,minutes and seconds in UTC time, 
% then you will get Observtory of SYSU ocean tide correction
% in radial,south,west ---> ITRF x,y,z in meters.

addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
EOP_data = get_eop(YR,MONTH,DATE,hour,min,sec);
addpath(genpath('E:\潮汐修正\test file\hardisp')) %路径
[F,P, TAMP, IDD1]=libiers_tdfrph_call(EOP_data(1),EOP_data(3));
% Load BLQ.txt ZHUHAI 
cto = importdata('inputOTL.txt');
cto_rsw = libiers_hardisp(EOP_data(1),EOP_data(3),cto,F,P, TAMP, IDD1);
% [u s w] = cto_rsw;
% ENU = [-w -s u];
NEU = [-cto_rsw(2) -cto_rsw(3) cto_rsw(1)];
lat = 22.3464/180*pi;
lon = 113.5543/180*pi;
coslon = cos(lon);
sinlon = sin(lon);
coslat = cos(lat);
sinlat = sin(lat);
R = [-sinlat*coslon, -sinlat*sinlon, coslat
-sinlon, coslon, 0
coslat*coslon, coslat*sinlon, sinlat];
% ITRF <--- N E U
dxotide = (inv(R) * NEU')';

end
