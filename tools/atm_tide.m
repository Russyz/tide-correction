function dxatide = atm_tide(YR,MONTH,DATE,hour,min,sec)
% Input YR,MONTH,DATE,hour,minutes and seconds in UTC time, 
% then you will get Observtory of SYSU atmosphere tide correction
% in radial,south,west ---> ITRF x,y,z in meters.
fid = fopen("E:\潮汐修正\IERS模型\海洋潮\s1_s2.txt");
s1_s2 = fscanf(fid,'%f %f %f %f',[4 Inf]);
fclose(fid);
addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
EOP_data = get_eop(YR,MONTH,DATE,hour,min,sec);
% t is UT1 in days
t = (hour + min/60 + (sec + EOP_data(3))/3600)/24;
omega_1 = 2 * pi;
omega_2 = 4 * pi;
% row 1: radial displacement (mm): dr(1), dr(2), dr(3), dr(4)
% row 2: tangential NS displacement (mm): dn(1), dn(2), dn(3), dn(4)
% row 3: tangential EW displacement (mm): de(1), de(2), de(3), de(4).
coeff_A_d1_u = (s1_s2(1,1));
coeff_B_d1_u = (s1_s2(2,1));
coeff_A_d2_u = (s1_s2(3,1));
coeff_B_d2_u = (s1_s2(4,1));

coeff_A_d1_n = (s1_s2(1,2));
coeff_B_d1_n = (s1_s2(2,2));
coeff_A_d2_n = (s1_s2(3,2));
coeff_B_d2_n = (s1_s2(4,2));

coeff_A_d1_e = (s1_s2(1,3));
coeff_B_d1_e = (s1_s2(2,3));
coeff_A_d2_e = (s1_s2(3,3));
coeff_B_d2_e = (s1_s2(4,3));

de = coeff_A_d1_e * cos(omega_1 * t)...
+ coeff_B_d1_e * sin(omega_1 * t)...
+ coeff_A_d2_e * cos(omega_2 * t)...
+ coeff_B_d2_e * sin(omega_2 * t);
dn = coeff_A_d1_n * cos(omega_1 * t)...
+ coeff_B_d1_n * sin(omega_1 * t)...
+ coeff_A_d2_n * cos(omega_2 * t)...
+ coeff_B_d2_n * sin(omega_2 * t);
du = coeff_A_d1_u * cos(omega_1 * t)...
+ coeff_B_d1_u  * sin(omega_1 * t)...
+ coeff_A_d2_u * cos(omega_2 * t)...
+ coeff_B_d2_u * sin(omega_2 * t);
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
dxatide = (inv(R) * [dn,de,du]' * 1e-3)';
end
