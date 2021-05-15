function GC2IT = GCRS2ITRS1(YR,MONTH,DATE,hour,min,sec)

% YR = year(datetime);
% MONTH = month(datetime);
% DATE = day(datetime)-2;
% [hour,  min,  sec] = deal(0);
utc1 = juliandate(YR,MONTH,DATE);
utc2 = (hour * 3600 + min * 60 + sec)/86400;
% addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
[MJD_UTC, TAI_UTC, dt, dx, dy] = get_eop1(YR,MONTH,DATE,hour,min,sec);
% [MJD_UTC, TAI_UTC, UT1_UTC, x_p, y_p, dX,dY] = get_eop1(YR,MONTH,DATE,hour,min,sec)
a2r = pi/(180 * 3600);
x_p = dx * a2r;
y_p = dy * a2r;
[UT1,UT2] = iauUtcut1(utc1,utc2,dt);
[TAI1,TAI2] = iauUtctai(utc1,utc2);
[TT1,TT2] = iauTaitt(TAI1,TAI2);
s = iauSp00(TT1,TT2);
ERA = iauEra00(UT1,UT2);
rc2t = iauC2i06a(TT1,TT2);
rpom = iauPom00(x_p,y_p,s);

% GC2IT = iauC2tcio(rc2t,ERA,rpom);
%  Given:
%     tta,ttb           TT as a 2-part Julian Date (Note 1)
%     uta,utb           UT1 as a 2-part Julian Date (Note 1)
%     xp,yp             coordinates of the pole (radians, Note 2)
GC2IT = iauC2t06a(TT1,TT2, UT1,UT2, x_p, y_p);
% GC2IT = inv(Q * rotation(3, -ERA)  * W);
end 


