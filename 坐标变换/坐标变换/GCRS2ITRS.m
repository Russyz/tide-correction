% 输入

function GC2IT = GCRS2ITRS(YR,MONTH,DATE,hour,min,sec)

% YR = year(datetime);
% MONTH = month(datetime);
% DATE = day(datetime)-2;
% [hour,  min,  sec] = deal(0);
utc1 = juliandate(YR,MONTH,DATE);
utc2 = (hour * 3600 + min * 60 + sec)/86400;
% addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
[MJD_UTC, TAI_UTC, dt, dx, dy] = get_eop1(YR,MONTH,DATE,hour,min,sec);
% MJD_UTC = EOP_data(1);
% TAI_UTC = EOP_data(2);
% dt = EOP_data(3);
a2r = pi/(180 * 3600);
x_p = dx * a2r;
y_p = dy * a2r;
[UT1,UT2] = iauUtcut1(utc1,utc2,dt);
[TAI1,TAI2] = iauUtctai(utc1,utc2);
[TT1,TT2] = iauTaitt(TAI1,TAI2);

s = py.pysofa.sp00(TT1,TT2);
% W = py.pysofa.rz(-s) * py.pysofa.ry(x_p) * py.pysofa.rx(y_p);
ERA = py.pysofa.era00(UT1,UT2);
% R = py.pysofa.rz(-ERA);
% Q = 
rc2t = py.pysofa.c2i06a(TT1,TT2);
rpom = py.pysofa.pom00(x_p,y_p,s);
GC2IT = py.pysofa.c2tcio(rc2t,ERA,rpom);



end 

