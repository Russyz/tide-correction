function IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour,min,sec)

% YR = year(datetime);
% MONTH = month(datetime);
% DATE = day(datetime)-2;
% [hour,  min,  sec] = deal(0);
utc1 = juliandate(YR,MONTH,DATE);
utc2 = (hour * 3600 + min * 60 + sec)/86400;
% addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
EOP_data = get_eop(YR,MONTH,DATE,hour,min,sec);
MJD_UTC = EOP_data(1);
TAI_UTC = EOP_data(2);
dt = EOP_data(3);
a2r = pi/(180 * 3600);
x_p = EOP_data(4) * a2r;
y_p = EOP_data(5) * a2r;
[UT1,UT2] = iauUtcut1(utc1,utc2,dt);
[TAI1,TAI2] = iauUtctai(utc1,utc2);
[TT1,TT2] = iauTaitt(TAI1,TAI2);
s = iauSp00(TT1,TT2);
ERA = iauEra00(UT1,UT2);
rc2t = iauC2i06a(TT1,TT2);
rpom = iauPom00(x_p,y_p,s);
IT2GC = double(inv(double(iauC2tcio(rc2t,ERA,rpom))));

end 