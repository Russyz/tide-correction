function IT2GC = ITRS2GCRS1(IY,IM,ID,IH,MIN,SEC)
% ===========================================
% IAU 2006/2000A, CIO based, using X,Y series
% ===========================================


constants

% UTC
% IY = 2007;
% IM = 4;
% ID = 5;
% IH = 12;
% MIN = 0;
% SEC = 0;
[MJD_UTC, TAI_UTC, UT1_UTC, x_p, y_p, dX,dY] = get_eop2(IY,IM,ID,IH,MIN,SEC);
% Polar motion (arcsec->radians).
XP = x_p * AS2R;
YP = y_p * AS2R;
% UT1-UTC (s).
DUT1 = UT1_UTC;
% % Nutation corrections wrt IAU 1976/1980 (mas->radians).
% DDP80 = -55.0655 * AS2R/1000;
% DDE80 = -6.3580 * AS2R/1000;
% % CIP offsets wrt IAU 2000A (mas->radians).
% DX00 = 0.1725 * AS2R/1000;
% DY00 = -0.2650 * AS2R/1000;
% CIP offsets wrt IAU 2006/2000A (mas->radians).
DX06 = dX * AS2R;
DY06 = dY * AS2R;

% TT (MJD).
[DJMJD0 DATE] = iauCal2jd(IY, IM, ID);
TIME = ( 60*(60*IH + MIN) + SEC ) / 86400;
UTC = DATE + TIME;
DAT = iauDat ( IY, IM, ID, TIME );
TAI = UTC + DAT/86400;
TT = TAI + 32.184/86400;
% UT1.
TUT = TIME + DUT1/86400;
UT1 = DATE + TUT;

% CIP and CIO, IAU 2006/2000A.
[X Y] = iauXy06 (DJMJD0, TT);
S = iauS06 (DJMJD0, TT, X, Y);
% Add CIP corrections.
X = X + DX06;
Y = Y + DY06;
% GCRS to CIRS matrix.
RC2I = iauC2ixys (X, Y, S);
% Earth rotation angle
ERA = iauEra00 (DJMJD0+DATE, TUT );
% Form celestial-terrestrial matrix (no polar motion yet)
RC2TI = iauCr(RC2I);
RC2TI = iauRz(ERA, RC2TI);
% Polar motion matrix (TIRS->ITRS, IERS 2003)
RPOM = iauPom00 (XP, YP, iauSp00(DJMJD0,TT));
% Form celestial-terrestrial matrix (including polar motion)
RC2IT = iauRxr(RPOM,RC2TI);
IT2GC = inv(RC2IT);
end