%  - - - - - - - - - -
%   i a u G s t 0 0 a
%  - - - - - - - - - -
%
%  Greenwich apparent sidereal time (consistent with IAU 2000
%  resolutions).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     uta,utb        UT1 as a 2-part Julian Date (Notes 1,2)
%     tta,ttb        TT as a 2-part Julian Date (Notes 1,2)
%
%  Returned (function value):
%                    Greenwich apparent sidereal time (radians)
%
%  Notes:
%  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
%     Julian Dates, apportioned in any convenient way between the
%     argument pairs.  For example, JD=2450123.7 could be expressed in
%     any of these ways, among others:
%
%            Part A        Part B
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable (in the case of UT;  the TT is not at all critical
%     in this respect).  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  For UT, the date & time
%     method is best matched to the algorithm that is used by the Earth
%     Rotation Angle function, called internally:  maximum precision is
%     delivered when the uta argument is for 0hrs UT1 on the day in
%     question and the utb argument lies in the range 0 to 1, or vice
%     versa.
%
%  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
%     and TT to predict the effects of precession-nutation.  If UT1 is
%     used for both purposes, errors of order 100 microarcseconds
%     result.
%
%  3) This GAST is compatible with the IAU 2000 resolutions and must be
%     used only in conjunction with other IAU 2000 compatible
%     components such as precession-nutation.
%
%  4) The result is returned in the range 0 to 2pi.
%
%  5) The algorithm is from Capitaine et al. (2003) and IERS
%     Conventions 2003.
%
%  Called:
%     iauGmst00    Greenwich mean sidereal time, IAU 2000
%     iauEe00a     equation of the equinoxes, IAU 2000A
%     iauAnp       normalize angle into range 0 to 2pi
%
%  References:
%     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
%     implement the IAU 2000 definition of UT1", Astronomy &
%     Astrophysics, 406, 1135-1149 (2003)
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2008 May 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gst = iauGst00a(uta, utb, tta, ttb)

gmst00 = iauGmst00(uta, utb, tta, ttb);
ee00a = iauEe00a(tta, ttb);
gst = iauAnp(gmst00 + ee00a);

