%  - - - - - - - - - -
%   i a u N u m 0 6 a
%  - - - - - - - - - -
%
%  Form the matrix of nutation for a given date, IAU 2006/2000A model.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1,date2             TT as a 2-part Julian Date (Note 1)
%
%  Returned:
%     rmatn                   nutation matrix
%
%  Notes:
%  1) The TT date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TT)=2450123.7 could be expressed in any of these ways,
%     among others:
%
%            date1          date2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable.  The J2000 method is best matched to the way
%     the argument is handled internally and will deliver the
%     optimum resolution.  The MJD method and the date & time methods
%     are both good compromises between resolution and convenience.
%
%  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
%     the p-vector V(true) is with respect to the true equatorial triad
%     of date and the p-vector V(mean) is with respect to the mean
%     equatorial triad of date.
%
%  Called:
%     iauObl06     mean obliquity, IAU 2006
%     iauNut06a    nutation, IAU 2006/2000A
%     iauNumat     form nutation matrix
%
%  Reference:
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992),
%     Section 3.222-3 (p114).
%
%  This revision:  2008 May 12
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rmatn = iauNum06a(date1, date2)

% Mean obliquity.
eps = iauObl06(date1, date2);

% Nutation components.
[dp, de] = iauNut06a(date1, date2);

% Nutation matrix.
rmatn = iauNumat(eps, dp, de);

