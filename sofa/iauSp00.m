%  - - - - - - - -
%   i a u S p 0 0
%  - - - - - - - -
%
%  The TIO locator s', positioning the Terrestrial Intermediate Origin
%  on the equator of the Celestial Intermediate Pole.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     date1,date2      TT as a 2-part Julian Date (Note 1)
%
%  Returned (function value):
%                      the TIO locator s' in radians (Note 2)
%
%  Notes:
%
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
%  2) The TIO locator s' is obtained from polar motion observations by
%     numerical integration, and so is in essence unpredictable.
%     However, it is dominated by a secular drift of about
%     47 microarcseconds per century, which is the approximation
%     evaluated by the present function.
%
%  Reference:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2008 May 24
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp = iauSp00(date1, date2)

constants

% Interval between fundamental epoch J2000.0 and current date (JC).
t = ((date1 - DJ00) + date2) / DJC;

% Approximate s'.
sp = -47e-6 * t * DAS2R;

