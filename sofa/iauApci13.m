%  - - - - - - - - - -
%   i a u A p c i 1 3
%  - - - - - - - - - -
%
%  For a terrestrial observer, prepare star-independent astrometry
%  parameters for transformations between ICRS and geocentric CIRS
%  coordinates.  The caller supplies the date, and SOFA models are used
%  to predict the Earth ephemeris and CIP/CIO.
%
%  The parameters produced by this function are required in the
%  parallax, light deflection, aberration, and bias-precession-nutation
%  parts of the astrometric transformation chain.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1        TDB as a 2-part...
%     date2        ...Julian Date (Note 1)
%
%  Returned:
%     astrom         star-independent astrometry parameters:
%      pmt           PM time interval (SSB, Julian years)
%      eb(3)         SSB to observer (vector, au)
%      eh(3)         Sun to observer (unit vector)
%      em            distance from Sun to observer (au)
%      v(3)          barycentric observer velocity (vector, c)
%      bm1           sqrt(1-|v|^2): reciprocal of Lorenz factor
%      bpn(3,3)      bias-precession-nutation matrix
%      along         unchanged
%      xpl           unchanged
%      ypl           unchanged
%      sphi          unchanged
%      cphi          unchanged
%      diurab        unchanged
%      eral          unchanged
%      refa          unchanged
%      refb          unchanged
%      eo             equation of the origins (ERA-GST)
%
%  Notes:
%  1) The TDB date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TDB)=2450123.7 could be expressed in any of these ways, among
%     others:
%
%            date1          date2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in cases
%     where the loss of several decimal digits of resolution is
%     acceptable.  The J2000 method is best matched to the way the
%     argument is handled internally and will deliver the optimum
%     resolution.  The MJD method and the date & time methods are both
%     good compromises between resolution and convenience.  For most
%     applications of this function the choice will not be at all
%     critical.
%
%     TT can be used instead of TDB without any significant impact on
%     accuracy.
%
%  2) All the vectors are with respect to BCRS axes.
%
%  3) In cases where the caller wishes to supply his own Earth
%     ephemeris and CIP/CIO, the function iauApci can be used instead
%     of the present function.
%
%  4) This is one of several functions that inserts into the astrom
%     structure star-independent parameters needed for the chain of
%     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
%
%     The various functions support different classes of observer and
%     portions of the transformation chain:
%
%          functions         observer        transformation
%
%       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
%       iauApci iauApci13    terrestrial     ICRS <-> CIRS
%       iauApco iauApco13    terrestrial     ICRS <-> observed
%       iauApcs iauApcs13    space           ICRS <-> GCRS
%       iauAper iauAper13    terrestrial     update Earth rotation
%       iauApio iauApio13    terrestrial     CIRS <-> observed
%
%     Those with names ending in "13" use contemporary SOFA models to
%     compute the various ephemerides.  The others accept ephemerides
%     supplied by the caller.
%
%     The transformation from ICRS to GCRS covers space motion,
%     parallax, light deflection, and aberration.  From GCRS to CIRS
%     comprises frame bias and precession-nutation.  From CIRS to
%     observed takes account of Earth rotation, polar motion, diurnal
%     aberration and parallax (unless subsumed into the ICRS <-> GCRS
%     transformation), and atmospheric refraction.
%
%  5) The context structure astrom produced by this function is used by
%     iauAtciq* and iauAticq*.
%
%  Called:
%     iauEpv00     Earth position and velocity
%     iauPnm06a    classical NPB matrix, IAU 2006/2000A
%     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
%     iauS06       the CIO locator s, given X,Y, IAU 2006
%     iauApci      astrometry parameters, ICRS-CIRS
%     iauEors      equation of the origins, given NPB matrix and s
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [astrom, eo] = iauApci13(date1, date2)

% Earth barycentric & heliocentric position/velocity (au, au/d).
[ehpv, ebpv] = iauEpv00(date1, date2);

% Form the equinox based BPN matrix, IAU 2006/2000A.
r = iauPnm06a(date1, date2);

% Extract CIP X,Y.
[x, y] = iauBpn2xy(r);

% Obtain CIO locator s.
s = iauS06(date1, date2, x, y);

% Compute the star-independent astrometry parameters.
astrom = iauApci(date1, date2, ebpv, ehpv(1,:), x, y, s);

% Equation of the origins.
eo = iauEors(r, s);

