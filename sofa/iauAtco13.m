%  - - - - - - - - - -
%   i a u A t c o 1 3
%  - - - - - - - - - -
%
%  ICRS RA,Dec to observed place.  The caller supplies UTC, site
%  coordinates, ambient air conditions and observing wavelength.
%
%  SOFA models are used for the Earth ephemeris, bias-precession-
%  nutation, Earth orientation and refraction.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc,dc     ICRS right ascension at J2000.0 (radians, Note 1)
%     pr        RA proper motion (radians/year; Note 2)
%     pd        Dec proper motion (radians/year)
%     px        parallax (arcsec)
%     rv        radial velocity (km/s, +ve if receding)
%     utc1      UTC as a 2-part...
%     utc2      ...quasi Julian Date (Notes 3-4)
%     dut1      UT1-UTC (seconds, Note 5)
%     elong     longitude (radians, east +ve, Note 6)
%     phi       latitude (geodetic, radians, Note 6)
%     hm        height above ellipsoid (m, geodetic, Notes 6,8)
%     xp,yp     polar motion coordinates (radians, Note 7)
%     phpa      pressure at the observer (hPa = mB, Note 8)
%     tc        ambient temperature at the observer (deg C)
%     rh        relative humidity at the observer (range 0-1)
%     wl        wavelength (micrometers, Note 9)
%
%  Returned:
%     aob       observed azimuth (radians: N=0,E=90)
%     zob       observed zenith distance (radians)
%     hob       observed hour angle (radians)
%     dob       observed declination (radians)
%     rob       observed right ascension (CIO-based, radians)
%     eo        equation of the origins (ERA-GST)
%
%  Notes:
%  1)  Star data for an epoch other than J2000.0 (for example from the
%      Hipparcos catalog, which has an epoch of J1991.25) will require
%      a preliminary call to iauPmsafe before use.
%
%  2)  The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
%
%  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
%      convenient way between the two arguments, for example where utc1
%      is the Julian Day Number and utc2 is the fraction of a day.
%
%      However, JD cannot unambiguously represent UTC during a leap
%      second unless special measures are taken.  The convention in the
%      present function is that the JD day represents UTC days whether
%      the length is 86399, 86400 or 86401 SI seconds.
%
%      Applications should use the function iauDtf2d to convert from
%      calendar date and time of day into 2-part quasi Julian Date, as
%      it implements the leap-second-ambiguity convention just
%      described.
%
%  4)  The warning status "dubious year" flags UTCs that predate the
%      introduction of the time scale or that are too far in the
%      future to be trusted.  See iauDat for further details.
%
%  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
%      one second at the end of each positive UTC leap second,
%      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
%      practice is under review, and in the future UT1-UTC may grow
%      essentially without limit.
%
%  6)  The geographical coordinates are with respect to the WGS84
%      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
%      longitude required by the present function is east-positive
%      (i.e. right-handed), in accordance with geographical convention.
%
%  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
%      values are the coordinates (in radians) of the Celestial
%      Intermediate Pole with respect to the International Terrestrial
%      Reference System (see IERS Conventions 2003), measured along the
%      meridians 0 and 90 deg west respectively.  For many
%      applications, xp and yp can be set to zero.
%
%  8)  If hm, the height above the ellipsoid of the observing station
%      in meters, is not known but phpa, the pressure in hPa (=mB),
%      is available, an adequate estimate of hm can be obtained from
%      the expression
%
%            hm = -29.3 * tsl * log ( phpa / 1013.25 );
%
%      where tsl is the approximate sea-level air temperature in K
%      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
%      52).  Similarly, if the pressure phpa is not known, it can be
%      estimated from the height of the observing station, hm, as
%      follows:
%
%            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
%
%      Note, however, that the refraction is nearly proportional to
%      the pressure and that an accurate phpa value is important for
%      precise work.
%
%  9)  The argument wl specifies the observing wavelength in
%      micrometers.  The transition from optical to radio is assumed to
%      occur at 100 micrometers (about 3000 GHz).
%
%  10) The accuracy of the result is limited by the corrections for
%      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
%      Providing the meteorological parameters are known accurately and
%      there are no gross local effects, the predicted observed
%      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
%      (radio) for a zenith distance of less than 70 degrees, better
%      than 30 arcsec (optical or radio) at 85 degrees and better
%      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
%
%      Without refraction, the complementary functions iauAtco13 and
%      iauAtoc13 are self-consistent to better than 1 microarcsecond
%      all over the celestial sphere.  With refraction included,
%      consistency falls off at high zenith distances, but is still
%      better than 0.05 arcsec at 85 degrees.
%
%  11) "Observed" Az,ZD means the position that would be seen by a
%      perfect geodetically aligned theodolite.  (Zenith distance is
%      used rather than altitude in order to reflect the fact that no
%      allowance is made for depression of the horizon.)  This is
%      related to the observed HA,Dec via the standard rotation, using
%      the geodetic latitude (corrected for polar motion), while the
%      observed HA and RA are related simply through the Earth rotation
%      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
%      means the position that would be seen by a perfect equatorial
%      with its polar axis aligned to the Earth's axis of rotation.
%
%  12) It is advisable to take great care with units, as even unlikely
%      values of the input parameters are accepted and processed in
%      accordance with the models used.
%
%  Called:
%     iauApco13    astrometry parameters, ICRS-observed, 2013
%     iauAtciq     quick ICRS to CIRS
%     iauAtioq     quick ICRS to observed
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aob, zob, hob, dob, rob, eo] = iauAtco13(rc, dc, pr, pd, px, ...
          rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)

astrom.pmt = 0;         % PM time interval (SSB, Julian years)
astrom.eb = zeros(3,1); % SSB to observer (vector, au)
astrom.eh = zeros(3,1); % Sun to observer (unit vector)
astrom.em = 0;          % distance from Sun to observer (au)
astrom.v = zeros(3,1);  % barycentric observer velocity (vector, c)
astrom.bm1 = 0;         % sqrt(1-|v|^2): reciprocal of Lorenz factor
astrom.bpn = zeros(3);  % bias-precession-nutation matrix
astrom.along = 0;       % longitude + s' + dERA(DUT) (radians)
astrom.phi = 0;         % geodetic latitude (radians)
astrom.xpl = 0;         % polar motion xp wrt local meridian (radians)
astrom.ypl = 0;         % polar motion yp wrt local meridian (radians)
astrom.sphi = 0;        % sine of geodetic latitude
astrom.cphi = 0;        % cosine of geodetic latitude
astrom.diurab = 0;      % magnitude of diurnal aberration vector
astrom.eral = 0;        % "local" Earth rotation angle (radians)
astrom.refa = 0;        % refraction constant A (radians)
astrom.refb = 0;        % refraction constant B (radians)

% Star-independent astrometry parameters.
[astrom, eo] = iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);

% Transform ICRS to CIRS.
[ri, di] = iauAtciq(rc, dc, pr, pd, px, rv, astrom);

% Transform CIRS to observed.
[aob, zob, hob, dob, rob] = iauAtioq(ri, di, astrom);

