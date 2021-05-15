%       SUBROUTINE ADMINT (AMPIN,IDTIN,PHIN,AMP,F,P,NIN,NOUT)
% *+
% *  - - - - - - - - - - -
% *   A D M I N T
% *  - - - - - - - - - - -
% *
% *  This routine is part of the International Earth Rotation and
% *  Reference Systems Service (IERS) Conventions software collection.
% *
% *  This subroutine returns the ocean loading displacement amplitude,
% *  frequency, and phase of a set of tidal constituents generated by
% *  the Bos-Scherneck website at http://www.oso.chalmers.se/~loading/.
% *  The variable nin is input as the number wanted, and the variable 
% *  nout is returned as the number provided.  The constituents used
% *  are stored in the arrays idd (Doodson number) and tamp
% *  (Cartwright-Edden amplitude).  The actual amp and phase of each
% *  of these are determined by spline interpolation of the real and
% *  imaginary part of the admittance, as specified at a subset of the
% *  constituents.
% *
% *  In general, Class 1, 2, and 3 models represent physical effects that
% *  act on geodetic parameters while canonical models provide lower-level
% *  representations or basic computations that are used by Class 1, 2, or
% *  3 models.
% * 
% *  Status:  Class 1 model
% *
% *     Class 1 models are those recommended to be used a priori in the
% *     reduction of raw space geodetic data in order to determine
% *     geodetic parameter estimates.
% *     Class 2 models are those that eliminate an observational
% *     singularity and are purely conventional in nature.
% *     Class 3 models are those that are not required as either Class
% *     1 or 2.
% *     Canonical models are accepted as is and cannot be classified as
% *     a Class 1, 2, or 3 model.
% *
% *  Given:
% *     AMPIN       d      Cartwright-Edden amplitude of tidal constituents
% *     IDTIN       i      Doodson number of tidal constituents
% *     PHIN        d      Phase of tidal constituents
% *     NIN         i      Number of harmonics used
% *
% *  Returned:
% *     AMP         d      Amplitude due to ocean loading
% *     F           d      Frequency due to ocean loading
% *     P           d      Phase due to ocean loading
% *     NOUT        i      Number of harmonics returned
% *
% *  Notes:
% *
% *  1) The phase is determined for a time set in COMMON block /date/ in
% *     the subroutine TDFRPH.
% *  
% *  2) The arrays F and P must be specified as double precision. 
% *
% *  Called:
% *     TDFRPH             Returns frequency and phase of a tidal
% *                        constituent with given Doodson number            
% *     SPLINE             Sets up array for cubic spline interpolation
% *     EVAL               Performs cubic spline interpolation 
% *     SHELLS             Sorts an array using Shell Sort
% *
% *  Test case:
% *     Test cases are provided in the main program HARDISP.F.
% *
% *  References:
% *
% *     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
% *     IERS Technical Note No. 32, BKG (2004)
% *
% *  Revisions:  
% *  2009 June 17 B.E. Stetzler  Initial changes to header
% *  2009 June 18 B.E. Stetzler  Used IMPLICIT NONE, declared more variables,
% *                              and added D0 to DOUBLE PRECISION variables 
% *  2009 August 19 B.E.Stetzler Capitalized all variables for FORTRAN 77
% *                              compatibility
% *-----------------------------------------------------------------------
% 
%       IMPLICIT NONE
%       INTEGER I,IDD,IDTIN,J,K,KEY,LL,NCON,NDI,NIN,NLP,NOUT,NSD,NT,II,KK
% 
% *+----------------------------------------------------------------------
% *  The parameters below set the number of harmonics used in the prediction
% *  (nt; This must also be set in the main program) and the number of
% *  constituents whose amp and phase may be specified (ncon)
% *-----------------------------------------------------------------------
%       PARAMETER (NT=342)
%       PARAMETER (NCON=20)
% 
%       REAL AIM,AMP,AMPIN,DI,DR,DTR,PHIN,RF,RL,SCR,SDI,SDR,TAMP,ZDI,ZDR,
%      .     EVAL,AM,RE,SF
%       DOUBLE PRECISION F,FR,P,PR
% 
%       DIMENSION AMPIN(*),IDTIN(6,*),PHIN(*)
%       DIMENSION AMP(*),F(*),P(*)
% 
% *  Arrays containing information about all stored constituents
%       DIMENSION IDD(6,NT),TAMP(NT)
% 
% *  Arrays containing information about the subset whose amp and phase may
% *  be specified, and scratch arrays for the spline routines for which
% *  at most ncon constituents may be specified.
% 
%       DIMENSION RL(NCON),AIM(NCON),RF(NCON),KEY(NCON),SCR(NCON),
%      . ZDI(NCON),ZDR(NCON),DI(NCON),DR(NCON),SDI(NCON),SDR(NCON)
%       DATA DTR/.01745329252/
%       DATA RL/NCON*0.0/,AIM/NCON*0.0/,RF/NCON*0.0/
%       DATA ZDI/NCON*0.0/,ZDR/NCON*0.0/,DI/NCON*0.0/,DR/NCON*0.0/
%       DATA SDI/NCON*0.0/,SDR/NCON*0.0/
%       DATA TAMP/
%      .  .632208, .294107, .121046, .079915, .023818,-.023589, .022994,
%      .  .019333,-.017871, .017192, .016018, .004671,-.004662,-.004519,
%      .  .004470, .004467, .002589,-.002455,-.002172, .001972, .001947,
%      .  .001914,-.001898, .001802, .001304, .001170, .001130, .001061,
%      . -.001022,-.001017, .001014, .000901,-.000857, .000855, .000855,
%      .  .000772, .000741, .000741,-.000721, .000698, .000658, .000654,
%      . -.000653, .000633, .000626,-.000598, .000590, .000544, .000479,
%      . -.000464, .000413,-.000390, .000373, .000366, .000366,-.000360,
%      . -.000355, .000354, .000329, .000328, .000319, .000302, .000279,
%      . -.000274,-.000272, .000248,-.000225, .000224,-.000223,-.000216,
%      .  .000211, .000209, .000194, .000185,-.000174,-.000171, .000159,
%      .  .000131, .000127, .000120, .000118, .000117, .000108, .000107,
%      .  .000105,-.000102, .000102, .000099,-.000096, .000095,-.000089,
%      . -.000085,-.000084,-.000081,-.000077,-.000072,-.000067, .000066,
%      .  .000064, .000063, .000063, .000063, .000062, .000062,-.000060,
%      .  .000056, .000053, .000051, .000050, .368645,-.262232,-.121995,
%      . -.050208, .050031,-.049470, .020620, .020613, .011279,-.009530,
%      . -.009469,-.008012, .007414,-.007300, .007227,-.007131,-.006644,
%      .  .005249, .004137, .004087, .003944, .003943, .003420, .003418,
%      .  .002885, .002884, .002160,-.001936, .001934,-.001798, .001690,
%      .  .001689, .001516, .001514,-.001511, .001383, .001372, .001371,
%      . -.001253,-.001075, .001020, .000901, .000865,-.000794, .000788,
%      .  .000782,-.000747,-.000745, .000670,-.000603,-.000597, .000542,
%      .  .000542,-.000541,-.000469,-.000440, .000438, .000422, .000410,
%      . -.000374,-.000365, .000345, .000335,-.000321,-.000319, .000307,
%      .  .000291, .000290,-.000289, .000286, .000275, .000271, .000263,
%      . -.000245, .000225, .000225, .000221,-.000202,-.000200,-.000199,
%      .  .000192, .000183, .000183, .000183,-.000170, .000169, .000168,
%      .  .000162, .000149,-.000147,-.000141, .000138, .000136, .000136,
%      .  .000127, .000127,-.000126,-.000121,-.000121, .000117,-.000116,
%      . -.000114,-.000114,-.000114, .000114, .000113, .000109, .000108,
%      .  .000106,-.000106,-.000106, .000105, .000104,-.000103,-.000100,
%      . -.000100,-.000100, .000099,-.000098, .000093, .000093, .000090,
%      . -.000088, .000083,-.000083,-.000082,-.000081,-.000079,-.000077,
%      . -.000075,-.000075,-.000075, .000071, .000071,-.000071, .000068,
%      .  .000068, .000065, .000065, .000064, .000064, .000064,-.000064,
%      . -.000060, .000056, .000056, .000053, .000053, .000053,-.000053,
%      .  .000053, .000053, .000052, .000050,-.066607,-.035184,-.030988,
%      .  .027929,-.027616,-.012753,-.006728,-.005837,-.005286,-.004921,
%      . -.002884,-.002583,-.002422, .002310, .002283,-.002037, .001883,
%      . -.001811,-.001687,-.001004,-.000925,-.000844, .000766, .000766,
%      . -.000700,-.000495,-.000492, .000491, .000483, .000437,-.000416,
%      . -.000384, .000374,-.000312,-.000288,-.000273, .000259, .000245,
%      . -.000232, .000229,-.000216, .000206,-.000204,-.000202, .000200,
%      .  .000195,-.000190, .000187, .000180,-.000179, .000170, .000153,
%      . -.000137,-.000119,-.000119,-.000112,-.000110,-.000110, .000107,
%      . -.000095,-.000095,-.000091,-.000090,-.000081,-.000079,-.000079,
%      .  .000077,-.000073, .000069,-.000067,-.000066, .000065, .000064,
%      . -.000062, .000060, .000059,-.000056, .000055,-.000051/
%       DATA IDD/
%      .  2, 0, 0, 0, 0, 0,   2, 2,-2, 0, 0, 0,   2,-1, 0, 1, 0, 0,  
%      .  2, 2, 0, 0, 0, 0,   2, 2, 0, 0, 1, 0,   2, 0, 0, 0,-1, 0,  
%      .  2,-1, 2,-1, 0, 0,   2,-2, 2, 0, 0, 0,   2, 1, 0,-1, 0, 0,  
%      .  2, 2,-3, 0, 0, 1,   2,-2, 0, 2, 0, 0,   2,-3, 2, 1, 0, 0,  
%      .  2, 1,-2, 1, 0, 0,   2,-1, 0, 1,-1, 0,   2, 3, 0,-1, 0, 0,  
%      .  2, 1, 0, 1, 0, 0,   2, 2, 0, 0, 2, 0,   2, 2,-1, 0, 0,-1,  
%      .  2, 0,-1, 0, 0, 1,   2, 1, 0, 1, 1, 0,   2, 3, 0,-1, 1, 0,  
%      .  2, 0, 1, 0, 0,-1,   2, 0,-2, 2, 0, 0,   2,-3, 0, 3, 0, 0,  
%      .  2,-2, 3, 0, 0,-1,   2, 4, 0, 0, 0, 0,   2,-1, 1, 1, 0,-1,  
%      .  2,-1, 3,-1, 0,-1,   2, 2, 0, 0,-1, 0,   2,-1,-1, 1, 0, 1,  
%      .  2, 4, 0, 0, 1, 0,   2,-3, 4,-1, 0, 0,   2,-1, 2,-1,-1, 0,  
%      .  2, 3,-2, 1, 0, 0,   2, 1, 2,-1, 0, 0,   2,-4, 2, 2, 0, 0,  
%      .  2, 4,-2, 0, 0, 0,   2, 0, 2, 0, 0, 0,   2,-2, 2, 0,-1, 0,  
%      .  2, 2,-4, 0, 0, 2,   2, 2,-2, 0,-1, 0,   2, 1, 0,-1,-1, 0,  
%      .  2,-1, 1, 0, 0, 0,   2, 2,-1, 0, 0, 1,   2, 2, 1, 0, 0,-1,  
%      .  2,-2, 0, 2,-1, 0,   2,-2, 4,-2, 0, 0,   2, 2, 2, 0, 0, 0,  
%      .  2,-4, 4, 0, 0, 0,   2,-1, 0,-1,-2, 0,   2, 1, 2,-1, 1, 0,  
%      .  2,-1,-2, 3, 0, 0,   2, 3,-2, 1, 1, 0,   2, 4, 0,-2, 0, 0,  
%      .  2, 0, 0, 2, 0, 0,   2, 0, 2,-2, 0, 0,   2, 0, 2, 0, 1, 0,  
%      .  2,-3, 3, 1, 0,-1,   2, 0, 0, 0,-2, 0,   2, 4, 0, 0, 2, 0,  
%      .  2, 4,-2, 0, 1, 0,   2, 0, 0, 0, 0, 2,   2, 1, 0, 1, 2, 0,  
%      .  2, 0,-2, 0,-2, 0,   2,-2, 1, 0, 0, 1,   2,-2, 1, 2, 0,-1,  
%      .  2,-1, 1,-1, 0, 1,   2, 5, 0,-1, 0, 0,   2, 1,-3, 1, 0, 1,  
%      .  2,-2,-1, 2, 0, 1,   2, 3, 0,-1, 2, 0,   2, 1,-2, 1,-1, 0,  
%      .  2, 5, 0,-1, 1, 0,   2,-4, 0, 4, 0, 0,   2,-3, 2, 1,-1, 0,  
%      .  2,-2, 1, 1, 0, 0,   2, 4, 0,-2, 1, 0,   2, 0, 0, 2, 1, 0,  
%      .  2,-5, 4, 1, 0, 0,   2, 0, 2, 0, 2, 0,   2,-1, 2, 1, 0, 0,  
%      .  2, 5,-2,-1, 0, 0,   2, 1,-1, 0, 0, 0,   2, 2,-2, 0, 0, 2,  
%      .  2,-5, 2, 3, 0, 0,   2,-1,-2, 1,-2, 0,   2,-3, 5,-1, 0,-1,  
%      .  2,-1, 0, 0, 0, 1,   2,-2, 0, 0,-2, 0,   2, 0,-1, 1, 0, 0,  
%      .  2,-3, 1, 1, 0, 1,   2, 3, 0,-1,-1, 0,   2, 1, 0, 1,-1, 0,  
%      .  2,-1, 2, 1, 1, 0,   2, 0,-3, 2, 0, 1,   2, 1,-1,-1, 0, 1,  
%      .  2,-3, 0, 3,-1, 0,   2, 0,-2, 2,-1, 0,   2,-4, 3, 2, 0,-1,  
%      .  2,-1, 0, 1,-2, 0,   2, 5, 0,-1, 2, 0,   2,-4, 5, 0, 0,-1,  
%      .  2,-2, 4, 0, 0,-2,   2,-1, 0, 1, 0, 2,   2,-2,-2, 4, 0, 0,  
%      .  2, 3,-2,-1,-1, 0,   2,-2, 5,-2, 0,-1,   2, 0,-1, 0,-1, 1,  
%      .  2, 5,-2,-1, 1, 0,   1, 1, 0, 0, 0, 0,   1,-1, 0, 0, 0, 0,  
%      .  1, 1,-2, 0, 0, 0,   1,-2, 0, 1, 0, 0,   1, 1, 0, 0, 1, 0,  
%      .  1,-1, 0, 0,-1, 0,   1, 2, 0,-1, 0, 0,   1, 0, 0, 1, 0, 0,  
%      .  1, 3, 0, 0, 0, 0,   1,-2, 2,-1, 0, 0,   1,-2, 0, 1,-1, 0,  
%      .  1,-3, 2, 0, 0, 0,   1, 0, 0,-1, 0, 0,   1, 1, 0, 0,-1, 0,  
%      .  1, 3, 0, 0, 1, 0,   1, 1,-3, 0, 0, 1,   1,-3, 0, 2, 0, 0,  
%      .  1, 1, 2, 0, 0, 0,   1, 0, 0, 1, 1, 0,   1, 2, 0,-1, 1, 0,  
%      .  1, 0, 2,-1, 0, 0,   1, 2,-2, 1, 0, 0,   1, 3,-2, 0, 0, 0,  
%      .  1,-1, 2, 0, 0, 0,   1, 1, 1, 0, 0,-1,   1, 1,-1, 0, 0, 1,  
%      .  1, 4, 0,-1, 0, 0,   1,-4, 2, 1, 0, 0,   1, 0,-2, 1, 0, 0,  
%      .  1,-2, 2,-1,-1, 0,   1, 3, 0,-2, 0, 0,   1,-1, 0, 2, 0, 0,  
%      .  1,-1, 0, 0,-2, 0,   1, 3, 0, 0, 2, 0,   1,-3, 2, 0,-1, 0,  
%      .  1, 4, 0,-1, 1, 0,   1, 0, 0,-1,-1, 0,   1, 1,-2, 0,-1, 0,  
%      .  1,-3, 0, 2,-1, 0,   1, 1, 0, 0, 2, 0,   1, 1,-1, 0, 0,-1,  
%      .  1,-1,-1, 0, 0, 1,   1, 0, 2,-1, 1, 0,   1,-1, 1, 0, 0,-1,  
%      .  1,-1,-2, 2, 0, 0,   1, 2,-2, 1, 1, 0,   1,-4, 0, 3, 0, 0,  
%      .  1,-1, 2, 0, 1, 0,   1, 3,-2, 0, 1, 0,   1, 2, 0,-1,-1, 0,  
%      .  1, 0, 0, 1,-1, 0,   1,-2, 2, 1, 0, 0,   1, 4,-2,-1, 0, 0,  
%      .  1,-3, 3, 0, 0,-1,   1,-2, 1, 1, 0,-1,   1,-2, 3,-1, 0,-1,  
%      .  1, 0,-2, 1,-1, 0,   1,-2,-1, 1, 0, 1,   1, 4,-2, 1, 0, 0,  
%      .  1,-4, 4,-1, 0, 0,   1,-4, 2, 1,-1, 0,   1, 5,-2, 0, 0, 0,  
%      .  1, 3, 0,-2, 1, 0,   1,-5, 2, 2, 0, 0,   1, 2, 0, 1, 0, 0,  
%      .  1, 1, 3, 0, 0,-1,   1,-2, 0, 1,-2, 0,   1, 4, 0,-1, 2, 0,  
%      .  1, 1,-4, 0, 0, 2,   1, 5, 0,-2, 0, 0,   1,-1, 0, 2, 1, 0,  
%      .  1,-2, 1, 0, 0, 0,   1, 4,-2, 1, 1, 0,   1,-3, 4,-2, 0, 0,  
%      .  1,-1, 3, 0, 0,-1,   1, 3,-3, 0, 0, 1,   1, 5,-2, 0, 1, 0,  
%      .  1, 1, 2, 0, 1, 0,   1, 2, 0, 1, 1, 0,   1,-5, 4, 0, 0, 0,  
%      .  1,-2, 0,-1,-2, 0,   1, 5, 0,-2, 1, 0,   1, 1, 2,-2, 0, 0,  
%      .  1, 1,-2, 2, 0, 0,   1,-2, 2, 1, 1, 0,   1, 0, 3,-1, 0,-1,  
%      .  1, 2,-3, 1, 0, 1,   1,-2,-2, 3, 0, 0,   1,-1, 2,-2, 0, 0,  
%      .  1,-4, 3, 1, 0,-1,   1,-4, 0, 3,-1, 0,   1,-1,-2, 2,-1, 0,  
%      .  1,-2, 0, 3, 0, 0,   1, 4, 0,-3, 0, 0,   1, 0, 1, 1, 0,-1,  
%      .  1, 2,-1,-1, 0, 1,   1, 2,-2, 1,-1, 0,   1, 0, 0,-1,-2, 0,  
%      .  1, 2, 0, 1, 2, 0,   1, 2,-2,-1,-1, 0,   1, 0, 0, 1, 2, 0,  
%      .  1, 0, 1, 0, 0, 0,   1, 2,-1, 0, 0, 0,   1, 0, 2,-1,-1, 0,  
%      .  1,-1,-2, 0,-2, 0,   1,-3, 1, 0, 0, 1,   1, 3,-2, 0,-1, 0,  
%      .  1,-1,-1, 0,-1, 1,   1, 4,-2,-1, 1, 0,   1, 2, 1,-1, 0,-1,  
%      .  1, 0,-1, 1, 0, 1,   1,-2, 4,-1, 0, 0,   1, 4,-4, 1, 0, 0,  
%      .  1,-3, 1, 2, 0,-1,   1,-3, 3, 0,-1,-1,   1, 1, 2, 0, 2, 0,  
%      .  1, 1,-2, 0,-2, 0,   1, 3, 0, 0, 3, 0,   1,-1, 2, 0,-1, 0,  
%      .  1,-2, 1,-1, 0, 1,   1, 0,-3, 1, 0, 1,   1,-3,-1, 2, 0, 1,  
%      .  1, 2, 0,-1, 2, 0,   1, 6,-2,-1, 0, 0,   1, 2, 2,-1, 0, 0,  
%      .  1,-1, 1, 0,-1,-1,   1,-2, 3,-1,-1,-1,   1,-1, 0, 0, 0, 2,  
%      .  1,-5, 0, 4, 0, 0,   1, 1, 0, 0, 0,-2,   1,-2, 1, 1,-1,-1,  
%      .  1, 1,-1, 0, 1, 1,   1, 1, 2, 0, 0,-2,   1,-3, 1, 1, 0, 0,  
%      .  1,-4, 4,-1,-1, 0,   1, 1, 0,-2,-1, 0,   1,-2,-1, 1,-1, 1,  
%      .  1,-3, 2, 2, 0, 0,   1, 5,-2,-2, 0, 0,   1, 3,-4, 2, 0, 0,  
%      .  1, 1,-2, 0, 0, 2,   1,-1, 4,-2, 0, 0,   1, 2, 2,-1, 1, 0,  
%      .  1,-5, 2, 2,-1, 0,   1, 1,-3, 0,-1, 1,   1, 1, 1, 0, 1,-1,  
%      .  1, 6,-2,-1, 1, 0,   1,-2, 2,-1,-2, 0,   1, 4,-2, 1, 2, 0,  
%      .  1,-6, 4, 1, 0, 0,   1, 5,-4, 0, 0, 0,   1,-3, 4, 0, 0, 0,  
%      .  1, 1, 2,-2, 1, 0,   1,-2, 1, 0,-1, 0,   0, 2, 0, 0, 0, 0,  
%      .  0, 1, 0,-1, 0, 0,   0, 0, 2, 0, 0, 0,   0, 0, 0, 0, 1, 0,  
%      .  0, 2, 0, 0, 1, 0,   0, 3, 0,-1, 0, 0,   0, 1,-2, 1, 0, 0,  
%      .  0, 2,-2, 0, 0, 0,   0, 3, 0,-1, 1, 0,   0, 0, 1, 0, 0,-1,  
%      .  0, 2, 0,-2, 0, 0,   0, 2, 0, 0, 2, 0,   0, 3,-2, 1, 0, 0,  
%      .  0, 1, 0,-1,-1, 0,   0, 1, 0,-1, 1, 0,   0, 4,-2, 0, 0, 0,  
%      .  0, 1, 0, 1, 0, 0,   0, 0, 3, 0, 0,-1,   0, 4, 0,-2, 0, 0,  
%      .  0, 3,-2, 1, 1, 0,   0, 3,-2,-1, 0, 0,   0, 4,-2, 0, 1, 0,  
%      .  0, 0, 2, 0, 1, 0,   0, 1, 0, 1, 1, 0,   0, 4, 0,-2, 1, 0,  
%      .  0, 3, 0,-1, 2, 0,   0, 5,-2,-1, 0, 0,   0, 1, 2,-1, 0, 0,  
%      .  0, 1,-2, 1,-1, 0,   0, 1,-2, 1, 1, 0,   0, 2,-2, 0,-1, 0,  
%      .  0, 2,-3, 0, 0, 1,   0, 2,-2, 0, 1, 0,   0, 0, 2,-2, 0, 0,  
%      .  0, 1,-3, 1, 0, 1,   0, 0, 0, 0, 2, 0,   0, 0, 1, 0, 0, 1,  
%      .  0, 1, 2,-1, 1, 0,   0, 3, 0,-3, 0, 0,   0, 2, 1, 0, 0,-1,  
%      .  0, 1,-1,-1, 0, 1,   0, 1, 0, 1, 2, 0,   0, 5,-2,-1, 1, 0,  
%      .  0, 2,-1, 0, 0, 1,   0, 2, 2,-2, 0, 0,   0, 1,-1, 0, 0, 0,  
%      .  0, 5, 0,-3, 0, 0,   0, 2, 0,-2, 1, 0,   0, 1, 1,-1, 0,-1,  
%      .  0, 3,-4, 1, 0, 0,   0, 0, 2, 0, 2, 0,   0, 2, 0,-2,-1, 0,  
%      .  0, 4,-3, 0, 0, 1,   0, 3,-1,-1, 0, 1,   0, 0, 2, 0, 0,-2,  
%      .  0, 3,-3, 1, 0, 1,   0, 2,-4, 2, 0, 0,   0, 4,-2,-2, 0, 0,  
%      .  0, 3, 1,-1, 0,-1,   0, 5,-4, 1, 0, 0,   0, 3,-2,-1,-1, 0,  
%      .  0, 3,-2, 1, 2, 0,   0, 4,-4, 0, 0, 0,   0, 6,-2,-2, 0, 0,  
%      .  0, 5, 0,-3, 1, 0,   0, 4,-2, 0, 2, 0,   0, 2, 2,-2, 1, 0,  
%      .  0, 0, 4, 0, 0,-2,   0, 3,-1, 0, 0, 0,   0, 3,-3,-1, 0, 1,  
%      .  0, 4, 0,-2, 2, 0,   0, 1,-2,-1,-1, 0,   0, 2,-1, 0, 0,-1,  
%      .  0, 4,-4, 2, 0, 0,   0, 2, 1, 0, 1,-1,   0, 3,-2,-1, 1, 0,  
%      .  0, 4,-3, 0, 1, 1,   0, 2, 0, 0, 3, 0,   0, 6,-4, 0, 0, 0/
% 
% *  Initialize variables.
%       K   = 0
%       NLP = 0
%       NDI = 0
%       NSD = 0
% 
%       DO LL=1,NIN
% *  See if Doodson numbers match
%          DO KK=1,NT
%             II = 0
%             DO I=1,6
%                II = II + IABS(IDD(I,KK)-IDTIN(I,LL))
%             ENDDO
%             IF(II.EQ.0) GO TO 5
%          ENDDO
% *  If you have a match, put line into array
%  5       IF(II.EQ.0.AND.K.LT.NCON) THEN
%             K = K + 1
%             RL(K) = AMPIN(LL)*COS(DTR*PHIN(LL))/ABS(TAMP(KK))
%             AIM(K)= AMPIN(LL)*SIN(DTR*PHIN(LL))/ABS(TAMP(KK))
% *+---------------------------------------------------------------------
% *  Now have real and imaginary parts of admittance, scaled by Cartwright-
% *  Edden amplitude. Admittance phase is whatever was used in the original
% *  expression. (Usually phase is given relative to some reference,
% *  but amplitude is in absolute units). Next get frequency.
% *----------------------------------------------------------------------
%             CALL TDFRPH(IDD(1,KK),FR,PR)
%             RF(K) = FR
%          ENDIF
%       ENDDO
% *+---------------------------------------------------------------------
% *  Done going through constituents; there are k of them.
% *  Have specified admittance at a number of points. Sort these by frequency
% *  and separate diurnal and semidiurnal, recopying admittances to get them
% *  in order using Shell Sort.
% *----------------------------------------------------------------------
% 
%       CALL SHELLS(RF,KEY,K)
%       DO I=1,K
%          IF(RF(I).LT.0.5) NLP = NLP + 1
%          IF(RF(I).LT.1.5.AND.RF(I).GT.0.5) NDI = NDI + 1
%          IF(RF(I).LT.2.5.AND.RF(I).GT.1.5) NSD = NSD + 1
%          SCR(I) = RL(KEY(I))
%       ENDDO
%       DO I=1,K
%          RL(I) = SCR(I)
%          SCR(I) = AIM(KEY(I))
%       ENDDO
%       DO I=1,K
%          AIM(I) = SCR(I)
%       ENDDO
% *+---------------------------------------------------------------------
% *  now set up splines (8 cases - four species, each real and imaginary)
% *  We have to allow for the case when there are no constituent amplitudes
% *  for the long-period tides.
% *----------------------------------------------------------------------
%       IF(NLP.NE.0) CALL SPLINE(NLP,RF,RL,ZDR,SCR)
%       IF(NLP.NE.0) CALL SPLINE(NLP,RF,AIM,ZDI,SCR)
%       CALL SPLINE(NDI,RF(NLP+1),RL(NLP+1),DR,SCR)
%       CALL SPLINE(NDI,RF(NLP+1),AIM(NLP+1),DI,SCR)
%       CALL SPLINE(NSD,RF(NLP+NDI+1),RL(NLP+NDI+1),SDR,SCR)
%       CALL SPLINE(NSD,RF(NLP+NDI+1),AIM(NLP+NDI+1),SDI,SCR)
% *  Evaluate all harmonics using the interpolated admittance
%       J = 1
%       DO I=1,NT
%          IF(IDD(1,I).EQ.0.AND.NLP.EQ.0) GO TO 11
%          CALL TDFRPH(IDD(1,I),F(J),P(J))
% *  Compute phase corrections to equilibrium tide using function EVAL
%          IF(IDD(1,I).EQ.0) P(J) = P(J) + 180.
%          IF(IDD(1,I).EQ.1) P(J) = P(J) + 90.
%          SF = F(J)
%          IF(IDD(1,I).EQ.0) RE = EVAL(SF,NLP,RF,RL,ZDR)
%          IF(IDD(1,I).EQ.0) AM = EVAL(SF,NLP,RF,AIM,ZDI)
%          IF(IDD(1,I).EQ.1) RE = EVAL(SF,NDI,RF(NLP+1),RL(NLP+1),DR)
%          IF(IDD(1,I).EQ.1) AM = EVAL(SF,NDI,RF(NLP+1),AIM(NLP+1),DI)
%          IF(IDD(1,I).EQ.2) RE =
%      .      EVAL(SF,NSD,RF(NLP+NDI+1),RL(NLP+NDI+1),SDR)
%          IF(IDD(1,I).EQ.2) AM =
%      .      EVAL(SF,NSD,RF(NLP+NDI+1),AIM(NLP+NDI+1),SDI)
%          AMP(J) = TAMP(I)*SQRT(RE**2+AM**2)
%          P(J) = P(J) + ATAN2(AM,RE)/DTR
%          IF(P(J).GT.180) P(J)=P(J)-360.
%          J = J + 1
%  11      CONTINUE
%       ENDDO
%       NOUT = J - 1
%       RETURN
% 
% *  Finished.
% 
% *+----------------------------------------------------------------------
% *
% *  Copyright (C) 2008
% *  IERS Conventions Center
% *
% *  ==================================
% *  IERS Conventions Software License
% *  ==================================
% *
% *  NOTICE TO USER:
% *
% *  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
% *  WHICH APPLY TO ITS USE.
% *
% *  1. The Software is provided by the IERS Conventions Center ("the
% *     Center").
% *
% *  2. Permission is granted to anyone to use the Software for any
% *     purpose, including commercial applications, free of charge,
% *     subject to the conditions and restrictions listed below.
% *
% *  3. You (the user) may adapt the Software and its algorithms for your
% *     own purposes and you may distribute the resulting "derived work"
% *     to others, provided that the derived work complies with the
% *     following requirements:
% *
% *     a) Your work shall be clearly identified so that it cannot be
% *        mistaken for IERS Conventions software and that it has been
% *        neither distributed by nor endorsed by the Center.
% *
% *     b) Your work (including source code) must contain descriptions of
% *        how the derived work is based upon and/or differs from the
% *        original Software.
% *
% *     c) The name(s) of all modified routine(s) that you distribute
% *        shall be changed.
% * 
% *     d) The origin of the IERS Conventions components of your derived
% *        work must not be misrepresented; you must not claim that you
% *        wrote the original Software.
% *
% *     e) The source code must be included for all routine(s) that you
% *        distribute.  This notice must be reproduced intact in any
% *        source distribution. 
% *
% *  4. In any published work produced by the user and which includes
% *     results achieved by using the Software, you shall acknowledge
% *     that the Software was used in obtaining those results.
% *
% *  5. The Software is provided to the user "as is" and the Center makes
% *     no warranty as to its use or performance.   The Center does not
% *     and cannot warrant the performance or results which the user may
% *     obtain by using the Software.  The Center makes no warranties,
% *     express or implied, as to non-infringement of third party rights,
% *     merchantability, or fitness for any particular purpose.  In no
% *     event will the Center be liable to the user for any consequential,
% *     incidental, or special damages, including any lost profits or lost
% *     savings, even if a Center representative has been advised of such
% *     damages, or for any claim by any third party.
% *
% *  Correspondence concerning IERS Conventions software should be
% *  addressed as follows:
% *
% *                     Gerard Petit
% *     Internet email: gpetit[at]bipm.org
% *     Postal address: IERS Conventions Center
% *                     Time, frequency and gravimetry section, BIPM
% *                     Pavillon de Breteuil
% *                     92312 Sevres  FRANCE
% *
% *     or
% *
% *                     Brian Luzum
% *     Internet email: brian.luzum[at]usno.navy.mil
% *     Postal address: IERS Conventions Center
% *                     Earth Orientation Department
% *                     3450 Massachusetts Ave, NW
% *                     Washington, DC 20392
% *
% *
% *-----------------------------------------------------------------------
%       END
%%

% Input
%   RF      frequency of the main 11 tides
%   RL, AIM real and imaginary parts of admittance, scaled by Cartwright-
%           Edden amplitude of the main tides
%   F       frequencies of tidal constituents (342)
%   P       phases of tidal constituents (342)
%   TAMP    Cartwright-Edden amplitudes (342) 
%   IDD1    first digit of Doodson number (342)


% Output:
%   AMP        Amplitude due to ocean loading
%   P          Corrected phase (implementation of 342 tides) due to ocean loading [deg]

%   Coded for VieVS: 
%   15 Jul 2011 by Hana Spicakova


 function [AMP,P] = libiers_admint(RF,RL,AIM,F,P,TAMP,IDD1)

% *+----------------------------------------------------------------------
% *  The parameters NT set the number of harmonics used in the prediction
% *  and the number of constituents whose amp and phase may be specified (ncon)
% *-----------------------------------------------------------------------
DTR=pi()/180;
SCR=AIM;

% *+---------------------------------------------------------------------
% *  now set up splines (8 cases - four species, each real and imaginary)
% *  We have to allow for the case when there are no constituent amplitudes
% *  for the long-period tides.
% *----------------------------------------------------------------------
NLP = 0;
NDI = 0;
NSD = 0;
NIN=size(RF,2);
for I = 1:NIN % K - length(RF)
    if (RF(I) < 0.5); NLP = NLP + 1; end % long periods
    if (RF(I) < 1.5 && RF(I) > 0.5); NDI = NDI + 1; end % diurnal
    if (RF(I) < 2.5 && RF(I) > 1.5); NSD = NSD + 1; end % semidiurnal
end
% NLP=3; NDI=4; NSD=4;

if (NLP ~= 0)
    ZDR = libiers_spline(NLP,RF,RL,SCR);
    ZDI = libiers_spline(NLP,RF,AIM,SCR);
end

DR = libiers_spline(NDI,RF(NLP+1:end),RL(NLP+1:end),SCR);
DI = libiers_spline(NDI,RF(NLP+1:end),AIM(NLP+1:end),SCR);
SDR = libiers_spline(NSD,RF(NLP+NDI+1:end),RL(NLP+NDI+1:end),SCR);
SDI = libiers_spline(NSD,RF(NLP+NDI+1:end),AIM(NLP+NDI+1:end),SCR);


% *  Evaluate all harmonics using the interpolated admittance
% *  Compute phase corrections to equilibrium tide using function EVAL
SF=F(IDD1==0);
RE_lp = libiers_eval(SF,NLP,RF,RL,ZDR);
AM_lp = libiers_eval(SF,NLP,RF,AIM,ZDI);

clear SF    
SF=F(IDD1==1);
RE_diu = libiers_eval(SF,NDI,RF(NLP+1:end),RL(NLP+1:end),DR);      
AM_diu = libiers_eval(SF,NDI,RF(NLP+1:end),AIM(NLP+1:end),DI);

clear SF
SF=F(IDD1==2);
RE_smd = libiers_eval(SF,NSD,RF(NLP+NDI+1:end),RL(NLP+NDI+1:end),SDR);
AM_smd = libiers_eval(SF,NSD,RF(NLP+NDI+1:end),AIM(NLP+NDI+1:end),SDI);


AMP = [TAMP(IDD1==2).*sqrt(RE_smd.^2 + AM_smd.^2)' ...
    TAMP(IDD1==1).*sqrt(RE_diu.^2 + AM_diu.^2)' ...
    TAMP(IDD1==0).*sqrt(RE_lp.^2  + AM_lp.^2)'] ;

dP =  [ atan2(AM_smd,RE_smd)./DTR; atan2(AM_diu,RE_diu)./DTR; atan2(AM_lp,RE_lp)./DTR];

NT = length(TAMP);
for I =1:NT %342
     if (IDD1(I) == 0); P(I) = P(I) + 180; end
     if (IDD1(I) == 1); P(I) = P(I) + 90; end
end

 
P=P+dP';

id=find(P > 180);
P(id)=P(id)-360;



