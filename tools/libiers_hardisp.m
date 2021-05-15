function [cto_rsw] = libiers_hardisp(mjd,leap,cto,F_lsr,P_lsr, TAMP_lsr, IDD1_lsr)

%  (Cartwright-Tayler) Doodson numbers of tides used in Scherneck lists:
%      M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa

IDT = [2, 0, 0, 0, 0, 0;
    2, 2,-2, 0, 0, 0;
    2,-1, 0, 1, 0, 0;
    2, 2, 0, 0, 0, 0;
    1, 1, 0, 0, 0, 0;
    1,-1, 0, 0, 0, 0;
    1, 1,-2, 0, 0, 0;
    1,-2, 0, 1, 0, 0;
    0, 2, 0, 0, 0, 0;
    0, 1, 0,-1, 0, 0;
    0, 0, 2, 0, 0, 0];

% Cartwright-Edden amplitudes of tides used in Scherneck lists (extracted from function libiers_tdfrph_call)
TAMPT = [ 0.6322    0.2941    0.1210    0.0799    0.3686   -0.2622   -0.1220   -0.0502   -0.0666   -0.0352   -0.0310];

%---------------------------------------------------------------------
%  Read in amplitudes and phases, in standard "Scherneck" form, from
%  standard input
%----------------------------------------------------------------------

TAMP = cto(1:3,1:end);
TPH = cto(4:6,1:end);

% Change sign for phase, to be negative for lags
TPH = -TPH;


% Get real and imaginary parts of admittance, scaled by Cartwright-
% Edden amplitude and frequency
[RF RL AIM] = libiers_admint_part1(TAMP, TPH, TAMPT, IDT, mjd, leap);

% *+---------------------------------------------------------------------
% *
% *  Find amplitudes and phases for all constituents, for each of the
% *  three displacements.
% *
% *  BLQ format order is vertical, horizontal EW, horizontal NS
% *
% *----------------------------------------------------------------------
% vertical
[AZ PZ] = libiers_admint (RF,RL(1,:),AIM(1,:),F_lsr,P_lsr, TAMP_lsr,IDD1_lsr);
% West
[AW PW] = libiers_admint (RF,RL(2,:),AIM(2,:),F_lsr,P_lsr, TAMP_lsr,IDD1_lsr);
% South
[AS PS] = libiers_admint (RF,RL(3,:),AIM(3,:),F_lsr,P_lsr, TAMP_lsr,IDD1_lsr);

DR=pi()/180;
PZ = DR.*PZ;
PS = DR.*PS;
PW = DR.*PW;

% We don't compute the displacement recursively, we want only the
% displacement at the time of the scan.

DZ = sum(AZ.*cos(PZ));
DS = sum(AS.*cos(PS));
DW = sum(AW.*cos(PW));


cto_rsw = [DZ DS DW];