function eop_data = get_eop(YR,MONTH,DATE,hour,min,sec)
% Input YR,MONTH,DATE,hour,min and sec in UTC time,
% then you will get EOP parameters such as
% MJD_UTC, TAI_UTC, UT1_UTC, dx (xpole), dy(ypole)
eopdata_name = dir('E:\潮汐修正\test file\坐标变换\aeroiersdata20*');

MJD_UTC = juliandate(YR,MONTH,DATE,hour,min,sec) - 2400000.5;

if length(eopdata_name) == 0
    aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt')
else
    load(eopdata_name.name);
    %     if mjd(length(dxy)) < MJD_UTC
    if mjd(length(dxy)) < MJD_UTC
        delete(eopdata_name.name)
        addpath(genpath('E:\潮汐修正\test file\坐标变换\'))
        aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt')
        % aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/EOP_C01_IAU2000_1900-now.txt')
    end
%     if length(find(char(pmip) == 'I')) < MJD_UTC
%         delete(eopdata_name.name)
%         addpath(genpath('E:\潮汐修正\test file\坐标变换\'))
%         aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt')
%         % aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/EOP_C01_IAU2000_1900-now.txt')
%     end
    
end

% MJD_UTC = juliandate(YR,MONTH,DATE,hour,min,sec) - 2400000.5;
% if mjd(length(dxy)) < MJD_UTC
%     delete(eopdata_name.name)
%     addpath(genpath('E:\潮汐修正\test file\坐标变换\'))
%     aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt')
%     % aeroReadIERSData(pwd,'url','https://datacenter.iers.org/data/latestVersion/EOP_C01_IAU2000_1900-now.txt')
% end
% MJD_UTC = juliandate(YR,MONTH,DATE,hour,min,sec);
dx=lagint(mjd(1:length(pm(:,1))),pm(:,1),MJD_UTC,4);
dy=lagint(mjd(1:length(pm(:,1))),pm(:,2),MJD_UTC,4);
UT1_UTC=lagint(mjd(1:length(ut1utc)),ut1utc,MJD_UTC,10);
addpath(genpath('E:\潮汐修正\test file\坐标变换\'))
leapsec_data = importdata('CDFLeapSeconds.txt');
year = leapsec_data(:,1);
for i = 1:(length(year)-1)
    if (year(i) <YR & YR<year(i+1))
        TAI_UTC=leapsec_data(i,4);
    else
        TAI_UTC=leapsec_data(end,4);
    end
end

addpath(genpath('E:\潮汐修正\test file\星历\JPL_DE430')) %路径
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_TT = MJD_UTC + TT_UTC/86400;
eop_data = [MJD_UTC, TAI_UTC, UT1_UTC, dx, dy];
end
