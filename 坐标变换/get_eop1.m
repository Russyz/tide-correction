function [MJD_UTC, TAI_UTC, UT1_UTC, dx, dy] = get_eop1(YR,MONTH,DATE,hour,min,sec)
MJD_UTC = juliandate(YR,MONTH,DATE,hour,min,sec) - 2400000.5;
fullURL = 'ftp://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now';
filename = 'eop_data.txt';
if ~exist(filename)
    filename = urlwrite(fullURL,filename);
end
startRow = 13;
formatSpec = '%4f%4f%4f%7f%11f%11f%12f%12f%11f%11f%11f%11f%11f%11f%12f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);
fclose(fileID);
EOP_DATA = [dataArray{1:end-1}];
% delete(filename);
clearvars filename startRow formatSpec fileID dataArray ans fullURL;

% headers = ['Year', 'Month', 'Day', 'MJD', 'x', 'y', 'UT1-UTC', 'LOD', 'dX', 'dY', 'x Err',
%     'y Err', 'UT1-UTC Err', 'LOD Err', 'dX Err', 'dY Err']

MJD = EOP_DATA(:,4);
x = EOP_DATA(:,5);
y = EOP_DATA(:,6);
ut1_utc = EOP_DATA(:,7);
dx=lagint(MJD,x,MJD_UTC,4);
dy=lagint(MJD,y,MJD_UTC,4);
UT1_UTC=lagint(MJD,ut1_utc,MJD_UTC,4);
addpath(genpath('E:\潮汐修正\test file\坐标变换\'))
format long
leapsec_data = importdata('CDFLeapSeconds.txt');
year = leapsec_data(:,1);
mon = leapsec_data(:,2);
dat = leapsec_data(:,3);
leap_mjd = juliandate(year,mon,dat) - 2400000.5;
% for i = 1:(length(year)-1)
%     if (year(i) <YR & YR<year(i+1))
%         TAI_UTC=leapsec_data(i,4);
%     else
%         TAI_UTC=leapsec_data(end,4);
%     end
% end
w = find(MJD_UTC>leap_mjd);
TAI_UTC = leapsec_data(max(w),4);
end
