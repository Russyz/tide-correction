function [TDB1, TDB2, TT1,TT2]  = UTC2TDB(YR,MONTH,DATE,hour,min,sec);
% TT-TDB
fold = 'D:\谷歌下载\inpop19a_TDB_m100_p100_asc';
filename = dir([fold,'\','*pos_TT*']);
startRow = 3;
formatSpec = '%23f%23f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%[^\n\r]';
fileID = fopen([fold ,'\', filename.name],'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
inpop19aTDBm100p100ascposTT = [dataArray{1:end-1}];
clearvars fold filename startRow formatSpec fileID dataArray ans;
% YR = 2006;
% MONTH = 12 ;
% DATE = 20;
% hour = 0;
% min = 0;
% sec = 0;
% jd = juliandate(YR,MONTH,DATE,hour,min,sec);
UTC1 = juliandate(YR,MONTH,DATE);
UTC2 = hour/24 + min/24*60 + sec/86400;

[TAI1,TAI2] = iauUtctai(UTC1,UTC2);
[TT1, TT2] = iauTaitt(TAI1,TAI2);
TT = TT1 + TT2;
%dtr            TDB-TT in seconds
cpm = find(TT>inpop19aTDBm100p100ascposTT(:,1)&TT<inpop19aTDBm100p100ascposTT(:,2));
N = length(inpop19aTDBm100p100ascposTT(1,:))-2;
Ta = (inpop19aTDBm100p100ascposTT(cpm,1));
Tb = (inpop19aTDBm100p100ascposTT(cpm,2));
C = (inpop19aTDBm100p100ascposTT(cpm,3:end));
% TBD - TT in seconds
dtr = Cheb1D(TT, N, Ta, Tb, C);
[TDB1, TDB2] = iauTttdb(TT1, TT2, dtr);
% TDB = TDB1 + TDB2;
end
