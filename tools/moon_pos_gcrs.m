function moon_position_gcrs = moon_pos_gcrs(YR,MONTH,DATE,hour,min,sec)
% Input time in utc to get moon position in GCRS in kilometers.

fold = 'D:\谷歌下载\inpop19a_TDB_m100_p100_asc';
filename = dir([fold,'\','*pos_Moo*']);
startRow = 3;
formatSpec = '%23f%23f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%f%[^\n\r]';
fileID = fopen([fold ,'\', filename.name],'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
inpop19aTDBm100p100ascposMoo = [dataArray{1:end-1}];
clearvars fold filename startRow formatSpec fileID dataArray ans;

[TDB1, TDB2, TT1,TT2]  = UTC2TDB(YR,MONTH,DATE,hour,min,sec);
TT = TT1 + TT2;
cpm = find(TT>inpop19aTDBm100p100ascposMoo(:,1)&TT<inpop19aTDBm100p100ascposMoo(:,2));
N = length(inpop19aTDBm100p100ascposMoo(1,:))-2;
Ta = (inpop19aTDBm100p100ascposMoo(cpm(1),1));
Tb = (inpop19aTDBm100p100ascposMoo(cpm(1),2));
Cx = (inpop19aTDBm100p100ascposMoo(cpm(1),3:end));
Cy = (inpop19aTDBm100p100ascposMoo(cpm(2),3:end));
Cz = (inpop19aTDBm100p100ascposMoo(cpm(3),3:end));
TDB = TDB1 + TDB2;
moon_position_gcrs = Cheb3D(TDB, N, Ta, Tb, Cx, Cy, Cz);
end