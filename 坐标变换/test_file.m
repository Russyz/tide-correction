% YR = 2021;
% MONTH=3;
% DATE=17;
% hour=6;
% min=0;
% sec =0;
YR = 2004;
MONTH=4;
DATE=6;
hour = 7; min = 51; sec = 28.386009;
solid_tide_dehan(YR,MONTH,DATE,hour,min,sec)
% 
% EOP_data = get_eop(YR,MONTH,DATE,hour,min,sec);
% 
% addpath(genpath('E:\潮汐修正\test file')) %路径
% dxtide = solid_tide_dehan(YR,MONTH,DATE,hour,min,sec)
% dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour,min,sec)
% dxatide = atm_tide(YR,MONTH,DATE,hour,min,sec)
% dxptide = solid_ptide(YR,MONTH,DATE,hour,min,sec)
% dxoptide = ocean_ptide(YR,MONTH,DATE,hour,min,sec)
% addpath(genpath('E:\潮汐修正\test file\hardisp')) %路径
% [F,P, TAMP, IDD1]=libiers_tdfrph_call(EOP_data(1),EOP_data(3));
% cto_data = importdata('inputOTL.txt');
% coordinate_station=[-2358701.936, 5410608.308, 2410085.530]
% addpath(genpath('E:\潮汐修正\IERS模型\rebuilt.VLBI-develop\VLBI-develop\CODE\VIE_MOD'))
% % [cto]=ctocean(EOP_data(1),EOP_data(3),coordinate_station,cto_data)
