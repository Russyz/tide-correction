% clear all;clc
% YR = year(datetime);
% MONTH = month(datetime);
% DATE = day(datetime)-1;
YR = 2006;
MONTH = 12 ;
DATE = 20;
coordinate_station=[-2358691.210,5410611.484,2410087.607];
format long g
% Constant tidal displacements from Earth and Sun using h2 = 0.0476 and l2 = 0.0107.
% Array ∆X ∆Y ∆Z ∆R ∆East ∆North
% meters meters meters meters meters meters
% Apollo 11 0.489 0.048 0.001 0.467 -0.151 -0.004
% Apollo 14 0.538 -0.046 -0.010 0.526 0.118 0.024
% Apollo 15 0.459 0.006 0.044 0.431 -0.023 -0.163
% Lunokhod 1 0.203 0.044 -0.061 0.073 0.152 -0.135
% Lunokhod 2 0.315 -0.002 -0.002 0.242 -0.164 -0.120

% COORD of Reflectors on the surface of moon in LCRF PA MER

Apollo14 = [1652689.369 1652818.682
    520998.431 520454.587
    109729.869 110361.165];
Apollo15 = [1554678.104 1554937.504
    98094.498 98604.886
    765005.863 764412.810];
Lunokhod2 = [1339363.598 1339388.213
    801870.995 802310.527
    756359.260 755849.393];
Lunokhod1 = [1114291.452 1114958.865
    781299.273 780934.127
    1076059.049 1075632.692];
constants;
j = 2;

for hour = 0
    for min = 0
        for sec = 880:1:900
            % UTC 时间
            UTC1 = juliandate(YR,MONTH,DATE);
            UTC2 = hour/24 + min/24*60 + sec/86400;
            % 原子时
            [TAI1,TAI2] = iauUtctai(UTC1,UTC2);
            % 大地时
            [TT1, TT2] = iauTaitt(TAI1,TAI2);
            [YR, MONTH, DATE, fd]= iauJd2cal(TT1,TT2);
            
            i=4*j-3;
            k = 3*j-2;
            j = j+1;
            
            
            addpath(genpath('E:\潮汐修正\test file\坐标变换')) %路径
            [MJD_UTC, TAI_UTC, dt, dx, dy] = get_eop1(YR,MONTH,DATE,hour,min,sec);
            STAT_LONG=113.55421725;
            STAT_LAT=22.34639411111111;
            STAT_HEI=405.713;
            % PYTHON to get TT-TDB
            % peph = py.calcephpy.CalcephBin.open("d:/jpl/inpop19a_TDB_m100_p100_tt.dat");
            % EPH_TT_MINUS_TDB = peph.compute(TT1, TT2, 0, 16);
            
            [TDB1, TDB2, TT1,TT2]  = UTC2TDB(YR,MONTH,DATE,hour,min,sec);
            TDB = TDB1 + TDB2;
            [TCG1, TCG2] = iauTttcg(TT1,TT2);
            TCG = TCG1 + TCG2;
            moon_pos = planetEphemeris(TCG,'Earth','Moon',430);
            % moon euler angle
            moon_euler = moonLibration(TDB,430);
            
            
            addpath('E:\潮汐修正\test file')
            rot = rotation(3, -moon_euler(1))*rotation(1, -moon_euler(2))*rotation(3, -moon_euler(3));
            
            format long g
            BCRS = moon_pos' * 1e3;
            
            BCRS_a15 = moon_pos' * 1e3 + rot * Apollo15(:,1);
            moon_pos1=  moon_pos_gcrs(YR,MONTH,DATE,hour,min,sec);
            moon_eul_gcrs(YR,MONTH,DATE,hour,min,sec)
            
            %             GC2IT1 = GCRS2ITRS1(YR,MONTH,DATE,hour,min,sec);
            GC2IT2 = GCRS2ITRS2(YR,MONTH,DATE,hour,min,sec);
            
            cpfpre = importdata('C:\Users\39736\Downloads\cpf_2.00a\cpf_sample_code_v2.00a\cpf_llr_c\luncenter_cpf_061220_35501.utx',' ',3,1);
            pre_data = cpfpre.data;
            
            xyz = pre_data(i,6:8);
            cpfpre_a15 = importdata('C:\Users\39736\Downloads\cpf_2.00a\cpf_sample_code_v2.00a\cpf_llr_c\apollo15_cpf_061220_35501.utx',' ',3,1);
            pre_data_a15 = cpfpre_a15.data;
            
            xyz_a15 = pre_data_a15(k,6:8);
            format long g
            (double(GC2IT2))*BCRS-xyz';
            norm(xyz)-norm(BCRS)
            (double(GC2IT2))*moon_pos1'*1e3-xyz';
            GC2IT2 * BCRS_a15 - xyz_a15';
            
        end
    end
end
% fMJD_UTC = 53101.3274118751;
% X_itrs = [-1033.4793830, 7901.2952754, 6380.3565958]';
% X_gcrs = [5102.5089592, 6123.0114033, 6378.1369247]';
% GC2IT = IERS.GCRS2ITRS(fMJD_UTC,dt,du,xp,yp);