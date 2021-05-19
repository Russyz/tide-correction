%% 扣除系统延时，得到测量TOF，包含大气延时等还未修正的量
clear all;clc;close all
c=299792458;% m/s
load('rec_c1.mat')
load('rec_c2.mat')
load('rec_c3.mat')
load('rec_c4.mat')
% 探测器各通道延时，先不考虑所有通道的数据，
% 对单个通道进行数据精度进行统计
load('mu_channel1.mat')
load('mu_channel2.mat')
load('mu_channel3.mat')
load('mu_channel4.mat')
% 判断第1通道是否测量到在阈值范围的数据并赋值
if length(rec_c1)==0
    t_c1=[];TOF_c1=[];residual_c1=[];
else
    t_c1=rec_c1(:,1);
    TOF_c1=rec_c1(:,2)-mu_channel1*10^-9;
    residual_c1=rec_c1(:,4)-mu_channel1;
end
% 判断第2通道是否测量到在阈值范围的数据并赋值
if length(rec_c2)==0
    t_c2=[];TOF_c2=[];residual_c2=[]
else
    t_c2=rec_c2(:,1);
    TOF_c2=rec_c2(:,2)-mu_channel2*10^-9;
    residual_c2=rec_c2(:,4)-mu_channel2;
end
% 判断第3通道是否测量到在阈值范围的数据并赋值
if length(rec_c3)==0
    t_c3=[];TOF_c3=[];residual_c3=[];
else
    t_c3=rec_c3(:,1);
    TOF_c3=rec_c3(:,2)-mu_channel3*10^-9;
    residual_c3=rec_c3(:,4)-mu_channel3;
end
% 判断第4通道是否测量到在阈值范围的数据并赋值
if length(rec_c4)==0
    t_c4=[];TOF_c4=[];residual_c4=[];
else
    t_c4=rec_c4(:,1);
    TOF_c4=rec_c4(:,2)-mu_channel4*10^-9;
    residual_c4=rec_c4(:,4)-mu_channel4;
end
% 扣除系统延时得到TOF和O-C残差
t0=[t_c1;t_c2;t_c3;t_c4];
TOF0=[TOF_c1;TOF_c2;TOF_c3;TOF_c4];
residual0=[residual_c1;residual_c2;residual_c3;residual_c4];


% 对信号进行识别
oo = 3;   %设定搜索框的大小：ns
[t1 p1] = sort([t0]); %对记录时刻进行排序
% pos = find(t1>4.7e4);
% t1 = t1(pos);
residual1 = (residual0(p1));
% residual1 = residual1(pos);
TOF1 = TOF0(p1);
% TOF1 = TOF1(pos);
w = 0;
%一阶直线搜索回波信号
for i = -0.1:0.001:0.1 %斜率范围
    for j=-50:0.1:50 %截距
        y=i*(t1-min(t1))+j;
        w1 = [];
        w1 = find(residual1>y&residual1<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %斜率
            q = j; %截距
            w2 = w1; %记录位置信息
        end
    end
end
t2 = t1(w2);
residual2 = residual1(w2);
TOF2 = TOF1(w2);

y=p*(t1-min(t0))+q;


% 对搜索网格内的残差进行2.5σ滤波n次，得到认为是回波信号的TOF3
n=10; %滤波次数
for k=1:n
    residual20=[];
    [p2,S2]=polyfit(t2-min(t2),residual2,1);
    [y_fit2,delta2]=polyval(p2,t2-min(t2),S2);
    for l=1:length(delta2)
        if abs(residual2(l)-y_fit2(l))<2.5*delta2(l)
            residual20(l)=residual2(l);
        end
    end
    w2=find(residual20);
    t3=t2(w2);
    residual3=residual2(w2);
    TOF3=TOF2(w2);
end
figure
plot(t1-min(t0),residual1,'.','MarkerSize',12);
hold on
plot(t3-min(t0),residual3,'.r','MarkerSize',12);
ylim([-100 100])
set(gca,'FontSize',40);
xlabel('Time [s]','fontsize',40);
ylabel('O-C residual [ns]','fontsize',40);
start_hour = floor(min(t3)/3600);
start_second = floor((min(t3) - 3600 * start_hour)/60);
end_hour = floor(max(t3)/3600);
end_second = floor((max(t3) - 3600 * start_hour)/60);

% 对回波信号数据TOF3进行多项式拟合，当n与n+1阶拟合标准差之差小于1mm时，n阶进行拟合，得到标准差
degree=[];
for n=1:100
    [p3,S3]=polyfit(t3,TOF3,n);
    [y_fit3,delta3]=polyval(p3,t3,S3);
    rms(n)=mean(delta3)*c/2;
    degree=[degree;n rms(n)*100];
    if n>1&abs(rms(n)-rms(n-1))<1e-3%判断当n与n+1阶拟合标准差的差值是否小于
        break
    end
end

figure

plot(t3,y_fit3,'-r','MarkerSize',20)
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('RTT/s','fontsize',40);
legend('Observed RRT', 'Polynomial fitting')
% TOF-拟合后曲线得到残差residual4
residual4 = (TOF3 - y_fit3)*1e9;%ns
t4 = t3;
TOF4 = TOF3;
figure


plot(t4,residual4,'.r','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time/seconds','fontsize',40);
ylabel('Residual/ns','fontsize',40);
% figure
% histogram(residual4,'BinWidth',0.02) %柱状图bin大小选取20ps
% set(gca,'FontSize',40);
% xlabel('O-C residual/s','fontsize',40);
% ylabel('Photon counts','fontsize',40);

% %% 生成frd文本，格式需要修改% fid=fopen([num2str(YR),num2str(MONTH),num2str(DATE),temp(end-2:end),'_','_residual.txt'],'w');%写入文件路径
% 
% [r,l]=size(spc);            % 得到矩阵的行数和列数
% for i=1:r
%     for j=1:l
%         fprintf(fid,'%5.12f\t',spc(i,j));
%     end
%     fprintf(fid,'\r\n');
% end
% fclose(fid);

% B = [t4 TOF4];
% % data_name = dir("*.a1*");
% % if length(data_name) == 0
% %     data_name = dir("*.l1*") 
% % end
% %     endsWith(data_name.name,".a11")
% fid=fopen(['d:\LLR\','L21.txt'],'w');%写入文件路径
% 
% [r,c]=size(B);            % 得到矩阵的行数和列数
% for i=1:r
%     for j=1:c
%         fprintf(fid,'%5.12f\t',B(i,j));
%     end
%     fprintf(fid,'\r\n');
% end
% fclose(fid);

%% %% 算出标准点NP
load('pre_data.mat')
c=299792458;%光速
sl=120;%标准点步长选取
gnt=[];
gnd=[];
% NP_xy=[];
for i=1:ceil((max(t4)-min(t4))/sl)
    w5=find((min(t4)+sl*(i-1))<=t4&t4<(min(t4)+sl*i));
    if length(w5) == 0
        break
    else
        t5=t4(w5);
        residual5=residual4(w5);
        TOF5=TOF4(w5);
        [p5,S5]=polyfit(t5,residual5,1);
        [y_fit5,delta5]=polyval(p5,t5,S5);
        RMS1_NP(i)=std(residual5)*c*1e-9/2;
        t_me=mean(t5);
        ct=[];
        for j=1:length(w5)
            ct(j)=abs(t5(j)-t_me);
        end
        w6=find(ct==min(ct));%离平均时刻最近点
        gnt=[gnt;t5(w6(1)) TOF5(w6(1)) mean(delta5)*c*1e-9/2];
        t_np=gnt(:,1);
        
        x=pre_data(2:end,9)*24*3600;
        y=pre_data(2:end,10);
        addpath('/home/z8/文档/潮汐修正/test file')
        delta=(gnt(:,2)-lagint(x,y,gnt(:,1),10))/2*c;
        h=floor(t5(w6)/3600)+8;
    end    
end
%% 检核 @GCRS  -- version 2
% Load P.T.H. data
pth_file = dir('*.pth');
fid = fopen(pth_file.name,'r');
pthdata = fscanf(fid,'%f %f %f %f');
wavelength = 1.064; %um
latitude = 22.34639411111111; % °
height = 0.405713; % km
baroPressure = pthdata(2) * 1e2; % Pa
temperature = pthdata(3); % ℃
RH = pthdata(4); %
data_and_reflector=dir('*.la1');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.la2');
end
temp = data_and_reflector.name;
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));
% coordinate_station_itrs=[-2358691.210,5410611.484,2410087.607];

% Read SP3 data
sp3name=dir('*.sp3');
T=importdata(sp3name.name,' ',22,0);

for i=1:length(T.data(:,1))/3
    sp3t(i,:)=T.data(3*i-2,:);
    sp3xyz(i,:)=T.data(3*i-1,:);
    sp3v_xyz(i,:)=T.data(3*i,:);
end
sp3yyyy=sp3t(:,1);sp3mt=sp3t(:,2);sp3d=sp3t(:,3);sp3h=sp3t(:,4);sp3m=sp3t(:,5);sp3s=sp3t(:,6);  %年月日时分秒
sp3t_s=sp3h*3600+sp3m*60+sp3s; %转化为日积秒
sp3x=sp3xyz(:,1);sp3y=sp3xyz(:,2);sp3z=sp3xyz(:,3); %卫星（x,y,z）坐标

w7 = find(DATE == sp3d);
sp3x_w7=sp3x(w7);sp3y_w7=sp3y(w7);sp3z_w7=sp3z(w7); sp3t_s_w7=sp3t_s(w7);%测量卫星（x,y,z）坐标,时间（t）
sp3t_h_w7=sp3h(w7);sp3t_m_w7=sp3m(w7);sp3t_ss_w7=sp3s(w7);
% Coordinates of satellite in ITRS
xsate_itrs = [sp3x_w7,sp3y_w7,sp3z_w7];
xsate_gcrs = [];
t_up_down = [];
pos_sat_bound = [];
o_minus_c = [];
% if coeff
% coeff = [0 0 ];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    % Calculation tide correction in meter
    addpath('/home/z8/文档/潮汐修正/test file')
    dxtide = solid_tide_dehan(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxptide = solid_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxatide = atm_tide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxoptide = ocean_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607];

    coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
    
    % 以TOF/2作为初始值
    t_up = gnt(i,2)/2;
    for j = 1:4
        % Time bias
        delta_t = 0;
        addpath('/home/z8/文档/潮汐修正/test file')
        sp3_IntpX_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,1),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpY_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,2),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpZ_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,3),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpXYZ_ITRS = [sp3_IntpX_ITRS;sp3_IntpY_ITRS;sp3_IntpZ_ITRS];
        
        IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i));

        format long g
        sp3_IntpXYZ_GCRS= IT2GC * sp3_IntpXYZ_ITRS;

        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        
        dxyz_itrs = sp3_IntpXYZ_ITRS - coordinate_station_itrs;
        
        dxyz_gcrs = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        %         norm(dxyz_itrs) - norm(dxyz_gcrs)
        dxt = dxyz_itrs(1);
        dyt = dxyz_itrs(2);
        dzt = dxyz_itrs(3);
        
        STAT_LONG=113.55421725;
        STAT_LAT=22.34639411111111;
        STAT_HEI=405.713;
        STAT_LONGrad = deg2rad(STAT_LONG);
        STAT_LATrad = deg2rad(STAT_LAT);
        sp3R   = norm(dxyz_gcrs);
        gvs(1) =  cos(STAT_LATrad)* cos(STAT_LONGrad);			% Station X unit vector
        gvs(2) =  cos(STAT_LATrad)* sin(STAT_LONGrad);			% Station Y unit vector
        gvs(3) =  sin(STAT_LATrad);					% Station Z unit vector
        dstn   =  norm(gvs);	% Normalise the unit vectors
        czd    =  (dxt*gvs(1)+dyt*gvs(2)+dzt*gvs(3))/(sp3R*dstn);		% Zenith height component of SAT->STAT vector / vector range
        altc   =  asin(czd)*360.0/(2.0*pi);
        
        correction = atmospheric_refraction(altc, baroPressure, temperature);
        elevationAngle = altc + correction;
        addpath('/home/z8/文档/潮汐修正/大气延时和引力延时修正')
        atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
        addpath('/home/z8/文档/潮汐修正/大气延时和引力延时修正')
        gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
        
        com_corr = 0.251;
        
        t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        t_up = t_up0;
        if (abs(t_up - t_up0) > 1e-12)
            t_up = t_up0;
        else
            t_up = t_up0;
            break
        end
    end
    t_down = gnt(i,2)/2;
    
    for k = 1:4
        WGS84_EARTH_ANGULAR_VELOCITY = 7.292115e-5; % rad/s
        updt = WGS84_EARTH_ANGULAR_VELOCITY * t_up ;
        downdt = WGS84_EARTH_ANGULAR_VELOCITY * t_down;
        %         syms t_down
        IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_down+t_up));
        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        sp3R = norm(dxyz);
        t_down0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        
        f = (sp3R + atm_delay + gravity_delay - com_corr)/c - t_down;
        %solve
        delta_t_down = abs(t_down0 - t_down);
        %         t_down = double(solve(f));
        t_down = t_down0;
        if delta_t_down<1e-12
            break
        end
    end
    
    t_up_down = [t_up_down;t_down,t_up];
    resi = (gnt(i,2) - (t_down + t_up)) * c;
    o_minus_c = [o_minus_c;resi];
    pos_sat_bound = [pos_sat_bound;sp3_IntpXYZ_GCRS'];
end

std(o_minus_c)

%% 测站偏心修正计算
load('pre_data')
pre_utc=pre_data(:,9)*24*3600 - 0.5;
pre_TOF=pre_data(:,10);
dev_pre_TOF = [pre_data(2,10)-pre_data(1,10);pre_data(2:end,10) - pre_data(1:end-1,10)]/2;

vel_radi_meas = lagint(pre_utc,dev_pre_TOF,gnt(:,1),10) * c;
coeff = polyfit(vel_radi_meas,o_minus_c/2,1)
% delta_t = coeff(1)
% bias = coeff(2)
figure
% plot(vel_radi_meas,o_minus_c/2,'.b','MarkerSize',12)
errorbar(vel_radi_meas,o_minus_c/2,gnt(:,3),'.r','MarkerSize',12)
set(gca,'FontSize',40);
xlabel('Radial velocity [meters/second]','fontsize',40);
ylabel('O-C residual [meters]','fontsize',40);
%% 测站偏心修正
% Load P.T.H. data
pth_file = dir('*.pth');
fid = fopen(pth_file.name,'r');
pthdata = fscanf(fid,'%f %f %f %f');
wavelength = 1.064; %um
latitude = 22.34639411111111; % °
height = 0.405713; % km
baroPressure = pthdata(2) * 1e2; % Pa
temperature = pthdata(3); % ℃
RH = pthdata(4); %
data_and_reflector=dir('*.la1');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.la2');
end
temp = data_and_reflector.name;
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));
% coordinate_station_itrs=[-2358691.210,5410611.484,2410087.607];

% Read SP3 data
sp3name=dir('*.sp3');
T=importdata(sp3name.name,' ',22,0);

for i=1:length(T.data(:,1))/3
    sp3t(i,:)=T.data(3*i-2,:);
    sp3xyz(i,:)=T.data(3*i-1,:);
    sp3v_xyz(i,:)=T.data(3*i,:);
end
sp3yyyy=sp3t(:,1);sp3mt=sp3t(:,2);sp3d=sp3t(:,3);sp3h=sp3t(:,4);sp3m=sp3t(:,5);sp3s=sp3t(:,6);  %年月日时分秒
sp3t_s=sp3h*3600+sp3m*60+sp3s; %转化为日积秒
sp3x=sp3xyz(:,1);sp3y=sp3xyz(:,2);sp3z=sp3xyz(:,3); %卫星（x,y,z）坐标

w7 = find(DATE == sp3d);
sp3x_w7=sp3x(w7);sp3y_w7=sp3y(w7);sp3z_w7=sp3z(w7); sp3t_s_w7=sp3t_s(w7);%测量卫星（x,y,z）坐标,时间（t）
sp3t_h_w7=sp3h(w7);sp3t_m_w7=sp3m(w7);sp3t_ss_w7=sp3s(w7);
% Coordinates of satellite in ITRS
xsate_itrs = [sp3x_w7,sp3y_w7,sp3z_w7];
sp3_IntpXYZ_itrs =[];
xsate_gcrs = [];
t_up_down = [];
pos_sat_bound = [];
o_minus_c = [];
% if coeff
tide_cor = [];
corr = [];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    % Calculation tide correction in meter
    dxtide = solid_tide_dehan(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxptide = solid_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxatide = atm_tide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxoptide = ocean_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    delta_x = 0;
    delta_y = 0;
    delta_z = 0;
    
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] + [delta_x, delta_y, delta_z];
    tide_cor = [tide_cor;dxtide + dxatide + dxptide + dxoptide];
    coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
    
    % 以TOF/2作为初始值
    t_up = gnt(i,2)/2;
    for j = 1:4
        % Time bias
        delta_t = 0 + coeff(1);
        sp3_IntpX_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,1),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpY_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,2),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpZ_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,3),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpXYZ_ITRS = [sp3_IntpX_ITRS;sp3_IntpY_ITRS;sp3_IntpZ_ITRS];
        sp3_IntpXYZ_itrs = [sp3_IntpXYZ_itrs;sp3_IntpX_ITRS,sp3_IntpY_ITRS,sp3_IntpZ_ITRS];
        IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i));

        format long g
        sp3_IntpXYZ_GCRS= IT2GC * sp3_IntpXYZ_ITRS;

        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        
        dxyz_itrs = sp3_IntpXYZ_ITRS - coordinate_station_itrs;
        
        dxyz_gcrs = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        %         norm(dxyz_itrs) - norm(dxyz_gcrs)
        dxt = dxyz_itrs(1);
        dyt = dxyz_itrs(2);
        dzt = dxyz_itrs(3);
        
        STAT_LONG=113.55421725;
        STAT_LAT=22.34639411111111;
        STAT_HEI=405.713;
        STAT_LONGrad = deg2rad(STAT_LONG);
        STAT_LATrad = deg2rad(STAT_LAT);
        sp3R   = norm(dxyz_gcrs);
        gvs(1) =  cos(STAT_LATrad)* cos(STAT_LONGrad);			% Station X unit vector
        gvs(2) =  cos(STAT_LATrad)* sin(STAT_LONGrad);			% Station Y unit vector
        gvs(3) =  sin(STAT_LATrad);					% Station Z unit vector
        dstn   =  norm(gvs);	% Normalise the unit vectors
        czd    =  (dxt*gvs(1)+dyt*gvs(2)+dzt*gvs(3))/(sp3R*dstn);		% Zenith height component of SAT->STAT vector / vector range
        altc   =  asin(czd)*360.0/(2.0*pi);
        
        correction = atmospheric_refraction(altc, baroPressure, temperature);
        elevationAngle = altc + correction;
        
        atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
        
        gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
        
        com_corr = 0.251;
        
        t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        t_up = t_up0;
        if (abs(t_up - t_up0) > 1e-12)
            t_up = t_up0;
        else
            t_up = t_up0;
            break
        end
    end
    t_down = gnt(i,2)/2;
    
    for k = 1:4
        WGS84_EARTH_ANGULAR_VELOCITY = 7.292115e-5; % rad/s
        updt = WGS84_EARTH_ANGULAR_VELOCITY * t_up ;
        downdt = WGS84_EARTH_ANGULAR_VELOCITY * t_down;
        %         syms t_down
        IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_down+t_up));
        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        sp3R = norm(dxyz);
        t_down0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        
        f = (sp3R + atm_delay + gravity_delay - com_corr)/c - t_down;
        %solve
        delta_t_down = abs(t_down0 - t_down);
        %         t_down = double(solve(f));
        t_down = t_down0;
        if delta_t_down<1e-12
            break
        end
    end
    
    t_up_down = [t_up_down;t_down,t_up];
    resi = (gnt(i,2) - (t_down + t_up)) * c;
    o_minus_c = [o_minus_c;resi];
    pos_sat_bound = [pos_sat_bound;sp3_IntpXYZ_GCRS'];
    corr = [corr;atm_delay + gravity_delay - com_corr];
    
    
end
%% station coordination calculation
A = [];
for i = 1:length(t_up_down(:,1))
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607];
    coordinate_station_itrs = (coordinate_station_itrs0 - tide_cor(i))';
    
    IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_up_down(i,1)));
    coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
    sp3_IntpXYZ_GCRS = pos_sat_bound;
    dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
    sp3R = norm(dxyz);
    A = [A;dxyz(1)/sp3R,dxyz(2)/sp3R,dxyz(3)/sp3R];
    
end
% b = regress(o_minus_c(1:i)/2,A)


%% station coordination calculation version--2
A = [];
B = [];
C = [];
D = [];
for i = 1:length(gnt(:,1))
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607];
    coordinate_station_itrs = (coordinate_station_itrs0 - tide_cor(i))';   
    
    
   
    dxyz = sp3_IntpXYZ_itrs(i,:) - coordinate_station_itrs';
    sp3R = norm(dxyz);
    A = [A;dxyz(1)/sp3R,dxyz(2)/sp3R,dxyz(3)/sp3R];
    B = [B;dxyz(1)/sp3R];
    C = [C;dxyz(2)/sp3R];
    D = [D;dxyz(3)/sp3R];
end
inv(A' * A) * A' * (o_minus_c/2)
% b = regress(o_minus_c(1:i)/2,A)
%% Gradient descent
delta_x = 0;
delta_y = 0;
delta_z = 0;
delta1 = 0.001;

lambda = 0.001;
F = [];
for n =1:100000
    E = norm(A*[delta_x; delta_y; delta_z] - o_minus_c/2);
    grad_E_x = (norm(A*[delta_x+delta1; delta_y; delta_z] - o_minus_c/2) - E)/delta1;
    grad_E_y = (norm(A*[delta_x; delta_y+delta1; delta_z] - o_minus_c/2) - E)/delta1;
    grad_E_z = (norm(A*[delta_x; delta_y; delta_z+delta1] - o_minus_c/2) - E)/delta1;
    if (lambda * grad_E_x) < 0
        E = E + lambda * grad_E_x;   
        delta_x = delta_x + lambda;
    else
        E = E - lambda * grad_E_x;
        delta_x = delta_x - lambda;
    end    
    if (lambda * grad_E_y) < 0
        E = E + lambda * grad_E_y;
        delta_y = delta_y + lambda;
    else
        E = E - lambda * grad_E_y;
        delta_y = delta_y - lambda;
    end    
    if (lambda * grad_E_z) < 0
        E = E + lambda * grad_E_z;
        delta_z = delta_z + lambda;
    else
        E = E - lambda * grad_E_z;
        delta_z = delta_z - lambda;
    end
%     delta_x = delta_x - lambda;
%     delta_y = delta_y - lambda;
%     delta_z = delta_z - lambda;
    F= [F;n, E];
%     if  n > 2 & abs(F(end-2)) < abs(E)        
%         
%         break
%         
%     end
    
    
    
end
minp = F(find(abs(F(:,2)) == min(abs(F(:,2)))),1);
delta_x = delta_x
delta_y = delta_y
delta_z = delta_z
% delta_x = 0.754;
% delta_y = 0.754;
% delta_z = 0.754;

%% 测站偏心修正
% Load P.T.H. data
pth_file = dir('*.pth');
fid = fopen(pth_file.name,'r');
pthdata = fscanf(fid,'%f %f %f %f');
wavelength = 1.064; %um
latitude = 22.34639411111111; % °
height = 0.405713; % km
baroPressure = pthdata(2) * 1e2; % Pa
temperature = pthdata(3); % ℃
RH = pthdata(4); %
data_and_reflector=dir('*.la1');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.la2');
end
temp = data_and_reflector.name;
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));
% coordinate_station_itrs=[-2358691.210,5410611.484,2410087.607];

% Read SP3 data
sp3name=dir('*.sp3');
T=importdata(sp3name.name,' ',22,0);

for i=1:length(T.data(:,1))/3
    sp3t(i,:)=T.data(3*i-2,:);
    sp3xyz(i,:)=T.data(3*i-1,:);
    sp3v_xyz(i,:)=T.data(3*i,:);
end
sp3yyyy=sp3t(:,1);sp3mt=sp3t(:,2);sp3d=sp3t(:,3);sp3h=sp3t(:,4);sp3m=sp3t(:,5);sp3s=sp3t(:,6);  %年月日时分秒
sp3t_s=sp3h*3600+sp3m*60+sp3s; %转化为日积秒
sp3x=sp3xyz(:,1);sp3y=sp3xyz(:,2);sp3z=sp3xyz(:,3); %卫星（x,y,z）坐标

w7 = find(DATE == sp3d);
sp3x_w7=sp3x(w7);sp3y_w7=sp3y(w7);sp3z_w7=sp3z(w7); sp3t_s_w7=sp3t_s(w7);%测量卫星（x,y,z）坐标,时间（t）
sp3t_h_w7=sp3h(w7);sp3t_m_w7=sp3m(w7);sp3t_ss_w7=sp3s(w7);
% Coordinates of satellite in ITRS
xsate_itrs = [sp3x_w7,sp3y_w7,sp3z_w7];
sp3_IntpXYZ_itrs =[];
xsate_gcrs = [];
t_up_down = [];
pos_sat_bound = [];
o_minus_c = [];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    % Calculation tide correction in meter
    dxtide = solid_tide_dehan(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxptide = solid_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxatide = atm_tide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxoptide = ocean_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     delta_x = 0;
%     delta_y = 0;
%     delta_z = 0;
    
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] - [delta_x, delta_y, delta_z];
    tide_cor = [tide_cor;dxtide + dxatide + dxptide + dxoptide];
    coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
    
    % 以TOF/2作为初始值
    t_up = gnt(i,2)/2;
    for j = 1:4
        % Time bias
        delta_t = 0;
        sp3_IntpX_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,1),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpY_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,2),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpZ_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,3),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpXYZ_ITRS = [sp3_IntpX_ITRS;sp3_IntpY_ITRS;sp3_IntpZ_ITRS];
        sp3_IntpXYZ_itrs = [sp3_IntpXYZ_itrs;sp3_IntpX_ITRS,sp3_IntpY_ITRS,sp3_IntpZ_ITRS];
        IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i));

        format long g
        sp3_IntpXYZ_GCRS= IT2GC * sp3_IntpXYZ_ITRS;

        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        
        dxyz_itrs = sp3_IntpXYZ_ITRS - coordinate_station_itrs;
        
        dxyz_gcrs = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        %         norm(dxyz_itrs) - norm(dxyz_gcrs)
        dxt = dxyz_itrs(1);
        dyt = dxyz_itrs(2);
        dzt = dxyz_itrs(3);
        
        STAT_LONG=113.55421725;
        STAT_LAT=22.34639411111111;
        STAT_HEI=405.713;
        STAT_LONGrad = deg2rad(STAT_LONG);
        STAT_LATrad = deg2rad(STAT_LAT);
        sp3R   = norm(dxyz_gcrs);
        gvs(1) =  cos(STAT_LATrad)* cos(STAT_LONGrad);			% Station X unit vector
        gvs(2) =  cos(STAT_LATrad)* sin(STAT_LONGrad);			% Station Y unit vector
        gvs(3) =  sin(STAT_LATrad);					% Station Z unit vector
        dstn   =  norm(gvs);	% Normalise the unit vectors
        czd    =  (dxt*gvs(1)+dyt*gvs(2)+dzt*gvs(3))/(sp3R*dstn);		% Zenith height component of SAT->STAT vector / vector range
        altc   =  asin(czd)*360.0/(2.0*pi);
        
        correction = atmospheric_refraction(altc, baroPressure, temperature);
        elevationAngle = altc + correction;
        
        atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
        
        gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
        
        com_corr = 0.251;
        
        t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        t_up = t_up0;
        if (abs(t_up - t_up0) > 1e-12)
            t_up = t_up0;
        else
            t_up = t_up0;
            break
        end
    end
    t_down = gnt(i,2)/2;
    
    for k = 1:4
        WGS84_EARTH_ANGULAR_VELOCITY = 7.292115e-5; % rad/s
        updt = WGS84_EARTH_ANGULAR_VELOCITY * t_up ;
        downdt = WGS84_EARTH_ANGULAR_VELOCITY * t_down;
        %         syms t_down
        IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_down+t_up));
        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        sp3R = norm(dxyz);
        t_down0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        
        f = (sp3R + atm_delay + gravity_delay - com_corr)/c - t_down;
        %solve
        delta_t_down = abs(t_down0 - t_down);
        %         t_down = double(solve(f));
        t_down = t_down0;
        if delta_t_down<1e-12
            break
        end
    end
    
    t_up_down = [t_up_down;t_down,t_up];
    resi = (gnt(i,2) - (t_down + t_up)) * c;
    o_minus_c = [o_minus_c;resi];  
    
end

%% 测站偏心修正计算
load('pre_data')
pre_utc=pre_data(:,9)*24*3600 - 0.5;
pre_TOF=pre_data(:,10);
dev_pre_TOF = [pre_data(2,10)-pre_data(1,10);pre_data(2:end,10) - pre_data(1:end-1,10)]/2;

vel_radi_meas = lagint(pre_utc,dev_pre_TOF,gnt(:,1),10) * c;
coeff = polyfit(vel_radi_meas,o_minus_c/2,1)
% delta_t = coeff(1)
% bias = coeff(2)
figure
% plot(vel_radi_meas,o_minus_c/2,'.b','MarkerSize',12)
errorbar(vel_radi_meas,o_minus_c/2,gnt(:,3),'.r','MarkerSize',12)
set(gca,'FontSize',40);
xlabel('Radial velocity [meters/second]','fontsize',40);
ylabel('O-C residual [meters]','fontsize',40);
%% 测站偏心修正
% Load P.T.H. data
pth_file = dir('*.pth');
fid = fopen(pth_file.name,'r');
pthdata = fscanf(fid,'%f %f %f %f');
wavelength = 1.064; %um
latitude = 22.34639411111111; % °
height = 0.405713; % km
baroPressure = pthdata(2) * 1e2; % Pa
temperature = pthdata(3); % ℃
RH = pthdata(4); %
data_and_reflector=dir('*.la1');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.la2');
end
temp = data_and_reflector.name;
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));
% coordinate_station_itrs=[-2358691.210,5410611.484,2410087.607];

% Read SP3 data
sp3name=dir('*.sp3');
T=importdata(sp3name.name,' ',22,0);

for i=1:length(T.data(:,1))/3
    sp3t(i,:)=T.data(3*i-2,:);
    sp3xyz(i,:)=T.data(3*i-1,:);
    sp3v_xyz(i,:)=T.data(3*i,:);
end
sp3yyyy=sp3t(:,1);sp3mt=sp3t(:,2);sp3d=sp3t(:,3);sp3h=sp3t(:,4);sp3m=sp3t(:,5);sp3s=sp3t(:,6);  %年月日时分秒
sp3t_s=sp3h*3600+sp3m*60+sp3s; %转化为日积秒
sp3x=sp3xyz(:,1);sp3y=sp3xyz(:,2);sp3z=sp3xyz(:,3); %卫星（x,y,z）坐标

w7 = find(DATE == sp3d);
sp3x_w7=sp3x(w7);sp3y_w7=sp3y(w7);sp3z_w7=sp3z(w7); sp3t_s_w7=sp3t_s(w7);%测量卫星（x,y,z）坐标,时间（t）
sp3t_h_w7=sp3h(w7);sp3t_m_w7=sp3m(w7);sp3t_ss_w7=sp3s(w7);
% Coordinates of satellite in ITRS
xsate_itrs = [sp3x_w7,sp3y_w7,sp3z_w7];
sp3_IntpXYZ_itrs =[];
xsate_gcrs = [];
t_up_down = [];
pos_sat_bound = [];
o_minus_c = [];
% if coeff
tide_cor = [];
corr = [];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    % Calculation tide correction in meter
    dxtide = solid_tide_dehan(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxptide = solid_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxatide = atm_tide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
    dxoptide = ocean_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     delta_x = 0;
%     delta_y = 0;
%     delta_z = 0;
    
    coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] - [delta_x, delta_y, delta_z];
    tide_cor = [tide_cor;dxtide + dxatide + dxptide + dxoptide];
    coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
    
    % 以TOF/2作为初始值
    t_up = gnt(i,2)/2;
    for j = 1:4
        % Time bias
        delta_t = 0 + coeff(1);
        sp3_IntpX_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,1),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpY_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,2),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpZ_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,3),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
        sp3_IntpXYZ_ITRS = [sp3_IntpX_ITRS;sp3_IntpY_ITRS;sp3_IntpZ_ITRS];
        sp3_IntpXYZ_itrs = [sp3_IntpXYZ_itrs;sp3_IntpX_ITRS,sp3_IntpY_ITRS,sp3_IntpZ_ITRS];
        IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i));

        format long g
        sp3_IntpXYZ_GCRS= IT2GC * sp3_IntpXYZ_ITRS;

        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        
        dxyz_itrs = sp3_IntpXYZ_ITRS - coordinate_station_itrs;
        
        dxyz_gcrs = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        %         norm(dxyz_itrs) - norm(dxyz_gcrs)
        dxt = dxyz_itrs(1);
        dyt = dxyz_itrs(2);
        dzt = dxyz_itrs(3);
        
        STAT_LONG=113.55421725;
        STAT_LAT=22.34639411111111;
        STAT_HEI=405.713;
        STAT_LONGrad = deg2rad(STAT_LONG);
        STAT_LATrad = deg2rad(STAT_LAT);
        sp3R   = norm(dxyz_gcrs);
        gvs(1) =  cos(STAT_LATrad)* cos(STAT_LONGrad);			% Station X unit vector
        gvs(2) =  cos(STAT_LATrad)* sin(STAT_LONGrad);			% Station Y unit vector
        gvs(3) =  sin(STAT_LATrad);					% Station Z unit vector
        dstn   =  norm(gvs);	% Normalise the unit vectors
        czd    =  (dxt*gvs(1)+dyt*gvs(2)+dzt*gvs(3))/(sp3R*dstn);		% Zenith height component of SAT->STAT vector / vector range
        altc   =  asin(czd)*360.0/(2.0*pi);
        
        correction = atmospheric_refraction(altc, baroPressure, temperature);
        elevationAngle = altc + correction;
        
        atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
        
        gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
        
        com_corr = 0.251;
        
        t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        t_up = t_up0;
        if (abs(t_up - t_up0) > 1e-12)
            t_up = t_up0;
        else
            t_up = t_up0;
            break
        end
    end
    t_down = gnt(i,2)/2;
    
    for k = 1:4
        WGS84_EARTH_ANGULAR_VELOCITY = 7.292115e-5; % rad/s
        updt = WGS84_EARTH_ANGULAR_VELOCITY * t_up ;
        downdt = WGS84_EARTH_ANGULAR_VELOCITY * t_down;
        %         syms t_down
        IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_down+t_up));
        coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
        dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
        sp3R = norm(dxyz);
        t_down0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        
        f = (sp3R + atm_delay + gravity_delay - com_corr)/c - t_down;
        %solve
        delta_t_down = abs(t_down0 - t_down);
        %         t_down = double(solve(f));
        t_down = t_down0;
        if delta_t_down<1e-12
            break
        end
    end
    
    t_up_down = [t_up_down;t_down,t_up];
    resi = (gnt(i,2) - (t_down + t_up)) * c;
    o_minus_c = [o_minus_c;resi];
    pos_sat_bound = [pos_sat_bound;sp3_IntpXYZ_GCRS'];
    corr = [corr;atm_delay + gravity_delay - com_corr];
    
    
end
figure
errorbar(vel_radi_meas,o_minus_c/2,gnt(:,3),'.r','MarkerSize',12)
set(gca,'FontSize',40);
xlabel('Radial velocity [meters/second]','fontsize',40);
ylabel('O-C residual [meters]','fontsize',40);
%% station coordination calculation version--3
% syms delta_x delta_y delta_z
% fun = inv(cov(A))*[delta_x; delta_y; delta_z]-cov(o_minus_c);
% inv(cov(A))*cov(o_minus_c)
% [delta_x, delta_y, delta_z] = solve(fun==0);
% [double(delta_x),double(delta_y),double(delta_z)]

%% 站址修正
% c=299792458;% m/s
% spc = [];
% 
% 
% for i = 1:length(t_up_down(:,1))
%     delta_RTT = 5;
%     delta_y = 0;
%     delta_z = 0;
%     for delta_x = -20:1:100        
%         
%         coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] + [delta_x, delta_y, delta_z];
%         coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)));
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_up0 = (sp3R + corr(i))/c;
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_up_down(i,1)+t_up_down(i,2)));
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_down0 = (sp3R + corr(i))/c;
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         delta_RTT0 = (t_up0 + t_down0 - gnt(i,2))*c;
%         if abs(delta_RTT0) < abs(delta_RTT)
%             delta_RTT = delta_RTT0;
%             delta_x0 = delta_x;
%         end
%     end
%     
%     
%     delta_x = delta_x0;
%     delta_z = 0;
%     delta_RTT = 6;
%     for delta_y = -20:1:20 
%         
%         
%         coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] + [delta_x, delta_y, delta_z];
%         coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)));
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_up0 = (sp3R + corr(i))/c;
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_up_down(i,1)+t_up_down(i,2)));
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_down0 = (sp3R + corr(i))/c;
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         delta_RTT0 = (t_up0 + t_down0 - gnt(i,2))*c;
%         if abs(delta_RTT0) < abs(delta_RTT)
%             delta_RTT = delta_RTT0;
%             delta_y0 = delta_y;
%         end
%     end
%     
%     
%     delta_x = delta_x0;
%     delta_y = delta_y0;
%     delta_RTT = 5;
%     for delta_z = -20:1:20 
%         
%         
%         coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] + [delta_x, delta_y, delta_z];
%         coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)));
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_up0 = (sp3R + corr(i))/c;
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_up_down(i,1)+t_up_down(i,2)));
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%         sp3R = norm(dxyz);
%         t_down0 = (sp3R + corr(i))/c;
%         sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%         delta_RTT0 = (t_up0 + t_down0 - gnt(i,2))*c;
%         if abs(delta_RTT0) < abs(delta_RTT)
%             delta_RTT = delta_RTT0;
%             delta_z0 = delta_z;
%         end
%     end
%     spc = [spc; delta_x0, delta_y0, delta_z0];
% end
% 
% fid=fopen([num2str(YR),num2str(MONTH),num2str(DATE),temp(end-2:end),'_','_residual.txt'],'w');%写入文件路径
% 
% [r,l]=size(spc);            % 得到矩阵的行数和列数
% for i=1:r
%     for j=1:l
%         fprintf(fid,'%2.3f\t',spc(i,j));
%     end
%     fprintf(fid,'\r\n');
% end
% fclose(fid);

%% 站址修正
% c=299792458;% m/s
% spc = [];
% for delta_x = -1:0.01:1
% %     for delta_y = -1:0.01:1
% %         for delta_z = -1:0.01:1
%             delta_RTT = 5;
%             for i = 1:length(t_up_down(:,1))
%                 %                 delta_RTT = 5;
%                 coordinate_station_itrs0=[-2358691.210,5410611.484,2410087.607] + [delta_x, delta_x, delta_x];%[delta_x, delta_y, delta_z];
%                 coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
%                 IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)));
%                 sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%                 coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%                 dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%                 sp3R = norm(dxyz);
%                 t_up0 = (sp3R + corr(i))/c;
%                 IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_up_down(i,1)+t_up_down(i,2)));
%                 coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%                 dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs';
%                 sp3R = norm(dxyz);
%                 t_down0 = (sp3R + corr(i))/c;
%                 sp3_IntpXYZ_GCRS = pos_sat_bound(i,:);
%                 delta_RTT0 = (t_up_down(i,1) + t_up_down(i,2) - gnt(i,2))*c;
%                 if delta_RTT0 < delta_RTT
%                     delta_RTT = delta_RTT0;
%                 end                
%             end
%             if abs(delta_RTT) < 5e-1                
% %                 spc = [spc;delta_RTT,delta_x, delta_y, delta_z];
%                 spc = [spc;delta_RTT,delta_x, delta_x, delta_x]
%             end
% %         end
% %     end
% %     delta_x
% end


% data_name = dir("*.a1*");
% if length(data_name) == 0
%     data_name = dir("*.l1*") 
% end
%     endsWith(data_name.name,".a11")

% fid=fopen([num2str(YR),num2str(MONTH),num2str(DATE),temp(end-2:end),'_','_residual.txt'],'w');%写入文件路径
% 
% [r,l]=size(spc);            % 得到矩阵的行数和列数
% for i=1:r
%     for j=1:l
%         fprintf(fid,'%5.12f\t',spc(i,j));
%     end
%     fprintf(fid,'\r\n');
% end
% fclose(fid);
%% 其他站检核
% % Example 1: YARL
% % importdata('D:\npt\la1\lageos1_20191116.npt', ' ', 12,0)
% c=299792458;% m/s
% std1 = 0;
% data = [11 21308.400550200000     0.047387701325 std1 2  120.0     44   35.0   0.088  -0.083      -1.0   7.33 0                                                                
% 11 21399.000549900000     0.048710360126 std1 2  120.0     24   34.0   0.024  -0.490      -1.0   4.00 0                                                            
% 11 21496.400550999999     0.050220552225 std1 2  120.0      1   35.0   0.000  -3.000      -1.0   0.17 0                                                             
% 11 21602.400549800001     0.051954488686 std1 2  120.0      2    9.0   0.000  -2.000      -1.0   0.33 0];
% gnt=[data(:,2) data(:,3) data(:,7)*c*1e-12/2];
% pthdata = [20 21308.401  981.30 310.40  23. 0  
% 20 21399.001  981.30 310.50  23. 0    
% 20 21496.401  981.30 310.50  23. 0  
% 20 21602.401  981.30 310.00  23. 0];
% wavelength = 0.532; %um
% latitude = -29.0464; % °
% height =  244e-3; % km
% 
% pth_time = pthdata(:,2);
% baroPressure_list = pthdata(:,3) * 1e2; % Pa
% temperature_list = pthdata(:,4)-273.15; % ℃
% RH_list = pthdata(:,5); %
% 
% % wavelength = 1.064; %um
% % latitude = 22.34639411111111; % °
% % height = 0.405713; % km
% data_and_reflector=dir('*.la1');
% if length(data_and_reflector) == 0
%     data_and_reflector=dir('*.la2');
% end
% temp = data_and_reflector.name;
% YR = str2num(temp(1:4));
% MONTH = str2num(temp(5:6));
% DATE = str2num(temp(7:8));
% 
% sp3name=dir('*.sp3');
% T=importdata(sp3name.name,' ',22,0);
% 
% for i=1:length(T.data(:,1))/3
%     sp3t(i,:)=T.data(3*i-2,:);
%     sp3xyz(i,:)=T.data(3*i-1,:);
%     sp3v_xyz(i,:)=T.data(3*i,:);
% end
% sp3yyyy=sp3t(:,1);sp3mt=sp3t(:,2);sp3d=sp3t(:,3);sp3h=sp3t(:,4);sp3m=sp3t(:,5);sp3s=sp3t(:,6);  %年月日时分秒
% sp3t_s=sp3h*3600+sp3m*60+sp3s; %转化为日积秒
% sp3x=sp3xyz(:,1);sp3y=sp3xyz(:,2);sp3z=sp3xyz(:,3); %卫星（x,y,z）坐标
% 
% w7 = find(DATE == sp3d);
% sp3x_w7=sp3x(w7);sp3y_w7=sp3y(w7);sp3z_w7=sp3z(w7); sp3t_s_w7=sp3t_s(w7);%测量卫星（x,y,z）坐标,时间（t）
% sp3t_h_w7=sp3h(w7);sp3t_m_w7=sp3m(w7);sp3t_ss_w7=sp3s(w7);
% % Coordinates of satellite in ITRS
% xsate_itrs = [sp3x_w7,sp3y_w7,sp3z_w7];
% xsate_gcrs = [];
% t_up_down = [];
% pos_sat_bound = [];
% o_minus_c = [];
% % coeff = [0 0 ];
% for i = 1:length(gnt(:,1))
%     %     syms delta_x delta_y delta_z
%     hour(i) = floor(gnt(i,1)/3600);
%     minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);  
%     sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
%     % Calculation tide correction in meter
%     dxtide = (YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     dxotide = ocean_tide_hardisp(YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     dxptide = solid_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     dxatide = atm_tide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
%     dxoptide = ocean_ptide(YR,MONTH,DATE,hour(i),minu(i),sec(i));
% %     coordinate_station_itrs0=[-2830744.471,4676580.282,3275072.816];
%     coordinate_station_itrs0=[-2389008,5043332,-3078526];
%     coordinate_station_itrs = (coordinate_station_itrs0 - dxtide - dxatide - dxptide - dxoptide)';
%     
%     % 以TOF/2作为初始值
%     t_up = gnt(i,2)/2;
%     for j = 1:4
%         % Time bias
%         delta_t = 0;
%         sp3_IntpX_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,1),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
%         sp3_IntpY_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,2),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
%         sp3_IntpZ_ITRS = lagint(sp3t_s_w7,xsate_itrs(:,3),gnt(i,1)+t_up+1e-7 + delta_t, 10)*1e3;
%         sp3_IntpXYZ_ITRS = [sp3_IntpX_ITRS;sp3_IntpY_ITRS;sp3_IntpZ_ITRS];
%         IT2GC = ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i));
% 
%         format long g
%         sp3_IntpXYZ_GCRS= IT2GC * sp3_IntpXYZ_ITRS;
% 
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         
%         dxyz_itrs = sp3_IntpXYZ_ITRS - coordinate_station_itrs;
%         
%         dxyz_gcrs = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
%         %         norm(dxyz_itrs) - norm(dxyz_gcrs)
%         dxt = dxyz_itrs(1);
%         dyt = dxyz_itrs(2);
%         dzt = dxyz_itrs(3);
%         
%         STAT_LONG=113.55421725;
%         STAT_LAT=22.34639411111111;
%         STAT_HEI=405.713;
%         STAT_LONGrad = deg2rad(STAT_LONG);
%         STAT_LATrad = deg2rad(STAT_LAT);
%         sp3R   = norm(dxyz_gcrs);
%         gvs(1) =  cos(STAT_LATrad)* cos(STAT_LONGrad);			% Station X unit vector
%         gvs(2) =  cos(STAT_LATrad)* sin(STAT_LONGrad);			% Station Y unit vector
%         gvs(3) =  sin(STAT_LATrad);					% Station Z unit vector
%         dstn   =  norm(gvs);	% Normalise the unit vectors
%         czd    =  (dxt*gvs(1)+dyt*gvs(2)+dzt*gvs(3))/(sp3R*dstn);		% Zenith height component of SAT->STAT vector / vector range
%         altc   =  asin(czd)*360.0/(2.0*pi);
%         baroPressure = interp1(pth_time,baroPressure_list,gnt(i,1));
%         temperature = interp1(pth_time,temperature_list,gnt(i,1));
%         RH = interp1(pth_time,RH_list,gnt(i,1));
%         correction = atmospheric_refraction(altc, baroPressure, temperature);
%         elevationAngle = altc + correction;
%         
%         atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
%         
%         gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
%         
%         com_corr = 0.251;
%         
%         t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
%         t_up = t_up0;
%         if (abs(t_up - t_up0) > 1e-9)
%             t_up = t_up0;
%         else
%             t_up = t_up0;
%             break
%         end
%     end
%     t_down = gnt(i,2)/2;
%     
%     for k = 1:4
%         WGS84_EARTH_ANGULAR_VELOCITY = 7.292115e-5; % rad/s
%         updt = WGS84_EARTH_ANGULAR_VELOCITY * t_up ;
%         downdt = WGS84_EARTH_ANGULAR_VELOCITY * t_down;
%         %         syms t_down
%         IT2GC = (ITRS2GCRS(YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_down+t_up));
%         coordinate_station_gcrs = IT2GC * coordinate_station_itrs;
%         dxyz = sp3_IntpXYZ_GCRS - coordinate_station_gcrs;
%         sp3R = norm(dxyz);
%         t_down0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
%         
%         f = (sp3R + atm_delay + gravity_delay - com_corr)/c - t_down;
%         %solve
%         delta_t_down = abs(t_down0 - t_down);
%         %         t_down = double(solve(f));
%         t_down = t_down0;
%         if delta_t_down<1e-9
%             break
%         end
%     end
%     
%     t_up_down = [t_up_down;t_down,t_up];
%     resi = (gnt(i,2) - (t_down + t_up)) * c;
%     o_minus_c = [o_minus_c;resi];
%     pos_sat_bound = [pos_sat_bound;sp3_IntpXYZ_GCRS'];
% end
% 
% std(o_minus_c)

