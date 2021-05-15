%% 读数据并存储
clc;clear all

emission = dir('*.c00');
receive1 = dir('*.c01');
receive2 = dir('*.c02');
receive3 = dir('*.c03');
receive4 = dir('*.c04');

Te=importdata(emission.name,' ',2,0) ;
Tc1=importdata(receive1.name, ' ', 2,0) ;
Tc2=importdata(receive2.name, ' ', 2,0) ;
Tc3=importdata(receive3.name, ' ', 2,0) ;
Tc4=importdata(receive4.name, ' ', 2,0) ;
% 各通道回波时刻值
channel1=Tc1.data(:,1)+Tc1.data(:,2);
channel10=Tc1.data(:,1);channel11=Tc1.data(:,2);
channel2=Tc2.data(:,1)+Tc2.data(:,2);
channel20=Tc2.data(:,1);channel21=Tc2.data(:,2);
channel3=Tc3.data(:,1)+Tc3.data(:,2);
channel30=Tc3.data(:,1);channel31=Tc3.data(:,2);
% channel3=[];
% channel30=[];channel31=[];
channel4=Tc4.data(:,1)+Tc4.data(:,2);
channel40=Tc4.data(:,1);channel41=Tc4.data(:,2);
save('channel1.mat','channel1')
save('channel10.mat','channel10')
save('channel11.mat','channel11')
save('channel2.mat','channel2')
save('channel20.mat','channel20')
save('channel21.mat','channel21')
save('channel3.mat','channel3')
save('channel30.mat','channel30')
save('channel31.mat','channel31')
save('channel4.mat','channel4')
save('channel40.mat','channel40')
save('channel41.mat','channel41')
% 主波时刻
emission_time=Te.data(:,1)+Te.data(:,2);
emission_time0=Te.data(:,1);
emission_time1=Te.data(:,2);
save('emission_time.mat','emission_time')
save('emission_time0.mat','emission_time0')
save('emission_time1.mat','emission_time1')
p_dot = find(emission.name == '.');
sate_name = emission.name(1:p_dot(2)-1);
predict = dir(sate_name); %修改预报卫星简称
P=importdata(predict.name, ' ', 2,0);
pre_data=P.data;
save('pre_data.mat','pre_data')
%% 主回波匹配
% % 读取预报
% load('pre_data.mat')
% % 读取主波
% load('emission_time.mat')
% load('emission_time0.mat')
% load('emission_time1.mat')
% % 读取回波1
% load('channel1.mat')
% load('channel10.mat')
% load('channel11.mat')
% % 读取回波2
% load('channel2.mat')
% load('channel20.mat')
% load('channel21.mat')
% % 读取回波3
% load('channel3.mat')
% load('channel30.mat')
% load('channel31.mat')
% % 读取回波4
% load('channel4.mat')
% load('channel40.mat')
% load('channel41.mat')

% 预报TOF
pre_utc=pre_data(:,9)*24*3600;
pre_TOF=pre_data(:,10);
rec_c1=[];
rec_c2=[];
rec_c3=[];
rec_c4=[];
w=[];
% 预报插值
pre_TOF=lagint(pre_utc,pre_TOF,emission_time,10)';
t_pre_receive = emission_time + pre_TOF;
% 匹配在阈值内残差
for m=1:length(emission_time)
    t1=[];
    w1=find(abs(channel1-t_pre_receive(m))<=1e-6);
    if ~isempty(w1)
        residual1=channel10(w1)-emission_time0(m)+channel11(w1)-emission_time1(m)-pre_TOF(m); %残差s
        TOF1=residual1 + pre_TOF(m); %飞行时间s
        for n=1:length(w1)
            rec_c1=[rec_c1;emission_time(m) TOF1(n) pre_TOF(m) residual1(n)*1e9];
        end
    end
    w2=find(abs(channel2-t_pre_receive(m))<=1e-6);
    if ~isempty(w2)
        residual2=(channel20(w2)-emission_time0(m)+channel21(w2)-emission_time1(m)-pre_TOF(m));% 残差ns
        TOF2=residual2 + pre_TOF(m);% 飞行时间s
        for n=1:length(w2)
            rec_c2=[rec_c2;emission_time(m) TOF2(n) pre_TOF(m) residual2(n)*1e9];
        end
    end
    w3=find(abs(channel3-t_pre_receive(m))<=1e-6);
    if ~isempty(w3)
        residual3=channel30(w3)-emission_time0(m)+channel31(w3)-emission_time1(m)-pre_TOF(m);% 残差s
        TOF3=residual3+pre_TOF(m);% 飞行时间s
        for n=1:length(w3)
            rec_c3=[rec_c3;emission_time(m) TOF3(n) pre_TOF(m) residual3(n)*1e9];
        end
    end
    w4=find(abs(channel4-t_pre_receive(m))<=1e-6);
    if ~isempty(w4)
        t4=emission_time(m);% 日积秒 s
        residual4=channel40(w4)-emission_time0(m)+channel41(w4)-emission_time1(m)-pre_TOF(m);% 残差ns
        TOF4=residual4 + pre_TOF(m);% 飞行时间s
        for n=1:length(w4)
            rec_c4=[rec_c4;emission_time(m) TOF4(n) pre_TOF(m) residual4(n)*1e9];
        end
    end
end
save('rec_c1.mat','rec_c1')
save('rec_c2.mat','rec_c2')
save('rec_c3.mat','rec_c3')
save('rec_c4.mat','rec_c4')
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
figure
histogram(residual0,'BinWidth',2)
set(gca,'FontSize',40);
xlabel('O-C Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
figure
plot(t0,residual0,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('O-C Residual/ns','fontsize',40);

figure
histogram(t0,'BinWidth',0.5)
set(gca,'FontSize',40);
xlabel('O-C Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
figure
plot(t0,residual0,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('O-C Residual/ns','fontsize',40);

figure
plot(t0-min(t0),residual0,'.','MarkerSize',12);
ylim([-100 100])
set(gca,'FontSize',40);
xlabel('Time [s]','fontsize',40);
ylabel('O-C Residual [ns]','fontsize',40);
start_hour = floor(min(t0)/3600);
start_second = floor((min(t0) - 3600 * start_hour)/60);
end_hour = floor(max(t0)/3600);
end_second = floor((max(t0) - 3600 * start_hour)/60);
% fprintf('开始于%2.0f时%2.0f分，结束于%2.0f时%2.0f分(UTC)\n',start_hour, start_second, end_hour, end_second);
%%

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
hold on
plot(t2,residual2,'.b','MarkerSize',12);
y=p*(t1-min(t0))+q;
% hold on
% plot(t1,y,'m')
% hold on
% plot(t1,y+oo,'m')
% set(gca,'FontSize',40);
% xlabel('Time/s','fontsize',40);
% ylabel('O-C Residual/ns','fontsize',40);

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

BW = 1e-2;%Binwidth:50ps
figure
histogram(residual4,'BinWidth',BW)
set(gca,'FontSize',40);
xlabel('O-C Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
num4=round((max(residual4)-min(residual4))/(BW));
N=length(residual4);%photon counts
x = min(residual4):(max(residual4)-min(residual4))/(num4-1):max(residual4);
y = hist(residual4,num4);
xx = x(:); %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w4 = gfit.c1/sqrt(2); % round trip sigma :ns
xc4 = gfit.b1;% 系统延时期望值
%计算半峰全宽（FWHM）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc4)/gfit.c1).^2);
% hold on
% figure
% plot(xxx,yyy)
FWHM4 = w4*2.355;
noise = (length(residual0)-length(residual4))/(2000*(max(t0)-min(t0)))*(max(residual4)-min(residual4))*(max(t3)-min(t3));%噪声水平
format long g
y=[w4*1e3,FWHM4*1e3,round(length(residual4)-noise)]; %双程标准差：ps； FWHM:ps； 扣除噪声后的有效回波数：counts
round(length(residual4)-noise)/(max(t4)-(min(t4)))
fprintf('开始于%2.0f时%2.0f分，结束于%2.0f时%2.0f分(UTC)\n',start_hour, start_second, end_hour, end_second);
fprintf('共计光子数为%3.0f，扣除噪声水平%3.0f，有效回波率为%1.3fcps\n',length(residual4),noise,round(length(residual4)-noise)/(max(t4)-(min(t4))))


    

%% 生成frd文本，格式需要修改
load('pre_data')
file_name = dir();
folder = file_name.folder;
pre_utc = pre_data(:,9)*24*3600;
pre_azi = pre_data(:,4);
pre_ele = pre_data(:,5);
azi = lagint(pre_utc,pre_azi,t3,10)';
ele = lagint(pre_utc,pre_ele,t3,10)';
B = [t3 azi ele residual3];
target_name = dir("*.a15").name;
if length(target_name) == 0
    target_name = dir("*.l21").name
    if length(target_name) == 0
        target_name = dir("*.l17").name
        if length(target_name) == 0
            target_name = dir("*.a14").name
            if length(target_name) == 0
                target_name = dir("*.a11").name
            end
        end
    end
end
target = target_name(end-2:end)
% if length(data_name) == 0
%     data_name = dir("*.l1*") 
% end
%     endsWith(data_name.name,".a11")
save_fold_name = ["e:\LLR_zhuh" + "(" + folder(14:21) + ")"];
if ~exist(save_fold_name,'dir')
    mkdir(save_fold_name);
end
fid=fopen([save_fold_name+'\'+target+'_azi_ele_residual.txt'],'w');%写入文件路径

[r,c]=size(B);            % 得到矩阵的行数和列数
for i=1:r
    for j=1:c
        fprintf(fid,'%5.12f\t',B(i,j));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);

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
        delta_t = 0 + coeff(1);
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
end

std(o_minus_c)

%% 测站偏心修正
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

%% 其他站检核
% Example 1: YARL
% importdata('D:\npt\la1\lageos1_20191116.npt', ' ', 12,0)
c=299792458;% m/s
std1 = 0;
data = [11 21308.400550200000     0.047387701325 std1 2  120.0     44   35.0   0.088  -0.083      -1.0   7.33 0                                                                
11 21399.000549900000     0.048710360126 std1 2  120.0     24   34.0   0.024  -0.490      -1.0   4.00 0                                                            
11 21496.400550999999     0.050220552225 std1 2  120.0      1   35.0   0.000  -3.000      -1.0   0.17 0                                                             
11 21602.400549800001     0.051954488686 std1 2  120.0      2    9.0   0.000  -2.000      -1.0   0.33 0];
gnt=[data(:,2) data(:,3) data(:,7)*c*1e-12/2];
pthdata = [20 21308.401  981.30 310.40  23. 0  
20 21399.001  981.30 310.50  23. 0    
20 21496.401  981.30 310.50  23. 0  
20 21602.401  981.30 310.00  23. 0];
wavelength = 0.532; %um
latitude = -29.0464; % °
height =  244e-3; % km

pth_time = pthdata(:,2);
baroPressure_list = pthdata(:,3) * 1e2; % Pa
temperature_list = pthdata(:,4)-273.15; % ℃
RH_list = pthdata(:,5); %

% wavelength = 1.064; %um
% latitude = 22.34639411111111; % °
% height = 0.405713; % km
data_and_reflector=dir('*.la1');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.la2');
end
temp = data_and_reflector.name;
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));

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
% coeff = [0 0 ];
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
%     coordinate_station_itrs0=[-2830744.471,4676580.282,3275072.816];
    coordinate_station_itrs0=[-2389008,5043332,-3078526];
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
        baroPressure = interp1(pth_time,baroPressure_list,gnt(i,1));
        temperature = interp1(pth_time,temperature_list,gnt(i,1));
        RH = interp1(pth_time,RH_list,gnt(i,1));
        correction = atmospheric_refraction(altc, baroPressure, temperature);
        elevationAngle = altc + correction;
        
        atm_delay = MPAtmDelayModelCal(wavelength,latitude,height,baroPressure,temperature,RH,elevationAngle);
        
        gravity_delay = gra_delay(YR,MONTH,DATE,hour(i),minu(i),sec(i),sp3_IntpXYZ_ITRS);
        
        com_corr = 0.251;
        
        t_up0 = (sp3R + atm_delay + gravity_delay - com_corr)/c;
        t_up = t_up0;
        if (abs(t_up - t_up0) > 1e-9)
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
        if delta_t_down<1e-9
            break
        end
    end
    
    t_up_down = [t_up_down;t_down,t_up];
    resi = (gnt(i,2) - (t_down + t_up)) * c;
    o_minus_c = [o_minus_c;resi];
    pos_sat_bound = [pos_sat_bound;sp3_IntpXYZ_GCRS'];
end

std(o_minus_c)

