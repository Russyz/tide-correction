%% 读数据并存储
tic
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
% channel4 = [];
% channel40 = [];
% channel41 =[];
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

% 主回波匹配
rec_c1=[];
rec_c2=[];
rec_c3=[];
rec_c4=[];
w=[];
% 预报插值
c=299792458;% m/s
% 远地靶
pre_TOF=2*2720.801/c;
% 镜筒地靶 
% pre_TOF=2*170/c;
% 近地靶
% pre_TOF=2*(1.136+0.51473)/c;
t_pre_receive = emission_time + pre_TOF;
% 匹配在阈值内残差
for m=1:length(emission_time)
    t1=[];
    w1=find(abs(channel1-t_pre_receive(m))<=1e-6);     
    if ~isempty(w1)
        residual1=channel10(w1)-emission_time0(m)+channel11(w1)-emission_time1(m)-pre_TOF; %残差s
        TOF1=residual1 + pre_TOF; %飞行时间s      
        for n=1:length(w1)
            rec_c1=[rec_c1;emission_time(m) TOF1(n) pre_TOF residual1(n)*1e9];
        end
    end
    w2=find(abs(channel2-t_pre_receive(m))<=1e-6);
    if ~isempty(w2)
        residual2=(channel20(w2)-emission_time0(m)+channel21(w2)-emission_time1(m)-pre_TOF);% 残差ns
        TOF2=residual2 + pre_TOF;% 飞行时间s
        for n=1:length(w2)
            rec_c2=[rec_c2;emission_time(m) TOF2(n) pre_TOF residual2(n)*1e9];
        end
    end
    w3=find(abs(channel3-t_pre_receive(m))<=1e-6);
    if ~isempty(w3)
        residual3=channel30(w3)-emission_time0(m)+channel31(w3)-emission_time1(m)-pre_TOF;% 残差s
        TOF3=residual3+pre_TOF;% 飞行时间s
        for n=1:length(w3)
            rec_c3=[rec_c3;emission_time(m) TOF3(n) pre_TOF residual3(n)*1e9];
        end
    end
    w4=find(abs(channel4-t_pre_receive(m))<=1e-6);
    if ~isempty(w4)
        t4=emission_time(m);% 日积秒 s
        residual4=channel40(w4)-emission_time0(m)+channel41(w4)-emission_time1(m)-pre_TOF;% 残差ns
        TOF4=residual4 + pre_TOF;% 飞行时间s
        for n=1:length(w4)
            rec_c4=[rec_c4;emission_time(m) TOF4(n) pre_TOF residual4(n)*1e9];
        end
    end
end
save('rec_c1.mat','rec_c1')
save('rec_c2.mat','rec_c2')
save('rec_c3.mat','rec_c3')
save('rec_c4.mat','rec_c4')
t_c1=rec_c1(:,1);

residual_c1=rec_c1(:,4);
t_c2=rec_c2(:,1);

residual_c2=rec_c2(:,4);
t_c3=rec_c3(:,1);

residual_c3=rec_c3(:,4);
t_c4=rec_c4(:,1);

residual_c4=rec_c4(:,4);
t0=[t_c1;t_c2;t_c3;t_c4];
residual0=[residual_c1;residual_c2;residual_c3;residual_c4];
figure
subplot(1,2,1);
plot(t0-min(t0),residual0,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time [s]','fontsize',40);
ylabel('O-C Residual [ns]','fontsize',40);
%set(gcf,'outerposition',get(0,'screensize'));	% 设置图形窗口位置和外尺寸为屏幕大小
%saveas(gcf,'residual.jpg');
subplot(1,2,2);
BW = 5;
histogram(residual0,'BinWidth',BW,'Orientation','horizontal') 
set(gca,'FontSize',40);
xlabel('O-C Residual [ns]','fontsize',40);
ylabel('Photon counts','fontsize',40);
set(gcf,'outerposition',get(0,'screensize'));	% 设置图形窗口位置和外尺寸为屏幕大小
saveas(gcf,'Residual&historgram.jpg');
close all

%% 扣除大气延时并计算系统延时
% 第一通道
clc;clear all;close all
load('rec_c1.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('F:\matlab')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt=(2720.801*ng-2720.801)/c*10^9; %3km地靶大气延时ns
t_c1=rec_c1(:,1);residual_c1=rec_c1(:,4)-2*delay_gt';ToF_c1=rec_c1(:,2)-2*delay_gt'*10^-9;
figure
plot(t_c1,rec_c1(:,4),'.')
oo = 1.5;   %设定离散程度ns
[t01 p1] = sort([t_c1]); %对记录时刻进行排序
residual01 = residual_c1(p1);
TOF01 = ToF_c1(p1);
w = 0;
for i = 0
    for j=240:0.1:270
        y=i*(t01-min(t01))+j;
        w1 = [];
        w1 = find(residual01>y&residual01<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %斜率
            q = j; %截距
            w2 = w1; %记录位置信息
        end
    end
end
t1 = t01(w2);
residual1 = residual01(w2);
TOF1 = TOF01(w2);
figure
plot(t1,residual1,'.')
% 2.5σ滤波
for k=1:10
    residual30=[];
    [p2,S2]=polyfit(t1-min(t1),residual1,1);
    [y_fit2,delta2]=polyval(p2,t1-min(t1),S2);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    t2=t1(w3);
    residual2=residual1(w3);
    TOF2=TOF1(w3);
end
t3=t2;
residual3=residual2;
figure
plot(t3,residual3,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('Residual/ns','fontsize',40);
figure
histogram(residual3)
TOF3=TOF2;
f_c1=[t3 TOF3 residual3];
save('f_c1.mat','f_c1')
P_F1=fitdist(abs(residual3),'Poisson');
mu_channel1=P_F1.lambda;
save('mu_channel1.mat','mu_channel1') %保存第一通道系统延时
std(residual3-mu_channel1)*15
%
% 第二通道
clc;clear all;close all
load('rec_c2.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('F:\matlab')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt=(2720.801*ng-2720.801)/c*10^9; %3km地靶大气延时ns
t_c1=rec_c2(:,1);residual_c1=rec_c2(:,4)-2*delay_gt';ToF_c1=rec_c2(:,2)-2*delay_gt'*10^-9;
oo = 1.5;   %设定离散程度ns
[t01 p1] = sort([t_c1]); %对记录时刻进行排序
residual01 = residual_c1(p1);
TOF01 = ToF_c1(p1);
w = 0;
for i = -0.001:0.0001:0.001
    for j=240:0.1:260
        y=i*(t01-min(t01))+j;
        w1 = [];
        w1 = find(residual01>y&residual01<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %斜率
            q = j; %截距
            w2 = w1; %记录位置信息
        end
    end
end
t1 = t01(w2);
residual1 = residual01(w2);
TOF1 = TOF01(w2);
% 2.5σ滤波
for k=1:10
    residual30=[];
    [p2,S2]=polyfit(t1-min(t1),residual1,1);
    [y_fit2,delta2]=polyval(p2,t1-min(t1),S2);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    t2=t1(w3);
    residual2=residual1(w3);
    TOF2=TOF1(w3);
end
t3=t2;
residual3=residual2;
figure
plot(t3,residual3,'.')
TOF3=TOF2;
f_c2=[t3 TOF3 residual3];
save('f_c2.mat','f_c2')
P_F1=fitdist(abs(residual3),'Poisson');
mu_channel2=P_F1.lambda;
save('mu_channel2.mat','mu_channel2') %保存第2通道系统延时

% 第三通道
clc;clear all;close all
load('rec_c3.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('F:\matlab')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt=(2720.801*ng-2720.801)/c*10^9; %3km地靶大气延时ns
t_c1=rec_c3(:,1);residual_c1=rec_c3(:,4)-2*delay_gt';ToF_c1=rec_c3(:,2)-2*delay_gt'*10^-9;
oo = 1.5;   %设定离散程度ns
[t01 p1] = sort([t_c1]); %对记录时刻进行排序
residual01 = residual_c1(p1);
TOF01 = ToF_c1(p1);
w = 0;
for i = -0.001:0.0001:0.001
    for j=240:0.1:260
        y=i*(t01-min(t01))+j;
        w1 = [];
        w1 = find(residual01>y&residual01<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %斜率
            q = j; %截距
            w2 = w1; %记录位置信息
        end
    end
end
t1 = t01(w2);
residual1 = residual01(w2);
TOF1 = TOF01(w2);
% 2.5σ滤波
for k=1:10
    residual30=[];
    [p2,S2]=polyfit(t1-min(t1),residual1,1);
    [y_fit2,delta2]=polyval(p2,t1-min(t1),S2);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    t2=t1(w3);
    residual2=residual1(w3);
    TOF2=TOF1(w3);
end
t3=t2;
residual3=residual2;
TOF3=TOF2;
f_c3=[t3 TOF3 residual3];
save('f_c3.mat','f_c3')
P_F1=fitdist(abs(residual3),'Poisson');
mu_channel3=P_F1.lambda;
save('mu_channel3.mat','mu_channel3') %保存第3通道系统延时

% 第四通道
clc;clear all;close all
load('rec_c4.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('F:\matlab')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt=(2720.801*ng-2720.801)/c*10^9; %3km地靶大气延时ns
t_c1=rec_c4(:,1);residual_c1=rec_c4(:,4)-2*delay_gt';ToF_c1=rec_c4(:,2)-2*delay_gt'*10^-9;
oo = 1.5;   %设定离散程度ns
[t01 p1] = sort([t_c1]); %对记录时刻进行排序
residual01 = residual_c1(p1);
TOF01 = ToF_c1(p1);
w = 0;
for i = -0.001:0.0001:0.001
    for j=240:0.1:260
        y=i*(t01-min(t01))+j;
        w1 = [];
        w1 = find(residual01>y&residual01<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %斜率
            q = j; %截距
            w2 = w1; %记录位置信息
        end
    end
end
t1 = t01(w2);
residual1 = residual01(w2);
TOF1 = TOF01(w2);
% 2.5σ滤波
for k=1:10
    residual30=[];
    [p2,S2]=polyfit(t1-min(t1),residual1,1);
    [y_fit2,delta2]=polyval(p2,t1-min(t1),S2);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    t2=t1(w3);
    residual2=residual1(w3);
    TOF2=TOF1(w3);
end
t3=t2;
residual3=residual2;
TOF3=TOF2;
f_c4=[t3 TOF3 residual3];
save('f_c4.mat','f_c4')
P_F1=fitdist(abs(residual3),'Poisson');
mu_channel4=P_F1.lambda;
save('mu_channel4.mat','mu_channel4') %保存第4通道系统延时
mean(delta2)*15
%% 
clc;clear all;close all

load('f_c1.mat')
load('f_c2.mat')
load('f_c3.mat')
load('f_c4.mat')
residual1 = f_c1(:,3);
residual2 = f_c2(:,3);
residual3 = f_c3(:,3);
residual4 = f_c4(:,3);
BW = 2e-2;%Binwidth:50ps
figure
subplot(2,2,1);
histogram(residual1,'BinWidth',BW) 
set(gca,'FontSize',20);
xlabel('Residual/ns','fontsize',20);
ylabel('Photon counts','fontsize',20);
num1=round((max(residual1)-min(residual1))/(BW));
N=length(residual1);%photon counts
x = min(residual1):(max(residual1)-min(residual1))/(num1-1):max(residual1);
y = hist(residual1,num1);
xx = x(:); %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w1 = gfit.c1/sqrt(2) % round trip sigma :ns
xc1 = gfit.b1;% 系统延时期望值
%计算半峰全宽（FWHM）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc1)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FWHM1 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));


subplot(2,2,2);
histogram(residual2,'BinWidth',BW) 
set(gca,'FontSize',20);
xlabel('Residual/ns','fontsize',20);
ylabel('Photon counts','fontsize',20);
num2=round((max(residual2)-min(residual2))/(BW));
N=length(residual2);%photon counts
x = min(residual2):(max(residual2)-min(residual2))/(num2-1):max(residual2);
y = hist(residual2,num2);
xx = x(:); %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w2 = gfit.c1/sqrt(2) % round trip sigma :ns
xc2 = gfit.b1;% 系统延时期望值
%计算半峰全宽（FWHM）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc2)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FWHM2 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));

subplot(2,2,3);
histogram(residual3,'BinWidth',BW) 
set(gca,'FontSize',20);
xlabel('Residual/ns','fontsize',20);
ylabel('Photon counts','fontsize',20);
num3=round((max(residual3)-min(residual3))/(BW));
N=length(residual3);%photon counts
x = min(residual3):(max(residual3)-min(residual3))/(num3-1):max(residual3);
y = hist(residual3,num3);
xx = x(:); %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w3 = gfit.c1/sqrt(2) % round trip sigma :ns
xc3 = gfit.b1;% 系统延时期望值
%计算半峰全宽（FWHM）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc3)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FWHM3 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));

subplot(2,2,4);
histogram(residual4,'BinWidth',BW) 
set(gca,'FontSize',20);
xlabel('Residual/ns','fontsize',20);
ylabel('Photon counts','fontsize',20);
num4=round((max(residual4)-min(residual4))/(BW));
N=length(residual4);%photon counts
x = min(residual4):(max(residual4)-min(residual4))/(num4-1):max(residual4);
y = hist(residual4,num4);
xx = x(:); %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w4 = gfit.c1/sqrt(2) % round trip sigma :ns
xc4 = gfit.b1;% 系统延时期望值
%计算半峰全宽（FWHM）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc4)/gfit.c1).^2);
hold on
plot(xxx,yyy)
saveas(gcf,'save.jpg');
FWHM4 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));

xc = [xc1,xc2,xc3,xc4];
w = [w1,w2,w3,w4]*1e3;
FWHM = [FWHM1,FWHM2,FWHM3,FWHM4]*1e3;
n_counts = [length(residual1),length(residual2),length(residual3),length(residual4)];
a = [xc;w;w*2.355;n_counts]
%%
close all
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
    residual_c1=rec_c1(:,4);
end
% 判断第2通道是否测量到在阈值范围的数据并赋值
if length(rec_c2)==0
    t_c2=[];TOF_c2=[];residual_c2=[]
else
    t_c2=rec_c2(:,1);
    TOF_c2=rec_c2(:,2)-mu_channel2*10^-9;
    residual_c2=rec_c2(:,4);
end
% 判断第3通道是否测量到在阈值范围的数据并赋值
if length(rec_c3)==0
    t_c3=[];TOF_c3=[];residual_c3=[];
else
    t_c3=rec_c3(:,1);
    TOF_c3=rec_c3(:,2)-mu_channel3*10^-9;
    residual_c3=rec_c3(:,4);
end
% 判断第4通道是否测量到在阈值范围的数据并赋值
if length(rec_c4)==0
    t_c4=[];TOF_c4=[];residual_c4=[];
else
    t_c4=rec_c4(:,1);
    TOF_c4=rec_c4(:,2)-mu_channel4*10^-9;
    residual_c4=rec_c4(:,4);
end
% 扣除系统延时得到TOF和O-C残差
t0=[t_c1;t_c2;t_c3;t_c4];
TOF0=[TOF_c1;TOF_c2;TOF_c3;TOF_c4];
residual0=[residual_c1;residual_c2;residual_c3;residual_c4];
t0 = t0(find(residual0>242&residual0<254));
residual0 = residual0(find(residual0>242&residual0<254));
figure
subplot(1,2,1);
plot(t0-min(t0),residual0,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('time [s]','fontsize',40);
ylabel('O-C residual [ns]','fontsize',40);
subplot(1,2,2);

BW = 3e-2;%Binwidth:30ps


histogram(residual0,'BinWidth',BW,'Orientation','horizontal') 

%ylim = ([242 254]);
set(gca,'FontSize',40);
xlabel('Photon counts','fontsize',40);
ylabel('O-C Residual [ns]','fontsize',40);
set(gcf,'outerposition',get(0,'screensize'));	% 设置图形窗口位置和外尺寸为屏幕大小
saveas(gcf,'Residual&historgram.jpg');