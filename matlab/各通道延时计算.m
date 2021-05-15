%% 读取主波时刻和各通道回波时刻值并保存
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

%% 各通道匹配 记录数据在rec_c1、rec_c2、rec_c3、rec_c4里面，每一个包含测量时刻(s)、实测TOF(s)、预报插值TOF(s)、残差(ns)
tic
clc;clear all
load('emission_time.mat')
load('channel1.mat')
c=299792458;%光速
tof1=2*2720.801/c;
rec_c1=[];
for m=1:length(emission_time)
    w1=find(emission_time(m)+tof1-1e-6<channel1&channel1<emission_time(m)+tof1+1e-6);
    if ~isempty(w1)
        t1=emission_time(m);% 日积秒 s
        residual1=(channel1(w1)-emission_time(m)-tof1)*1e9;% 残差ns
        ToF1=channel1(w1)-emission_time(m);% 飞行时间s
        for n=1:length(w1)
            rec_c1=[rec_c1;t1 ToF1(n) tof1 residual1(n)];
        end
    end    
end
save('rec_c1.mat','rec_c1')
toc
% 第二通道匹配
clc;clear all
load('emission_time.mat')
load('channel2.mat')
c=299792458;%光速
tof2=2*2720.801/c;
rec_c2=[];
for m=1:length(emission_time)
    w2=find(emission_time(m)+tof2-1e-6<channel2&channel2<emission_time(m)+tof2+1e-6);
    if ~isempty(w2)
        t2=emission_time(m);% 日积秒 s
        residual2=(channel2(w2)-emission_time(m)-tof2)*1e9;% 残差ns
        ToF2=channel2(w2)-emission_time(m);% 飞行时间s
        for n=1:length(w2)
            rec_c2=[rec_c2;t2 ToF2(n) tof2 residual2(n)];
        end
    end    
end
save('rec_c2.mat','rec_c2')
% 第三通道匹配
clc;clear all
load('emission_time.mat')
load('channel3.mat')
c=299792458;%光速
tof3=2*2720.801/c;
rec_c3=[];
for m=1:length(emission_time)
    w3=find(emission_time(m)+tof3-1e-6<channel3&channel3<emission_time(m)+tof3+1e-6);
    if ~isempty(w3)
        t3=emission_time(m);% 日积秒 s
        residual3=(channel3(w3)-emission_time(m)-tof3)*1e9;% 残差ns
        ToF3=channel3(w3)-emission_time(m);% 飞行时间s
        for n=1:length(w3)
            rec_c3=[rec_c3;t3 ToF3(n) tof3 residual3(n)];
        end
    end
end
save('rec_c3.mat','rec_c3')
% 第四通道匹配
clc;clear all
load('emission_time.mat')
load('channel4.mat')
c=299792458;%光速
tof4=2*2720.801/c;
rec_c4=[];
for m=1:length(emission_time)
    w4=find(emission_time(m)+tof4-1e-6<channel4&channel4<emission_time(m)+tof4+1e-6);
    if ~isempty(w4)
        t4=emission_time(m);% 日积秒 s
        residual4=(channel4(w4)-emission_time(m)-tof4)*1e9;% 残差ns
        ToF4=channel4(w4)-emission_time(m);% 飞行时间s
        for n=1:length(w4)
            rec_c4=[rec_c4;t4 ToF4(n) tof4 residual4(n)];
        end
    end
end
save('rec_c4.mat','rec_c4')
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
addpath(genpath('C:\Users\thinkpad\Desktop\数据处理\python')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt=(2720.801*ng-2720.801)/c*10^9; %3km地靶大气延时ns
t_c1=rec_c1(:,1);residual_c1=rec_c1(:,4)-2*delay_gt';ToF_c1=rec_c1(:,2)-2*delay_gt'*10^-9;
oo = 1.5;   %设定离散程度ns
[t01 p1] = sort([t_c1]); %对记录时刻进行排序
residual01 = residual_c1(p1);
TOF01 = ToF_c1(p1);
w = 0;
for i = -0.001:0.0001:0.001
    for j=-240:0.1:250
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
    
%%
% 第二通道
clc;clear all;close all
load('rec_c2.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('C:\Users\thinkpad\Desktop\数据处理\python')) %大气折射率计算公式路径
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
    for j=-240:0.1:250
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
addpath(genpath('C:\Users\thinkpad\Desktop\数据处理\python')) %大气折射率计算公式路径
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
    for j=-240:0.1:250
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

% 第二通道
clc;clear all;close all
load('rec_c4.mat')
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('C:\Users\thinkpad\Desktop\数据处理\python')) %大气折射率计算公式路径
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
    for j=-240:0.1:250
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
BW = 5e-2;%Binwidth:50ps
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

BW = 5e-2;%Binwidth:50ps
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
BW = 5e-2;%Binwidth:50ps
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
BW = 5e-2;%Binwidth:50ps
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
FWHM4 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));

xc = [xc1,xc2,xc3,xc4];
w = [w1,w2,w3,w4]*1e3;
FWHM = [FWHM1,FWHM2,FWHM3,FWHM4]*1e3;
n_counts = [length(residual1),length(residual2),length(residual3),length(residual4)];
a = [xc;w;FWHM;n_counts]
