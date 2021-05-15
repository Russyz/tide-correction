%% 读取数据
clc;clear all;close all
data = importdata('C:\Users\thinkpad\Documents\WeChat Files\zcs1103\FileStorage\File\2020-06\1106地靶数据.dat');
residual = data.data; %各通道匹配残差
c = 299792458; % 光速m/s
gt_distance = 2720.801; % 地靶与参考点之间距离m
tof_gt = 2 * gt_distance / c; % 地靶到参考点往返时间s

residual_c01 = residual(:,1);
w1 = find(~isnan(residual(:,1)));
residual_c01 = residual_c01(w1) - tof_gt; %第一通道匹配残差

residual_c02 = residual(:,2);
w2 = find(~isnan(residual(:,2)));
residual_c02 = residual_c02(w2) - tof_gt; %第二通道匹配残差

residual_c03 = residual(:,3);
w3 = find(~isnan(residual(:,3)));
residual_c03 = residual_c03(w3) - tof_gt; %第三通道匹配残差

residual_c04 = residual(:,4);
w4 = find(~isnan(residual(:,4)));
residual_c04 = residual_c04(w4) - tof_gt; %第四通道匹配残差

%% 扣除大气延时
c=299792458;%光速
pth_data = dir('*.pth'); 
pth=importdata(pth_data.name) ;
baroPressure=pth(2);%大气压mbar
temperature=pth(3);%温度℃
RH=pth(4);%相对湿度具体值%
addpath(genpath('C:\Users\thinkpad\Desktop\数据处理\python')) %大气折射率计算公式路径
wavelength = 1.064;%um
ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH); %大气折射率
delay_gt = gt_distance*(ng-1) / c; %3km地靶大气延时ns
% 各通道扣除大气延时后，系统延时
sys_delay_c01 = residual_c01 - 2 * delay_gt;
sys_delay_c02 = residual_c02 - 2 * delay_gt;
sys_delay_c03 = residual_c03 - 2 * delay_gt;
sys_delay_c04 = residual_c04 - 2 * delay_gt;

%% 第一通道延时
w = 0;
oo = 1.5e-9; %ns
for j=2.4e-7:1e-10:2.5e-7
    w1 = [];
    w1 = find(sys_delay_c01>j&sys_delay_c01<(j+oo));
    ww = length(w1);
    if w<ww
        w = ww;
        q = j; %截距
        w2 = w1; %记录位置信息
    end
end
residual1 = sys_delay_c01(w2);

% 2.5σ滤波10次
for k=1:10
    residual30=[];
    [p1,S1]=polyfit((1:length(residual1))',residual1,1);
    [y_fit2,delta2]=polyval(p1,1:length(residual1),S1);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    residual2=residual1(w3);
end
% 滤波后残差图
figure
plot(residual2,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time','fontsize',40);
ylabel('Residual/s','fontsize',40);

% 计算系统延时，画残差柱状图
BW = 5e-2;%Binwidth:50ps
figure
histogram(residual2*1e9,'BinWidth',BW) 
set(gca,'FontSize',40);
xlabel('Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
% 高斯拟合
num=round((max(residual2)-min(residual2))/(BW*1e-9));
N=length(residual2);%photon counts
x = min(residual2):(max(residual2)-min(residual2))/(num-1):max(residual2);
y = hist(residual2,num);
xx = x(:)*1e9; %残差ns
yy = y(:);
% 高斯拟合
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w = gfit.c1/sqrt(2) % round trip sigma :ns
xc = gfit.b1;% 系统延时期望值
%计算半峰全宽（FMHW）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FMHW = xxx(max(find(abs(yyy-max(yyy)/2)<0.1)))-xxx(min(find(abs(yyy -max(yyy)/2)<0.1)));

P_F1=fitdist(abs(residual2),'Poisson');
mu_channel1=P_F1.lambda; %系统延时s
% save('mu_channel1.mat','mu_channel1') %保存第1通道系统延时
%% 第二通道延时
w = 0;
oo = 1.5e-9; % 残差门宽选取 ns
for j=2.4e-7:1e-10:2.5e-7
    w1 = [];
    w1 = find(sys_delay_c02>j&sys_delay_c02<(j+oo));
    ww = length(w1);
    if w<ww
        w = ww;
        q = j; %截距
        w2 = w1; %记录位置信息
    end
end
residual1 = sys_delay_c02(w2); %扫描得到最多残差数据

% 2.5σ滤波10次
for k=1:10
    residual30=[];
    [p1,S1]=polyfit((1:length(residual1))',residual1,1);
    [y_fit2,delta2]=polyval(p1,1:length(residual1),S1);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    residual2=residual1(w3);
end
% 滤波后残差图
figure
plot(residual2,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time','fontsize',40);
ylabel('Residual/s','fontsize',40);
% 计算系统延时，画残差柱状图
BW = 5e-2;%Binwidth:50ps
figure
histogram(residual2*1e9,'BinWidth',BW) 
set(gca,'FontSize',40);
xlabel('Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
% 高斯拟合
num=round((max(residual2)-min(residual2))/(BW*1e-9));
N=length(residual2);%photon counts
x = min(residual2):(max(residual2)-min(residual2))/(num-1):max(residual2);
y = hist(residual2,num);
xx = x(:)*1e9; %残差ns
yy = y(:);
% 高斯拟合
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w = gfit.c1/sqrt(2) % round trip sigma :ns
xc = gfit.b1;% 系统延时期望值
%计算半峰全宽（FMHW）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FMHW = xxx(max(find(abs(yyy-max(yyy)/2)<0.1)))-xxx(min(find(abs(yyy -max(yyy)/2)<0.1)));
P_F1=fitdist(abs(residual2),'Poisson');
mu_channel2=P_F1.lambda;%系统延时s
% save('mu_channel2.mat','mu_channel2') %保存第2通道系统延时
%% 第三通道延时
w = 0;
oo = 1.5e-9; % 残差门宽选取 ns
for j=2.4e-7:1e-10:2.5e-7
    w1 = [];
    w1 = find(sys_delay_c03>j&sys_delay_c03<(j+oo));
    ww = length(w1);
    if w<ww
        w = ww;
        q = j; %截距
        w2 = w1; %记录位置信息
    end
end
residual1 = sys_delay_c03(w2); %扫描得到最多残差数据

% 2.5σ滤波10次
for k=1:10
    residual30=[];
    [p1,S1]=polyfit((1:length(residual1))',residual1,1);
    [y_fit2,delta2]=polyval(p1,1:length(residual1),S1);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    residual2=residual1(w3);
end
% 滤波后残差图
figure
plot(residual2,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time','fontsize',40);
ylabel('Residual/s','fontsize',40);
% 计算系统延时，画残差柱状图
BW = 5e-2;%Binwidth:50ps
figure
histogram(residual2*1e9,'BinWidth',BW) 
set(gca,'FontSize',40);
xlabel('Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
% 高斯拟合
num=round((max(residual2)-min(residual2))/(BW*1e-9));
N=length(residual2);%photon counts
x = min(residual2):(max(residual2)-min(residual2))/(num-1):max(residual2);
y = hist(residual2,num);
xx = x(:)*1e9; %残差ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w = gfit.c1/sqrt(2) % round trip sigma :ns
xc = gfit.b1;% 系统延时期望值
%计算半峰全宽（FMHW）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FMHW = xxx(max(find(abs(yyy-max(yyy)/2)<0.1)))-xxx(min(find(abs(yyy -max(yyy)/2)<0.1)));
P_F1=fitdist(abs(residual2),'Poisson');
mu_channel3=P_F1.lambda;%系统延时s
% save('mu_channel3.mat','mu_channel3') %保存第3通道系统延时
%% 第四通道延时
w = 0;
oo = 1.5e-9; % 残差门宽选取 ns
for j=2.4e-7:1e-10:2.5e-7
    w1 = [];
    w1 = find(sys_delay_c04>j&sys_delay_c04<(j+oo));
    ww = length(w1);
    if w<ww
        w = ww;
        q = j; %截距
        w2 = w1; %记录位置信息
    end
end
residual1 = sys_delay_c04(w2); %扫描得到最多残差数据

% 2.5σ滤波10次
for k=1:10
    residual30=[];
    [p1,S1]=polyfit((1:length(residual1))',residual1,1);
    [y_fit2,delta2]=polyval(p1,1:length(residual1),S1);
    for l=1:length(delta2)        
        if abs(residual1(l)-y_fit2(l))<2.5*delta2(l)
            residual30(l)=residual1(l);
        end
    end
    w3=find(residual30);
    residual2=residual1(w3);
end
% 滤波后残差图
figure
plot(residual2,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time','fontsize',40);
ylabel('Residual/s','fontsize',40);
% 计算系统延时，画残差柱状图
BW = 5e-2;%Binwidth:50ps
figure
histogram(residual2*1e9,'BinWidth',BW) 
set(gca,'FontSize',40);
xlabel('Residual/ns','fontsize',40);
ylabel('Photon counts','fontsize',40);
% 高斯拟合
num=round((max(residual2)-min(residual2))/(BW*1e-9));
N=length(residual2);%photon counts
x = min(residual2):(max(residual2)-min(residual2))/(num-1):max(residual2);
y = hist(residual2,num);
xx = x(:)*1e9; %残差ns
yy = y(:);
% 设定拟合
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w = gfit.c1/sqrt(2) % round trip sigma :ns
xc = gfit.b1;% 系统延时期望值
%计算半峰全宽（FMHW）
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc)/gfit.c1).^2);
hold on
plot(xxx,yyy)
FMHW = xxx(max(find(abs(yyy-max(yyy)/2)<0.1)))-xxx(min(find(abs(yyy -max(yyy)/2)<0.1)));
P_F1=fitdist(abs(residual2),'Poisson');
mu_channel4=P_F1.lambda;%系统延时s
% save('mu_channel4.mat','mu_channel4') %保存第4通道系统延时
