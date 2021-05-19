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
    t_c2=[];TOF_c2=[];residual_c2=[];
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


%% %% 算出标准点NP
load('pre_data.mat')
c=299792458;%光速
sl=100;%标准点步长选取
gnt=[];
gnd=[];

% NP_xy=[];
for i=1:ceil((max(t3)-min(t3))/sl)
    w5=find((min(t3)+sl*(i-1))<=t3&t3<(min(t3)+sl*i));
    if length(w5) == 0
        break
    else
        t5=t3(w5);
        residual5=residual3(w5);
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
%         addpath('/home/z8/文档/潮汐修正/test file')
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
data_and_reflector=dir('*.a15');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.l21');
        if length(data_and_reflector) == 0
            data_and_reflector=dir('*.l17');
            if length(data_and_reflector) == 0
                data_and_reflector=dir('*.a11');
                    if length(data_and_reflector) == 0
                        data_and_reflector=dir('*.a14');
                    end
            end
        end
end


temp = data_and_reflector.name;
target = temp(end-2:end);
YR = str2num(temp(1:4));
MONTH = str2num(temp(5:6));
DATE = str2num(temp(7:8));
% coordinate_station_itrs=[-2358691.210,5410611.484,2410087.607];
o_minus_c = [];
% if coeff
% coeff = [0 0 ];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    tof = lunar_target_cpf(target,YR,MONTH,DATE,hour(i),minu(i),sec(i));
    
    resi = (gnt(i,2) - tof) * c;
    o_minus_c = [o_minus_c;resi];
    
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
mean(o_minus_c/2)
%% 测站偏心修正
% Load P.T.H. data
t_bias = coeff(1);
o_minus_c = [];
for i = 1:length(gnt(:,1))
    %     syms delta_x delta_y delta_z
    hour(i) = floor(gnt(i,1)/3600);
    minu(i) = floor((gnt(i,1) - hour(i) * 3600)/60);
    sec(i) = gnt(i,1) - hour(i) * 3600 - minu(i) * 60;
    tof = lunar_target_cpf(target,YR,MONTH,DATE,hour(i),minu(i),sec(i)+t_bias);
    
    resi = (gnt(i,2) - tof) * c;
    o_minus_c = [o_minus_c;resi];
    
end

std(o_minus_c)


