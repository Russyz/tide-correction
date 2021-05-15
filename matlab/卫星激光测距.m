%% �����ݲ��洢
clc;clear all
predict = dir('*.a15'); %�޸�Ԥ�����Ǽ��
if length(predict) == 0
    predict = dir('*.l21');
end
if length(predict) == 0
    predict = dir('*.a11');
end    
if length(predict) == 0
    predict = dir('*.a14');
end
if length(predict) == 0
    predict = dir('*.l17');
end
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
% ��ͨ���ز�ʱ��ֵ
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
channel4 = [];
channel40 = [];
channel41 =[];
save('channel4.mat','channel4')
save('channel40.mat','channel40')
save('channel41.mat','channel41')
% ����ʱ��
emission_time=Te.data(:,1)+Te.data(:,2);
emission_time0=Te.data(:,1);
emission_time1=Te.data(:,2);
save('emission_time.mat','emission_time')
save('emission_time0.mat','emission_time0')
save('emission_time1.mat','emission_time1')
P=importdata(predict.name, ' ', 2,0);
pre_data=P.data;
save('pre_data.mat','pre_data')
%% ���ز�ƥ��
% % ��ȡԤ��
% load('pre_data.mat')
% % ��ȡ����
% load('emission_time.mat')
% load('emission_time0.mat')
% load('emission_time1.mat')
% % ��ȡ�ز�1
% load('channel1.mat')
% load('channel10.mat')
% load('channel11.mat')
% % ��ȡ�ز�2
% load('channel2.mat')
% load('channel20.mat')
% load('channel21.mat')
% % ��ȡ�ز�3
% load('channel3.mat')
% load('channel30.mat')
% load('channel31.mat')
% % ��ȡ�ز�4
% load('channel4.mat')
% load('channel40.mat')
% load('channel41.mat')

% Ԥ��ʱ���tof
tic
pre_utc=pre_data(:,9)*24*3600;
pre_tof=pre_data(:,10);
rec_c1=[];
rec_c2=[];
rec_c3=[];
rec_c4=[];
w=[];
% Ԥ����ֵ
pre_TOF=interp1(pre_utc,pre_tof,emission_time,'spline');
t_pre_receive = emission_time + pre_TOF;
% ƥ������ֵ�ڲв�
for m=1:length(emission_time)
    t1=[];
    w1=find(abs(channel1-t_pre_receive(m))<=1e-6);     
    if ~isempty(w1)
        residual1=channel10(w1)-emission_time0(m)+channel11(w1)-emission_time1(m)-pre_TOF(m); %�в�s
        TOF1=residual1 + pre_TOF(m); %����ʱ��s      
        for n=1:length(w1)
            rec_c1=[rec_c1;emission_time(m) TOF1(n) pre_TOF(m) residual1(n)*1e9];
        end
    end
    w2=find(abs(channel2-t_pre_receive(m))<=1e-6);
    if ~isempty(w2)
        residual2=(channel20(w2)-emission_time0(m)+channel21(w2)-emission_time1(m)-pre_TOF(m));% �в�ns
        TOF2=residual2 + pre_TOF(m);% ����ʱ��s
        for n=1:length(w2)
            rec_c2=[rec_c2;emission_time(m) TOF2(n) pre_TOF(m) residual2(n)*1e9];
        end
    end
    w3=find(abs(channel3-t_pre_receive(m))<=1e-6);
    if ~isempty(w3)
        residual3=channel30(w3)-emission_time0(m)+channel31(w3)-emission_time1(m)-pre_TOF(m);% �в�s
        TOF3=residual3+pre_TOF(m);% ����ʱ��s
        for n=1:length(w3)
            rec_c3=[rec_c3;emission_time(m) TOF3(n) pre_TOF(m) residual3(n)*1e9];
        end
    end
    w4=find(abs(channel4-t_pre_receive(m))<=1e-6);
    if ~isempty(w4)
        t4=emission_time(m);% �ջ��� s
        residual4=channel40(w4)-emission_time0(m)+channel41(w4)-emission_time1(m)-pre_TOF(m);% �в�ns
        TOF4=residual4 + pre_TOF(m);% ����ʱ��s
        for n=1:length(w4)
            rec_c4=[rec_c4;emission_time(m) TOF4(n) pre_TOF(m) residual4(n)*1e9];
        end
    end
end
save('rec_c1.mat','rec_c1')
save('rec_c2.mat','rec_c2')
save('rec_c3.mat','rec_c3')
save('rec_c4.mat','rec_c4')
toc
%% �۳�ϵͳ��ʱ���õ�����TOF������������ʱ�Ȼ�δ��������
close all
c=299792458;% m/s
load('rec_c1.mat')
load('rec_c2.mat')
load('rec_c3.mat')
load('rec_c4.mat')
% ̽������ͨ����ʱ���Ȳ���������ͨ�������ݣ�
% �Ե���ͨ���������ݾ��Ƚ���ͳ��
load('mu_channel1.mat')
load('mu_channel2.mat')
load('mu_channel3.mat')
load('mu_channel4.mat')
% �жϵ�1ͨ���Ƿ����������ֵ��Χ�����ݲ���ֵ
if length(rec_c1)==0
    t_c1=[];TOF_c1=[];residual_c1=[];
else
    t_c1=rec_c1(:,1);
    TOF_c1=rec_c1(:,2)-mu_channel1*10^-9;
    residual_c1=rec_c1(:,4)-mu_channel1;
end
% �жϵ�2ͨ���Ƿ����������ֵ��Χ�����ݲ���ֵ
if length(rec_c2)==0
    t_c2=[];TOF_c2=[];residual_c2=[]
else
    t_c2=rec_c2(:,1);
    TOF_c2=rec_c2(:,2)-mu_channel2*10^-9;
    residual_c2=rec_c2(:,4)-mu_channel2;
end
% �жϵ�3ͨ���Ƿ����������ֵ��Χ�����ݲ���ֵ
if length(rec_c3)==0
    t_c3=[];TOF_c3=[];residual_c3=[];
else
    t_c3=rec_c3(:,1);
    TOF_c3=rec_c3(:,2)-mu_channel3*10^-9;
    residual_c3=rec_c3(:,4)-mu_channel3;
end
% �жϵ�4ͨ���Ƿ����������ֵ��Χ�����ݲ���ֵ
if length(rec_c4)==0
    t_c4=[];TOF_c4=[];residual_c4=[];
else
    t_c4=rec_c4(:,1);
    TOF_c4=rec_c4(:,2)-mu_channel4*10^-9;
    residual_c4=rec_c4(:,4)-mu_channel4;
end
% �۳�ϵͳ��ʱ�õ�TOF��O-C�в�
t0=[t_c1;t_c2;t_c3;t_c4];
TOF0=[TOF_c1;TOF_c2;TOF_c3;TOF_c4];
residual0=[residual_c1;residual_c2;residual_c3;residual_c4];
figure
plot(t0-min(t0),residual0,'.','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('time [hour]','fontsize',40);
ylabel('O-C residual [ns]','fontsize',40);
% ���źŽ���ʶ��
oo = 2;   %�趨������Ĵ�С��ns
[t1 p1] = sort([t0]); %�Լ�¼ʱ�̽�������
% pos = find(t1>4.7e4);
% t1 = t1(pos);
residual1 = (residual0(p1));
% residual1 = residual1(pos);
TOF1 = TOF0(p1);
% www = find((11*3600+45*60)<t1&(12*3600+15*60)>t1);
% TOF1 = TOF1(www);
% t1=t1(www);
% residual1 = residual1(www);
% TOF1 = TOF1(pos);
w = 0;
%һ��ֱ�������ز��ź�
for i = -0.1:0.001:0.1 %б�ʷ�Χ
    
    for j=-100:1:100 %�ؾ�
        y=i*(t1-min(t0))+j;
        w1 = [];
        w1 = find(residual1>y&residual1<(y+oo));
        ww = length(w1);
        if w<ww
            w = ww;
            p = i; %б��
            q = j; %�ؾ�
            w2 = w1; %��¼λ����Ϣ
        end
    end
end
t2 = t1(w2);
residual2 = residual1(w2);
TOF2 = TOF1(w2);
hold on
plot(t2-min(t0),residual2,'.b','MarkerSize',12);


y=p*(t1-min(t0))+q;
hold on
plot(t1-min(t0),y,'m')
hold on
plot(t1-min(t0),y+oo,'m')
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('O-C Residual/ns','fontsize',40);

% �����������ڵĲв����2.5���˲�n�Σ��õ���Ϊ�ǻز��źŵ�TOF3
n=10; %�˲�����
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
    t2=t2(w2);
    residual2=residual2(w2);
    TOF2=TOF2(w2);
end
t3=t2;
residual3=residual2;
TOF3=TOF2;
hold on 
plot(t3-min(t0),residual3,'.r','MarkerSize',12);

figure
subplot(1,2,1);
plot(t0-min(t0),residual0,'.','MarkerSize',12);
hold on
plot(t2-min(t0),residual2,'.g','MarkerSize',12);
hold on
plot(t1-min(t0),y,'m')
hold on
plot(t1-min(t0),y+oo,'m')
set(gca,'FontSize',40);
plot(t3-min(t0),residual3,'.r','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time [s]','fontsize',40);
ylabel('O-C Residual [ns]','fontsize',40);
%set(gcf,'outerposition',get(0,'screensize'));	% ����ͼ�δ���λ�ú���ߴ�Ϊ��Ļ��С
%saveas(gcf,'residual.jpg');
subplot(1,2,2);

BW = 3e-2;%Binwidth:30ps

num3=round((max(residual3)-min(residual3))/(BW));
N=length(residual3);%photon counts
x0 = min(residual3):(max(residual3)-min(residual3))/(num3-1):max(residual3);
y0 = hist(residual3,num3);
xx = x0(:); %�в�ns
yy = y0(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w4 = gfit.c1/sqrt(2); % round trip sigma :ns
xc4 = gfit.b1;% ϵͳ��ʱ����ֵ
%������ȫ��FWHM��
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc4)/gfit.c1).^2);
histogram(residual3,'BinWidth',BW,'Orientation','horizontal') 
set(gca,'FontSize',40);
xlabel('Photon counts','fontsize',40);
ylabel('O-C Residual [ns]','fontsize',40);

hold on
plot(yyy,xxx)

set(gcf,'outerposition',get(0,'screensize'));	% ����ͼ�δ���λ�ú���ߴ�Ϊ��Ļ��С
saveas(gcf,'Residual&historgram.jpg');
% �Իز��ź�����TOF3���ж���ʽ��ϣ���n��n+1����ϱ�׼��С��1mmʱ��n�׽�����ϣ��õ���׼��
%
degree=[];
for n=1:15    
    [p3,S3]=polyfit(t3,TOF3,n);
    [y_fit3,delta3]=polyval(p3,t3,S3);
    rms(n)=mean(delta3)*c/2;
    degree=[degree;n rms(n)*100];
    if n>1&abs(rms(n)-rms(n-1))<2.5e-3%�жϵ�n��n+1����ϱ�׼��Ĳ�ֵ�Ƿ�С��
        break
    end
end
clear rms

h=floor(max(t3)/3600)+8;
if h>=24
    h=h-24;
end
m=floor(mod(max(t3),3600)/60);
s=((mod(max(t3),3600)/60)-floor(mod(max(t3),3600)/60))*60;
h0=floor(min(t3)/3600)+8;
if h0>=24
    h0=h0-24;
end
m0=floor(mod(min(t3),3600)/60);
s0=((mod(min(t3),3600)/60)-floor(mod(min(t3),3600)/60))*60;



figure
plot(t3,TOF3,'.r','MarkerSize',12)
hold on
plot(t3,y_fit3,'-r','MarkerSize',20)
set(gca,'FontSize',40);
xlabel('Time/s','fontsize',40);
ylabel('RTT/s','fontsize',40);
legend('Observed RTT', 'Polynomial fitting')
% TOF-��Ϻ����ߵõ��в�residual4
residual4 = (TOF3 - y_fit3)*1e9;%ns
figure
plot(t3,residual4,'.r','MarkerSize',12);
set(gca,'FontSize',40);
xlabel('Time/seconds','fontsize',40);
ylabel('Residual/ns','fontsize',40);
% figure
% histogram(residual4,'BinWidth',0.02) %��״ͼbin��Сѡȡ20ps
% set(gca,'FontSize',40);
% xlabel('O-C residual/s','fontsize',40);
% ylabel('Photon counts','fontsize',40);

BW = 5e-3;%Binwidth:50ps
figure
histogram(residual4,'BinWidth',BW) 
set(gca,'FontSize',40);
xlabel('O-C residual(ns)','fontsize',40);
ylabel('Photon counts','fontsize',40);
num4=round((max(residual4)-min(residual4))/(BW));
N=length(residual4);%photon counts
x = min(residual4):(max(residual4)-min(residual4))/(num4-1):max(residual4);
y = hist(residual4,num4);
xx = x(:); %�в�ns
yy = y(:);
% Set up fittype and options.
gfit = fit(xx,yy,'gauss1'); % gfit(x) =  a1*exp(-((x-b1)/c1)^2),sigma=c1/sqrt(2)
a = gfit.a1;
w4 = gfit.c1/sqrt(2); % round trip sigma :ns
xc4 = gfit.b1;% ϵͳ��ʱ����ֵ
%������ȫ��FWHM��
xxx=(min(xx)-0.1):0.0001:(max(xx)+0.1);
yyy=a*exp(-((xxx-xc4)/gfit.c1).^2);
hold on
plot(xxx,yyy)
set(gcf,'outerposition',get(0,'screensize'));	% ����ͼ�δ���λ�ú���ߴ�Ϊ��Ļ��С
saveas(gcf,'save.jpg');
FWHM4 = xxx(max(find(abs(yyy-max(yyy)/2)<2)))-xxx(min(find(abs(yyy -max(yyy)/2)<2)));
noise = (length(residual0)-length(residual4))/(2000*(max(t0)-min(t0)))*(max(residual4)-min(residual4))*(max(t3)-min(t3));%����ˮƽ
n_echo = round(length(residual4)-noise)
y=[w4*1e3,w4*2.355*1e3,n_echo]; %˫�̱�׼�ps�� FWHM:ps�� �۳����������Ч�ز�����counts
w4*0.15
save y.txt -ascii y
format long g
w4*1e-9*c/2/sqrt(N)*1e3
w4/sqrt(N)*1e3/2
[h0 m0 s0 h m s max(t3)-min(t3) length(residual3) n rms(residual4)*c/2*1e-7 max(TOF3)*c/2*1e-3 n_echo/(max(t3)-min(t3))]
%% ���ܹ������