%{
...
How to do a lagrange Interpolation?
It is used to find the polynomial approximation of the data points .
If we have 'n' data points then we can approximate by 'n-1' th degree
polynomial

Created on 07/02/2020 18:07 By Karthi 
...
%}
clearvars;clc

% x = [-2 0 2];% 
% y = [4 2 8];
% x = [0 0.5 1 1.5 2];
% y = [0 0.19 0.26 0.29 0.31];
% x = [1 3 4 6];
% y = [7 53 157 857];
% 
% 
% 
% [xVal,Yval] = LagrangeInterpol(x,y);
% plot(xVal,Yval)
% grid on 
% xlabel('\it{x values}')
% ylabel('\it{Interpolated y values}')

X = [0 pi/6 pi/4 pi/3 pi/2];
Y = sin(X);
x = linspace(0,pi/2,10000);
M = 3;
y = lagrange1( x, M,X, Y);
y1 = sin(x);
figure
errorbar(x,y,R,'.g')
hold on
plot(X, Y, 'or', x, y, '.k', x, y1, '-b');
legend('误差','样本点','拉格朗日插值估算','sin(x)');
          
x = sp3t_s_w7(1:2:end)/3600;
y = sp3x_w7(1:2:end)/1e3;
xi = sp3t_s_w7(2:2:end-1)/3600;
yi = (sp3x_w7(2:2:end-1))'/1e3;

Y = lagint( x, y,xi,10);
% 
% Y = (interp1(x, y,xi, 'spline'))';
figure
plot(xi,(Y-yi)*1e9,'.b','MarkerSize',12)
set(gca,'FontSize',40);
xlim([0,24])
xlabel('Time [hours]','fontsize',40);
ylabel('X coordiantes [millimeters]','fontsize',40);
figure
errorbar(xi,yi,Y-yi,'.g','MarkerSize',12)
hold on
plot(x, y, 'or', xi, Y, '.k', xi, yi, '-b','MarkerSize',8);
legend('误差','样本点','拉格朗日插值估算','sp3 data','fontsize',20);
set(gca,'FontSize',40);
xlabel('Time [hours]','fontsize',40);
ylabel('X coordiantes [Megameters]','fontsize',40);
xlim([0,24])


set(gca,'FontSize',40);
xlabel('Time [Hours]','fontsize',40);
ylabel('X coordiantes [Megameters]','fontsize',40);