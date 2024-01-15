%% REVOLUTION SURFACE AREA MINIMIZATION
clear
clc

% Draw catenary curve
x = 0:0.01:1;
a = 424169/500000;
y = a*cosh((x-0.5)/a);
plot(x,y)
hold on

% Plot optimum discretized curve that minimizes the revolution surface area
f = @integral;
U0 = ones(1,5);
[U_sD,S_sD,fX_sD,i_sD] = steepestDescent(f, U0, 'maxStep', 0.1, 'oneDimSearch', "dichoSearch");
x = linspace(0,1,numel(U0));
scatter(x, U_sD(i_sD-1,:))
ylim([0, 1])

U0 = ones(1,100);
[U_sD,S_sD,fX_sD,i_sD] = steepestDescent(f, U0, 'maxStep', 10, 'oneDimSearch', "dichoSearch");
x = linspace(0,1,numel(U0));
figure
scatter(x, U_sD(i_sD-1,:))
ylim([0, 1])