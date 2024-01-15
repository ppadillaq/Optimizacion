%% REVOLUTION SURFACE AREA MINIMIZATION
clear
clc

% Catenaria
xCat = 0:0.01:1;
a = 424169/500000;
yCat = a*cosh((xCat-0.5)/a);

f = @integral;
% n = núm. de puntos
n=[5 20 25 50 100 1000];
% maxStep = paso máximo para búsqueda dicotómica
maxStep = [0.5 10 10 10 10 10];
%times = vector de tiempos de ejecución
times=zeros(1,length(n));
%errores = vector de errores cometidos
errores=zeros(1,length(n));
%iters = vector de números de iteraciones empleadas 
iters=zeros(1,length(n));
%aprox = vector de aproximaciones a la integral
aprox=zeros(1,length(n));
for i=1:length(n)
U0 = ones(1,n(i));

%Optimizamos y tomamos tiempo
tic
[U_sD,S_sD,fX_sD,i_sD] = steepestDescent(f, U0,'maxStep', maxStep(i), 'oneDimSearch', "dichoSearch");
times(i)=toc;

%Iteraciones a graficar
iterToGraf=[1 floor(i_sD*0.25) floor(i_sD*0.5) floor(i_sD*0.75) i_sD-1];

%Grafica Evolucion
figure('Name',"n="+num2str(n(i)));
% Grafica Catenaria
plot(xCat,yCat)
hold on
% Graficas de las iteraciones
x = linspace(0,1,numel(U0));
for j=1:5
s=scatter(x, U_sD(iterToGraf(j),:),15);
s.MarkerEdgeColor=[iterToGraf(j)/i_sD 0 0];
end
ylim([0.8, 1])
legend(['Catenaria' "iter="+string(iterToGraf)])
title("Evolución n="+num2str(n(i)))
xlabel("X")
ylabel("Y")
hold off
errores(i)=((0.9536241052-fX_sD(i_sD-1))/0.9536241052)*100;
aprox(i)=fX_sD(i_sD-1);
iters(i)=i_sD;
end

%Grafica Tiempo 
figure('Name','Costo Temporal')
bar(categorical(n),times)
text(categorical(n),times,num2str(times'),'vert','bottom','horiz','center'); 
xlabel("n")
ylabel("Time (s)")
title('Costo Temporal')



