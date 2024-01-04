%% Parametros
n = 1000;  %numero de intervalos a calcular
eps = 1e-8;  %valor de la normal del salto
tol = 1e-4;   %tolerancia para budi
tic();

%% Inicializacion de Variables
X = linspace(0,1,n+1); %Divisiones del eje x
U = ones(1,n-1);       %Se empieza a iterar con todos los uj igual a 1. Es decir con la recta que une (0,1) y (1,1)

%% Incrementar n hasta llegar al objetivo
% Se utiliza n bajo para aproximar el valor rápidamente y tras converger se
% incrementa n interpolando valores para tener nueva solucion inicial más
% cercana a la que debe converger
n_obj = n;
n = 1; %n inicial
iter_tot = 1;
while n<=n_obj
    fprintf('n=%d\n',n);
    %interpolar U con nuevo n desde la U previa
    X_new = linspace(0,1,n+1);
    U_new = interp1(X,[1,U,1],X_new); %interpolar con ambos extremos fijados a 1
    U_new = U_new(2:end-1);           %eliminar los extremos fijados a 1
    X=X_new;
    U=U_new;
    dU = deriv(U);         %Cálculo del gradiente en el punto inicial
    FU = fobjetivo(U);     %Cálculo del valor de la funcion objetivo en el punto inicial
    Uold = U*Inf;          %Valor "viejo" para detener el bucle
    %% Método del gradiente
    iter = 1;
    while norm(U-Uold)>eps
        Uold = U; %Se guarda el vector de U actual para hacer la comparación de la norma
        tmin = budi(@fobjetivo, U, dU, 0.3, tol/10, tol);

        %% ejecutar movimiento en direccion gradiente negativo con el paso calculado
        U = U - tmin*dU;   %Se actualiza el vector U
        FU = fobjetivo(U); %Actualizar el valor de la funcion objetivo en el nuevo punto
        dU = deriv(U);     %Actualizar el valor del gradiente
        if mod(iter,1000)==0
            fprintf(1,'n=%d | iter=%d | tmin=%e | Fobj(U)=%e | norm(U-Uold)=%e\n',n,iter,tmin,FU,norm(U-Uold));
        end
        iter = iter + 1;
    end
    iter_tot = iter_tot + iter;
    n_factor = 2; %factor para calcular nueva n
    %salir del bucle si n ya es n objetivo
    if n==n_obj, break; end
    %calcular nuevo valor de n
    if n*n_factor>n_obj
        n=n_obj;
    else
        n=floor(n*n_factor);
    end
end
fprintf(1,'El número de iteraciones totales ha sido %d en %e segundos\n',iter_tot,toc());
% Añadir valores inicial y final fijados a 1 al resultado
U = [1,U,1];

%% Trazar resultado y función objetivo analítica para comparar
scatter(X,U); %Se dibujan los puntos del vector U
% Parámetros de la catenaria
a = 424169/500000; % Valor de a para que la catenaria pase por (0,1) y (1,1)
b = 0.5; %eje de simetría
x = X;
y = a * cosh((x-b) / a);  % Calcular y en función de x
% Trazar la gráfica
hold on
plot(x, y);
title('Gráfica de Catenaria');
xlabel('x');
ylabel('u');
grid on;
ylim([0, 1]);

%% Funcion objetivo
function FU = fobjetivo(U)
    n = numel(U)+1;
    Uk  = [1,U];
    Ukk = [U,1];
    FU = sum((Uk+Ukk).*sqrt(1+(n.*(Ukk-Uk)).^2)./(2*n));
end

%% Función de cálculo del gradiente simbólico
% La derivada de los extremos es innecesaria, no se utiliza para optimizar
% y el valor de la funcion en ambos puntos está fijado a 1.
function dFU = deriv(U)
    n = numel(U)+1;
    % Se añade puntos extremos con valor 1
    Ukm1 = [1,U(1:end-1)]; %U_(k-1)
    Uk   = U;              %U_k
    Ukp1 = [U(2:end),1];   %U_(k+1)
    dFU = 1/(2*n) * (sqrt(1+n^2*(Uk-Ukm1).^2)+(n^2*(Uk-Ukm1).*(Uk+Ukm1)./sqrt(1+n^2*(Uk-Ukm1).^2))) + 1/(2*n)*(sqrt(1+n^2*(Ukp1-Uk).^2)-(n^2*(Ukp1-Uk).*(Ukp1+Uk)./sqrt(1+n^2*(Ukp1-Uk).^2)));
end

%% función budi con descenso de gradiente integrado
function xmin = budi(f, U, dU, tini, eps, tol)
    ak = 0;
    bk = tini;
    while abs(bk-ak) >= tol
        x1 = (ak + bk)/2;
        xa = x1 - eps;
        xb = x1 + eps;
        fxa = f(U-dU*xa); %optimizar paso para descenso de gradiente
        fxb = f(U-dU*xb);
        if fxa < fxb
            bk = xb;
        elseif fxa >= fxb
            ak = xa;
        end
        abs(bk - ak);
    end
    xmin = (ak + bk)/2;
end
