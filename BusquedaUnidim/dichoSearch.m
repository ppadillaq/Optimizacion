function [x_opt, f_opt, x, intervalo, L, n, elapsedTime] = dichoSearch(func,xA, xB, delta, tol)
    %DICHOSEARCH Dichotomic Search Method
    %   Detailed explanation goes here
    %   Implementaci'on basada en el algoritmo expuesto en:
    %   ------------------------------------------------------------------
    %   Rao, S. S. (2019). Engineering optimization: theory and practice. 
    %   John Wiley & Sons.
    %   ------------------------------------------------------------------
    %
    %   Inputs:
    %     func: unidimensional function
    %     xA
    %     xB
    %     delta
    %     tol: tolerancia del error de aproximaci'on (valor porcentual)
    %
    %   Outputs:
    %     x_opt: optimum point
    %     f_opt: valor 'optimo de la funci'on (f(x_opt))
    %     x: array de pares de puntos evaluados (matriz n/2 x 2)
    %     intervalo: array de extremos de intervalo en cada iteraci'on
    %     L: array de longitudes de intervalo en cada iteraci'on
    %     n: n'umero de experimentos (evaluaci'on de la funcion)
    %     elapsedTime: tiempo de ejecuci'on
    
    timerVal = tic;
    
    % longitud del intervalo inicial
    L0 = xB - xA;

    % estimacion del numero de experimentos
    minimo = (L0 - delta)/(2*tol*L0 - delta);
    n = 2;
    valorTest = 2^(n/2);
    while valorTest < minimo
        n = n + 2;
        valorTest = 2^(n/2);
    end

    x = zeros(n/2,2);
    intervalo = zeros(n/2 + 1, 2);
    L = zeros(n/2 + 1, 1);
    % intervalo inicial
    intervalo(1,:) = [xA xB];
    L(1) = L0;
    for i = 1:(n/2)
        x(i,1) = intervalo(i,1) + L(i)/2 - delta/2;
        x(i,2) = intervalo(i,1) + L(i)/2 + delta/2;
        f1 = feval(func, x(i,1));
        f2 = feval(func, x(i,2));
        if f1 < f2
            intervalo(i+1,1) = intervalo(i,1);
            intervalo(i+1,2) = x(i,2);
        elseif f1 > f2
            intervalo(i+1,1) = x(i,1);
            intervalo(i+1,2) = intervalo(i,2);
        else
            intervalo(i+1,:) = [x(i,1) x(i,2)];
        end
        L(i+1) = intervalo(i+1,2) - intervalo(i+1,1);
    end

    % the optimus point is the mean point of final interval
    x_opt = mean(intervalo(end,:));
    f_opt = feval(func, x_opt);

    elapsedTime = toc(timerVal);
end