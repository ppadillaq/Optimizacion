function [X_data,S_data,fX_data,i] = steepestDescent(f, X0, varargin)
    %STEEPESTDESCENT Método del descenso más rápido (Cauchy).
    %   Implementación basada en el algoritmo expuesto en:
    %   ------------------------------------------------------------------
    %   Rao, S. S. (2019). Engineering optimization: theory and practice. 
    %   John Wiley & Sons.
    %   ------------------------------------------------------------------
    %   El uso del vector gradiente negativo como dirección de minimización
    %   fue primeramente empleado por Cauchy en 1847. En este método se
    %   comienza en un punto semilla inicial X1 e iterativamente se mueve a
    %   lo largo de direcciones de descenso más pronunciadas hasta que el
    %   punto óptimo es encontrado. Este método se puede resumir en los
    %   pasos siguientes:
    %
    %   1. Comenzar con un punto arbitrario initial X1. Establecer el
    %   número de iteración como i = 1.
    %   2. Encontrar la dirección de búsqueda S_i como
    %            S_i = -Delta f(Xi)
    %   3. Determinar la longitud de paso óptima lambda_i^* en la dirección
    %      S_i y establecer:
    %            Xi+1 = Xi  +  lambda_i^* · S_i
    %   4. Evaluar el nuevo punto Xi+1 para optimalidad (según criterios de 
    %      convergencia). Si Xi+1 es óptimo, detener el proceso. De otro 
    %      modo, ir al paso 5.
    %   5. Establecer i = i + 1 e ir al paso 2.
    %
    
    maxIter = 1e6; % maximum no. of iterations
    converged = [0 0 0]; % convergence flags (3 criteria)
    i = 1; % iteration counter

    % Tolerancias de convergencia
    eps1 = 1e-8; % Criterio 1
    eps2 = 1e-8; % Criterio 2
    eps3 = 1e-8; % Criterio 3

    % si argumentos opcionales son introducidos:
    if ~isempty(varargin)
        props = varargin(1:2:numel(varargin));
        values = varargin(2:2:numel(varargin));
        for j = 1:numel(props)
            switch props{j}
                case 'maxIter'
                    maxIter = values{j};
                case 'maxStep'
                    maxStep = values{j}; % paso máx. para búsqueda dico.
                case 'oneDimSearch'
                    % Método de minimización unidimensional
                    if values{j} == "dichoSearch"
                        % Búsqueda Dicotómica
                        delta = 1e-6; % margen de evaluación
                        tol = 1e-7; % tolerancia del error
                        oneDimSearch = @(f) dichoSearch(f,0, maxStep, delta, tol);
                    elseif values{j} == "quasiNewton1dim"
                        % Método Cuasi-Newton 1D
                        step = 0.01;
                        oneDimSearch = @(f) quasiNewton1dim(f, 1, step);
                    end
                otherwise
                    disp('Unknown function argument.')
            end
        end
    end

    % Matriz de iteraciones - Preallocation
    X_data = NaN(maxIter, numel(X0));
    S_data = NaN(maxIter, numel(X0));
    fX_data = NaN(maxIter, 1);
    
    [fX0, delta_f] = f(X0);
    
    while ~any(converged) && i < maxIter
        fprintf('Iteración %d...\n',i)
        %delta_f = [gradient{1}(X(i,:)) gradient{2}(X(i,:))];
        S = -1*delta_f;
        wrapper = @(lambda) wrapper1dim(f, lambda, X0, S);
        [lambda_opt, fX] = oneDimSearch(wrapper);
        X = X0 + lambda_opt * S;
        [~, delta_f] = f(X);
        converged = convergenceCriteria(X, X0, fX, fX0, delta_f, eps1, eps2, eps3);

        % graba resultados
        X_data(i,:) = X0;
        fX_data(i) = fX0;
        S_data(i,:) = S;
    
        i = i + 1;
        X0 = X;
        fX0 = fX;
    end
    
    if i > maxIter
        disp('La solución encontrada no ha alcanzado el criterio de convergencia en el número máximo de iteraciones seleccionado.')
    else
        idx = find(converged);
        for k = 1:numel(idx)
            fprintf('Criterio de convergencia %d alcanzado, en la iteración %d.\n', idx(k), i-1)
        end
    end

end

function converged = convergenceCriteria(X, X0, fX, fX0, delta_f, eps1, eps2, eps3)
    % Criterion 1
    c1 = abs((fX - fX0)/fX0) <= eps1;
    % Criterion 2
    c2 = all(abs(delta_f) <= eps2);
    % Criterion 3
    c3 = norm(X - X0) <= eps3;
    converged = [c1, c2, c3];
end