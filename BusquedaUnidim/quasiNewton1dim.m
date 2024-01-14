function [lambda_opt, f_opt, lambda, f, i, elapsedTime] = quasiNewton1dim(func, lambda0, step, varargin)
    %QUASINEWTON1DIM Quasi-Newton method
    %   Direct rood methods.
    %   Adaptaci'on del M'etodo de Newton para el caso en el que no es
    %   posible o es dif'icil derivar la funci'on. Emplea diferencias
    %   centrales para aproximar la primera y la segunda derivada.
    %
    %   Implementation based in the algorithm as exposed in:
    %   ------------------------------------------------------------------
    %   Rao, S. S. (2019). Engineering optimization: theory and practice. 
    %   John Wiley & Sons.
    %   ------------------------------------------------------------------
    %   Pedro Padilla Quesada. 2023
    %
    %   Outputs:
    %     lambda
    %     i
    %     elapsedTime
    %

    % Default parameters
    maxIter = 100;  % maximum number of iterations
    eps     = 1e-6; % convergence tolerance

    % is optional arguments are introduced:
    if ~isempty(varargin)
        props = varargin(1:2:numel(varargin));
        values = varargin(2:2:numel(varargin));
        for j = 1:numel(props)
            switch props{j}
                case 'maxIter'
                    maxIter = values{j};
                case 'eps'
                    eps = values{j};
                otherwise
                    disp('Unknown function argument.')
            end
        end
    end

    % Start clock
    timerVal = tic;

    % Pre-allocation
    lambda = zeros(maxIter,1);
    lambda(1) = lambda0;

    f     = feval(func,lambda(1));
    f_izq = feval(func,lambda(1) - step);
    f_dch = feval(func,lambda(1) + step);
    df = (f_dch - f_izq)/(2*step);
    ddf = (f_dch - 2*f + f_izq)/step^2;
    i = 2;
    while abs(df) > eps 
        lambda(i) = lambda(i-1) - df/ddf;
        f     = feval(func,lambda(i));
        f_izq = feval(func,lambda(i) - step);
        f_dch = feval(func,lambda(i) + step);
        df = (f_dch - f_izq)/(2*step);
        ddf = (f_dch - 2*f + f_izq)/step^2;
        i = i + 1;
        if i > maxIter
            disp('Maximum number of iterations has been reached. Quasi-Newton method does not converge.')
            break;
        end
    end
    lambda_opt = lambda(i-1);
    f_opt = feval(func, lambda_opt);
    elapsedTime = toc(timerVal);
end