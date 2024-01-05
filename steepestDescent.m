function [X,S,fX,i] = steepestDescent(f, X0, varargin)
    %STEEPESTDESCENT Summary of this function goes here
    %   Detailed explanation goes here
    %   Pedro Padilla Quesada
    
    maxIter = 1e4; % maximum no. of iterations
    converged = [0 0 0]; % convergence flags (3 criteria)
    i = 1; % iteration counter

    % One-dimensional minimization default method: Dichotomous Search
    delta = 1e-6; % eval margin
    tol = 1e-7; % error tolerance
    oneDimSearch = @(f) dichoSearch(f,0, 10, delta, tol);

    % is optional arguments are introduced:
    if ~isempty(varargin)
        props = varargin(1:2:numel(varargin));
        values = varargin(2:2:numel(varargin));
        for j = 1:numel(props)
            switch props{j}
                case 'maxIter'
                    maxIter = values{j};
                case 'oneDimSearch'
                    if values{j} == "dichoSearch"
                        % default method
                    elseif values{j} == "quasiNewton1dim"
                        step = 0.01;
                        oneDimSearch = @(f) quasiNewton1dim(f, 0, step);
                    end
                otherwise
                    disp('Unknown function argument.')
            end
        end
    end

    % Iterations matrix - Preallocation
    X = NaN(maxIter, numel(X0));
    S = NaN(maxIter, numel(X0));
    fX = NaN(maxIter, 1);
    X(i,:) = X0;
    [fX(i), delta_f] = f(X0);
    
    while ~any(converged) && i < maxIter
        fprintf('Iteration %d...\n',i)
        %delta_f = [gradient{1}(X(i,:)) gradient{2}(X(i,:))];
        S(i,:) = -1*delta_f;
        wrapper = @(lambda) wrapper1dim(f, lambda, X(i,:), S(i,:));
        [lambda_opt, fX(i+1)] = oneDimSearch(wrapper);
        X(i+1,:) = X(i,:) + lambda_opt * S(i,:);
        [~, delta_f] = f(X(i+1,:));
        converged = convergenceCriteria(X(i+1,:), X(i,:), fX(i+1),fX(i),delta_f, 1e-8,1e-8,1e-8);
    
        i = i + 1;
    end
    
    if i > maxIter
        disp('The solution found has not reached convergence criteria in the given maximum number of iterations.')
    else
        idx = find(converged);
        for k = 1:numel(idx)
            fprintf('Convergence criteria %d is met, at iteration %d.\n', idx(k), i-1)
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