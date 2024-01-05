function [I, grad] = integral(u)
    %INTEGRAL Integral of trapecia for revolution surface area calculation.
    %   It serves as a minimization cost function to calculate the minimal 
    %   surface.
    %
    %   Inputs:
    %     u      : discretized generatrix curve.
    %
    %   Outputs:
    %     I      : revolution surface area (scaled by 1/(2*pi))
    %     grad   : function gradient
    %

    % number of segments
    n = numel(u)-1;

    T = zeros(n,1);
    grad = zeros(1,n+1);
    for j = 1:n
        T(j) = ((u(j+1)+u(j))*sqrt(1 + n^2*(u(j+1)-u(j))^2))/(2*n);
    end
    I = sum(T);

    % Gradient calculation
    % a = sqrt(1+n^2*(u(2) - u(1))^2);
    % grad(1) = (a - n^2*(u(2) - u(1))*(u(2) + u(1))/a)/(2*n);
    for j = 2:n
        a = sqrt(1 + n^2*(u(j) - u(j-1))^2);
        b = sqrt(1 + n^2*(u(j+1) - u(j))^2);
        grad(j) = (n^2*(u(j) - u(j-1))*(u(j) + u(j-1))/a + a + b - n^2*(u(j+1) - u(j))*(u(j+1) + u(j))/b)/(2*n);
    end
    % a = sqrt(1+n^2*(u(n) - u(n-1))^2);
    % grad(n+1) = (n^2*(u(n) - u(n-1))*(u(n-1) + u(n))/a + a)/(2*n);
end