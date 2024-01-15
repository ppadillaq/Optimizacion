function out = wrapper1dim(funcionOriginal, lambda, X0, S)
    %WRAPPER1DIM Función envoltorio para emplear métodos de búsqueda
    %   unidimensional en problemas de optimización multidimensional.
    %   Realiza un cambio de variable en la función objetivo. Sustituye el
    %   vector de variables de decisión X por una variable escalar lambda,
    %   que representa el paso a optimizar en la dirección de descenso S.
    %
    %   Inputs:
    %     funcionOriginal: función objetivo original (multidimensional)
    %     X0: punto inicial
    %     S: dirección de descenso
    %     lambda: paso en la dirección S (variable de decisión)
    
    X = X0 + lambda * S;
    out = funcionOriginal(X);

end