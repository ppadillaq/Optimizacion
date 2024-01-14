function out = wrapper1dim(originalFunction, lambda, X0, S)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
X = X0 + lambda * S;
out = originalFunction(X);
end