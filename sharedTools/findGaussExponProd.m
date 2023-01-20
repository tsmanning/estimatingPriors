function [fxn,fmax] = findGaussExponProd(x,a,b,mu,sigma)
    
    % Product of exponential distribution and Gaussian and max

    fxn = (1/(sigma*sqrt(2*pi)))*exp(a*x + b - (1/(2*sigma^2))*(x-mu).^2);    
    
    fmax = a*sigma + mu;
    
end
