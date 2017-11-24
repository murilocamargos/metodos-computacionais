function [ I, ar, br ] = chebyshev( F, a, b, N )
%CHEBYSHEV Calcula a integral numérica de uma função F no intervalo [a,b]
% utilizando o método de Chebyshev.
%
% A função F precisa conter apenas operações bitwise.

%% Matriz de pesos da soma finita com primeiro e último pesos iguais a 0.5
w = [.5 ones(1, N-1) .5];

%% Cria iteradores de 0 a N
[s,r] = meshgrid(0:N, 0:N);

%% Calcula os valores dos coeficientes a_r
t = cos(s * pi / N) * (b - a)/2 + (b + a)/2;
ar = F(t) .* cos(pi * s .* r / N);
ar = (2/N) * ar*w';

%% Calcula os valores dos coeficientes b_r
r = 2:N;
br = (ar(r-1) - ar(r+1))./(2*(r-1)');

%% Calcula o valor da integral
I = 2*sum(br(1:2:end)) * ((b-a)/2);

end