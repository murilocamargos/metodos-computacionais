function [ x ] = rk4( f, h, x0, intv )
%RK4 Summary of this function goes here
%   Detailed explanation goes here

t = intv(1):h:intv(2);               % Intervalo da simula��o
x = zeros(length(t), length(x0));    % Pr�-aloca��o da matriz de resultados
x(1,:) = x0;                         % Condi��es iniciais

%% Encontra solu��es do sistema de eq. diferenciais utilizando o m�todo de Runge-Kutta de quarta ordem.
for i = 1:length(t)-1
    k1 = h * f(t(i), x(i, :));
    k2 = h * f(t(i) + h/2, x(i, :) + k1/2);
    k3 = h * f(t(i) + h/2, x(i, :) + k2/2);
    k4 = h * f(t(i) + h, x(i, :) + k3);
    x(i + 1, :) = x(i, :) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

end

