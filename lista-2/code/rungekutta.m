clear
clc
format long

addpath('../../tools');

%% Parâmetros iniciais da equação da bolha
Pl = 0.195091;
Pv = 0.166391;
el = 4.34802;
ev = 1.79634;

%% Parâmetros da equação posta na forma R'(t)^2 = -b/a + kR(t)^(-2a)
a = 1.2386;
b = 0.0058149;
k = b/a;
t0 = 0;

f = @(t,X) [X(2); (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]

h = 0.1;                            % Tamanho do passo
d = 2;                              % Quantidade de equações no sistema
t = 0:h:12.5;                       % Intervalo da simulação
x = zeros(d, length(t));            % Pré-alocação da matriz de resultados
x(:,1) = [1; 0.0001];               % Condições iniciais

% Caso a velocidade seja menor que 1e-3, considero que é zero e posso
% calcular o tempo de colapso.
if x(2,1) < 1e-3
    tc = sqrt(a)/sqrt(b) * beta((a+1)/(2*a), 0.5)/(2*a) * x(1,1);
    t = 0:h:tc;
end

%% Encontra soluções do sistema de eq. diferenciais utilizando o método de Runge-Kutta de quarta ordem.
for i = 1:length(t)-1
    k1 = h * f(t(i), x(:, i));
    k2 = h * f(t(i) + h/2, x(:, i) + k1/2);
    k3 = h * f(t(i) + h/2, x(:, i) + k2/2);
    k4 = h * f(t(i) + h, x(:, i) + k3);
    x(:, i+1) = x(:, i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%% Plota as variáveis analisadas
hold on
plot(t, x(1,:));
plot(t, x(2,:));
xlabel('Tempo (s)');
ylabel('Variável medida');
grid on;
legend('Raio', 'Velocidade');
hold off

h = 0.01;
R = @(x) sqrt(a*k./(b + a^3*k^3*(x - t0).^2));
t = @(tau) chebyshev(@(x) (R(x)).^(2*a+1), 0, tau, 22);
tau = 0;
ts = [];
rs = [];
while 1
    rs(end+1) = R(tau);
    ts(end+1) = t(tau);
    tau = tau + h;
    if ts(end) > tc
        break
    end
end


% R = @(tau) sqrt(a*k./(b + a^3*k^2*(tau - t0).^2));
% F21 = @(a,b,c,z) gamma(c)/(gamma(b)*gamma(c-b)) * double(int(x^(b-1)*(1-z*x)^(-a)*(1-x)^(-b+c-1), x, 0, 1));
% t = @(tau) F21(0.5, 1+1/(2*a), 3/2, 1-R(tau)^(-2*a)) * sqrt(R(tau)^(-2*a) - 1)/(sqrt(a)*sqrt(b));
% ts = [];
% rs = [];
% for i = 0:1000
%     rs(end+1) = R(tau);
%     ts(end+1) = t(tau);
% end