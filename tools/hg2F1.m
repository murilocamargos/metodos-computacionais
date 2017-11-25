function [ F ] = hg2F1( a, b, c, z, N )
%HG2F1 Função hipergeométrica 2F1
% o parâmetro N é opcional e é utilizado na integral de Chebyshev
if ~exist('N', 'var')
    N = 22;
end
F = gamma(c)/(gamma(b)*gamma(c-b)) * chebyshev(@(x) (x.^(b-1) .* (1-x).^(-b+c-1)) ./ ((1-z*x).^a), eps, 1-eps, N);
end

