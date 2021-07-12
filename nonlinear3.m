function [c, ceq] = nonlinear3(Dw)
xi_min = 0;
xi_max = 2.5;
alpha_min = 2;
alpha_max = 10;
R_max =  30; %kW


% Other Parameters
lambda = 20;
M_cursive = 110; 
R_cursive = 1400;

% Compute E[x]
x_min = 10;
x_max = 100;
E_x = .5*(x_max^2 - x_min^2)/(x_max - x_min);
E_x2 = (x_max^3 - x_min^3)/(2*x_max - 3*x_min);


% Compute E[u]
% E_u = integral(@(u)  u.*pdfu(u, Dw(1), Dw(2), x_max, x_min, alpha_max, alpha_min, xi_max, xi_min), -1, 100)
% E_1_u = integral(@(u)  (1./u).*pdfu(u, Dw(1), Dw(2), x_max, x_min, alpha_max, alpha_min, xi_max, xi_min), -1, 10)
% 
numofpts = 20000;
xj = unifrnd(x_min, x_max, numofpts,1);
xij = unifrnd(xi_min, xi_max, numofpts,1);
alphaj = unifrnd(alpha_min, alpha_max, numofpts,1);

uj = max(xij, -1./(2*Dw(1)) * (alphaj./xj) + Dw(2) );

E_u = mean(uj);
E_1_u = mean(1./max(xij, -1./(2*Dw(1)) * (alphaj./xj) + Dw(2)));

% Compute E[r] & E[\theta]
E_r = E_x*E_1_u;
E_r2 = mean(xj.^2./uj.^2);
% Compute delta(M)
num = M_cursive - lambda*E_u;
denom = 2*(lambda*E_u + (num)/3);
domain_logic1 = (M_cursive/lambda > E_u);
domain_logic2 = (M_cursive/lambda <= E_u);
deltaM = exp(-(num)^2/(denom))*domain_logic1 + 1*domain_logic2;
conf_M = 1 - deltaM;
c(1) = .30 - 1 + deltaM;

% Compute gamma(R)
m = linspace(ceil(R_cursive/R_max),floor(R_cursive/E_r),...
                -ceil(R_cursive/R_max) +floor(R_cursive/E_r) + 1);
R_num = -(R_cursive - m*E_r).^2;
R_denom = 2*(m*E_r2 + R_max*(R_cursive - m*E_r)/3);
sum_terms = exp(R_num./R_denom).*poisspdf(m,lambda);

R_domain_logic1 = (R_cursive/lambda > E_u*E_r);
R_domain_logic2 = (R_cursive/lambda <= E_u*E_r);
gammaR = min([1, sum(sum_terms)])*R_domain_logic1 + ...
                1*R_domain_logic2;
            
conf_R = 1 - gammaR;
c(2) = .75 - 1 + gammaR;
c(3) = Dw(1) - Dw(2)*R_max;
% T = table(conf_M, conf_R, E_u,  'VariableNames', {'1 - del(M)', '1-gam(R)', 'E[u]'})
ceq = [];
end