function [c, ceq] = nonlinear(R)
% PMF of Alpha
pmf_alpha = [5, 10, 15];
pmf = [1/3; 1/3; 1/3]';

% Compute E[x]
x_min = 5;
x_max = 60;
E_x = .5*(x_max^2 - x_min^2)/(x_max - x_min);

% Other Parameters
lambda = 30;
M = 60; %59
R_cursive = 1000; %990

% % Compute probabilities - Old
% p1 = sum(pmf.*(pmf_alpha > max(.05*R(1)*R(2)/(R(2) - R(1)), .1*R(1)*R(3)/(R(3) - R(1)) )));
% p2 = sum(pmf.*(pmf_alpha > .05*R(2)*R(3)/(R(3) - R(2))).*(pmf_alpha < -.05*R(1)*R(2)/(R(1) - R(2))) );
% p3 = sum(pmf.*(pmf_alpha < min(-.05*R(2)*R(3)/(R(2) - R(3)),...
%           -.1*R(1)*R(3)/(R(1) - R(3)) )))

% Compute probabilities    
p1 = sum(pmf.*(pmf_alpha < min(.05*R(1)*R(2)/(R(2) - R(1)), .1*R(1)*R(3)/(R(3) - R(1)) )));
p2 = sum(pmf.*(pmf_alpha < .05*R(2)*R(3)/(R(3) - R(2))).*(pmf_alpha > -.05*R(1)*R(2)/(R(1) - R(2))) );
p3 = sum(pmf.*(pmf_alpha > max(-.05*R(2)*R(3)/(R(2) - R(3)),...
          -.1*R(1)*R(3)/(R(1) - R(3)) )));
      
% Compute E[r] & E[\theta]
E_r = (R(1)*p1 + R(2)*p2  + R(3)*p3);
E_r2 = (R(1)^2*p1 + R(2)^2*p2  + R(3)^2*p3);

E_theta = E_x*(1/R(1)*p1 + 1/R(2)*p2  + 1/R(3)*p3);
R_max = max(R);

% Compute delta(M)
num = M - lambda*E_x*(1/R(1)*p1 + 1/R(2)*p2  + 1/R(3)*p3);
denom = 2*(lambda*E_x*(1/R(1)*p1 + 1/R(2)*p2 + 1/R(3)*p3) + (num)/3);
domain_logic1 = (M/lambda > E_theta);
domain_logic2 = (M/lambda <= E_theta);
deltaM = exp(-(num)^2/(denom))*domain_logic1 + 1*domain_logic2;
conf_M = 1 - deltaM;
c(1) = .85 - 1 + deltaM;

% Compute gamma(R)
m = linspace(ceil(R_cursive/R_max),floor(R_cursive/E_r),...
                -ceil(R_cursive/R_max) +floor(R_cursive/E_r) + 1);
R_num = -(R_cursive - m*E_r).^2;
R_denom = 2*(m*E_r2 + R_max*(R_cursive - m*E_r)/3);
sum_terms = exp(R_num./R_denom).*poisspdf(m,lambda);

R_domain_logic1 = (R_cursive/lambda > E_theta*E_r);
R_domain_logic2 = (R_cursive/lambda <= E_theta*E_r);
gammaR = min([1, sum(sum_terms)])*R_domain_logic1 + ...
                1*R_domain_logic2;
            
conf_R = 1 - gammaR;
c(2) = .85 - 1 + gammaR;
T = table(conf_M, conf_R, E_theta, p1, p2, p3, 'VariableNames', {'1 - del(M)', '1-gam(R)', 'E[theta]', 'P1', 'P2', 'P3'})
ceq = [];
end