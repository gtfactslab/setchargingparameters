%% Choose Charging Rate
clc; close all; clear all;
x_min = 5;
x_max = 50;

R_min = 10
R_max = 60;
% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 1
% P ( a_ < a_j < a^)
fun = @(R) 1*(1/R(1) + 1/R(2) + 1/R(3))
% fun = @(R) 1*(R(1)^2 + R(2)^2 + R(3)^2)

% R1 - R2 < 0
A = [1 -1 0; 
     0 1 -1;
     1 0 -1;
     1 0 0; 
     0 1 0; 
     0 0 1;
     -1 0 0; 
     0 -1 0; 
     0 0 -1;];
  
b = [ 0;
      0;
      0;
      R_max - 1e-5; 
      R_max - 1e-5;        
      R_max - 1e-5; 
      0; 
      0;
      0];

% lb = [0;0];
% ub = [R_max; R_max];
nonlcon = @nonlinear
[x,fval,exitflag,~] = fmincon(fun,[27;30;34], A, b, [], [], [],[], nonlcon)


%% Full Simulation - Choose Charging Rate
% This code section runs the approximate chance-constrained optimization 
% Program 1 from the CCTA 2021 paper titled "Pricing Parameter Design for 
% Electric Vehicle Charging"

clc; close all; clear all;
x_min = 5;
x_max = 50;

R_min = 35;
R_max = 40;

% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 1
% P ( a_ < a_j < a^)

%Declare Objective Function
fun = @(R) 1*(1/R(1) + 1/R(2) + 1/R(3))

% Set constraints A*x < b
A = [1 -1 0; 
     0 1 -1;
     1 0 -1;
     1 0 0; 
     0 1 0; 
     0 0 1;
     -1 0 0; 
     0 -1 0; 
     0 0 -1;];
  
b = [ -4; %-4
      -5; % -5
      0;
      R_max - 1e-5; 
      R_max - 1e-5;        
      R_max - 1e-5; 
      0; 
      0;
      0];

% Set nonlinear constraints
nonlcon = @nonlinear

X0_history = [0,0,0]; %Variable to store initial conditions
fval_history = [];
solution = [];
% Loop through Initial Conditions
for x1 = R_min:R_max
    for x2 = R_min:R_max
        if x1 ~= x2
            for x3 = R_min:R_max
                if x3 ~= x1 && x3 ~= x2
                      row_vals = sort([x1, x2, x3]);
                      if ~ismember(row_vals, X0_history, 'rows')
                          [x,fval,exitflag,~] = fmincon(fun,row_vals', A, b, [], [], [],[], nonlcon);
                          if exitflag == 1 || exitflag == 2
                             X0_history = [X0_history; sort([x1, x2, x3])];
                             fval_history = [fval_history; fval];
                             solution = [solution; x'];
                             [X0_history(2:end,:), solution, fval_history]
                          end
                      end

                end
            end
        end
    end
end


%% Choose Charging Price
clc; close all; clear all;
x_min = 5;
x_max = 50;

V_min = .20;
V_max = .50;
% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 
% P ( a_ < a_j < a^)
fun = @(V) -1*(V(1) + V(2) + V(3))

A = [1 -1 0; 
     0 1 -1;
     1 0 -1;
     1 0 0; 
     0 1 0; 
     0 0 1;
     -1 0 0; 
     0 -1 0; 
     0 0 -1;];
  
b = [ -.01; %-.01
      -.02; %-.02
      0;
      V_max - 1e-5; 
      V_max - 1e-5;        
      V_max - 1e-5; 
      0; 
      0;
      0];

% lb = [0;0];
% ub = [R_max; R_max];
nonlcon = @nonlinear2
val = fmincon(fun,[0.27; 0.31; 0.33], A, b, [], [], [],[], nonlcon)

%% Choose Charging Price
clc; close all; clear all;
x_min = 5;
x_max = 50;

V_min = .40;
V_max = .50;
% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 
% P ( a_ < a_j < a^)
fun = @(V) -1*(V(1) + V(2) + V(3))

A = [1 -1 0; 
     0 1 -1;
     1 0 -1;
     1 0 0; 
     0 1 0; 
     0 0 1;
     -1 0 0; 
     0 -1 0; 
     0 0 -1;];
  
b = [ -.05; %-.01
      -.04; %-.02
      0;
      V_max - 1e-5; 
      V_max - 1e-5;        
      V_max - 1e-5; 
      0; 
      0;
      0];

% lb = [0;0];
% ub = [R_max; R_max];
V0_history = [0,0,0]; %Variable to store initial conditions
fval_history = [];
solution = [];
nonlcon = @nonlinear2

for v01 = V_min:.01:V_max
    for v02 = V_min:.01:V_max
        if v01 ~= v02
            for v03 = V_min:.01:V_max
                if v02 ~= v03 && v01 ~= v03
                    [x,fval,exitflag,~] = fmincon(fun,[v01; v02; v03], A, b, [], [], [],[], nonlcon);
                    if exitflag == 1 || exitflag == 2
                        V0_history = [V0_history; sort([v01, v02, v03])];
                        fval_history = [fval_history; fval];
                        solution = [solution; x'];
                        [V0_history(2:end,:), solution, fval_history]
                    end
                end
            end
        end 
    end
end
%% Choose Surge price and Mean charging time
clc; close all; clear all;
x_min = 5;
x_max = 50;
R_max =  30; %kW
% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 
% P ( a_ < a_j < a^)
fun = @(Dw) -1*(1./Dw(1) + Dw(2))

A = [0 -1;
     0 1;
     -1 0;] ;
  
b = [ -x_max/R_max;
       5;
       -.0469];

% lb = [0;0];
% ub = [R_max; R_max];
nonlcon = @nonlinear3
val = fmincon(fun, [1,4], A, b, [], [], [],[], nonlcon)

%% Choose Surge price and Mean charging time
clc; close all; clear all;
x_min = 5;
x_max = 50;
R_max =  40; %kW
om_max = 5;
D_min = 0.0469;
D_max = 3;
% E[x/r] = E[x] * E[1/r] = p1 * 1/R1 + p2 
% P ( a_ < a_j < a^)
fun = @(Dw) -1*(1./Dw(1) + Dw(2))

A = [0 -1;
     0 1;
     -1 0;
     1  0] ;
  
b = [ -x_max/R_max;
       om_max;
       -D_min;
       D_max];
   
   
Dw0_history = [0,0]; %Variable to store initial conditions
fval_history = [];
solution = [];
% lb = [0;0];
% ub = [R_max; R_max];
nonlcon = @nonlinear3
for Dval = 2:.25:D_max
    for wval = 2.5:.5:om_max
        [x,fval,exitflag,~] = fmincon(fun, [Dval; wval], A, b, [], [], [],[], nonlcon);
        if exitflag == 1 || exitflag == 2
            Dw0_history = [Dw0_history; sort([Dval, wval])];
            fval_history = [fval_history; fval];
            solution = [solution; x'];
            [Dw0_history(2:end,:), solution, fval_history]
        end
    end
end
temp = [Dw0_history(2:end,:), solution, fval_history];
[~, idx] = sort(temp(:,5));

final = temp(idx,:)