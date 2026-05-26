clear; clc;
load('data.mat');

Pgmin = double(Pgmin);
Pgmax = double(Pgmax);
c2    = double(c2);
c1    = double(c1);
c0    = double(c0);
Pd    = double(Pd);
Pmax  = double(Pmax);

% Dimensiuni vector de optimizare
n_vars = ng + nb;

%% Functia obiectiv
obj_fun = @(x) sum(c2 .* x(1:ng).^2 + c1 .* x(1:ng)) + sum(c0);

%% Constrangeri de inegalitate 
A = [zeros(nl, ng),  full(Bf);
     zeros(nl, ng), -full(Bf)];
b = [Pmax; Pmax];

%% Constrangeri de egalitate 
% 1. Bilant putere in fiecare nod
Aeq_bal = [full(Cg), -full(B)];
beq_bal = Pd;

% 2. Fixare nod de referinta 
Aeq_slack = zeros(1, n_vars);
Aeq_slack(ng + i_slack) = 1;
beq_slack = 0;

Aeq = [Aeq_bal; Aeq_slack];
beq = [beq_bal; beq_slack];

%% Margini variabile 
lb = [Pgmin; -inf(nb, 1)];
ub = [Pgmax;  inf(nb, 1)];

%% Punct initial fezabil
Pg_init = Pgmax * sum(Pd) / sum(Pgmax);
theta_init = zeros(nb, 1);
x0 = [Pg_init; theta_init];

%% Configurare algoritm
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'OptimalityTolerance', 1e-7, ...
    'ConstraintTolerance', 1e-7);

%% Rulare 
tic;
[x_opt, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, [], options);
t_exec = toc;

%% Extragere si salvare date
Pg = x_opt(1:ng);
theta = x_opt(ng+1:end);
cost = fval;
Pflow = Bf * theta;

sol_fmincon.Pg = Pg;
sol_fmincon.theta = theta;
sol_fmincon.cost = cost;
sol_fmincon.Pflow = Pflow;
sol_fmincon.t_exec = t_exec;
sol_fmincon.iterations = output.iterations;

save('sol_fmincon.mat', 'sol_fmincon');
