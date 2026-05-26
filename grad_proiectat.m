% La fiecare iteratie:
%   1. Calculeaza gradientul: grad = 2*c2.*Pg + c1
%   2. Pas de gradient: Pg_tmp = Pg - alpha * grad
%   3. Proiecteaza pe C (plan-afin) prin cautare binara pe ν
%   4. Rezolva theta din ecuatia de bilant (pentru flux si cost)

clear; clc;
load('data.mat');

%% Parametri
alpha   = 1 / (2 * max(c2));   % pas = 1/L,  L = max eigenvalue of Hessian
eps     = 1e-7;                 % criteriu de oprire
max_it  = 2000;

%% Punct initial fezabil 
Pd_total = sum(Pd);
Pg = Pgmax * Pd_total / sum(Pgmax);

hist_cost = zeros(max_it, 1);
hist_stop = zeros(max_it, 1);

%% Algoritm
for k = 1:max_it

    grad_f = 2 * c2 .* Pg + c1;     % gradient

    Pg_tmp = Pg - alpha * grad_f;   % pas gradient

    Pg_nou = project_affine_box(Pg_tmp, Pgmin, Pgmax, Pd_total);

    %theta si flux
    idx          = setdiff(1:nb, i_slack);
    rhs          = Cg(idx,:) * Pg - Pd(idx);
    theta        = zeros(nb, 1);
    theta(idx)   = B(idx,idx) \ rhs;

    % Inregistreaza iteratia
    hist_cost(k) = sum(c2 .* Pg.^2 + c1 .* Pg) + sum(c0);
    hist_stop(k) = norm(Pg_nou - Pg);

    % Criteriu de oprire
    if hist_stop(k) < eps
        hist_cost = hist_cost(1:k);
        hist_stop = hist_stop(1:k);
        break;
    end

    Pg = Pg_nou;
end

% Theta si flux finale
idx          = setdiff(1:nb, i_slack);
theta        = zeros(nb, 1);
theta(idx)   = B(idx,idx) \ (Cg(idx,:)*Pg - Pd(idx));
Pflow        = Bf * theta;
cost         = sum(c2 .* Pg.^2 + c1 .* Pg) + sum(c0);

% Verifica bilantul
res_bilant = norm(Cg*Pg - B*theta - Pd);

%% Save pentru comparatia finala
sol_gp.Pg        = Pg;
sol_gp.theta     = theta;
sol_gp.cost      = cost;
sol_gp.Pflow     = Pflow;
sol_gp.hist_cost = hist_cost;
sol_gp.hist_stop = hist_stop;
save('sol_gp.mat', 'sol_gp');


%%
function x = project_affine_box(v, lb, ub, s)

    nu_lo = min(v - ub) - 1;
    nu_hi = max(v - lb) + 1;

    for iter = 1:200
        nu    = (nu_lo + nu_hi) / 2;
        x     = max(lb, min(ub, v - nu));
        if sum(x) > s
            nu_lo = nu;
        else
            nu_hi = nu;
        end
        if (nu_hi - nu_lo) < 1e-12
            break;
        end
    end
    x = max(lb, min(ub, v - nu));
end