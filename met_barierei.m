clear; clc;
load('data.mat');

Pgmin = double(Pgmin);
Pgmax = double(Pgmax);
c2    = double(c2);
c1    = double(c1);
c0    = double(c0);
Pd    = double(Pd);
Pmax  = double(Pmax);

%% Parametri
mu      = 5.0;    % parametru initial 
gamma   = 0.15;   % factor de reducere mu
n_outer = 50;     % cicluri exterioare max
n_inner = 80;     % pasi Newton maximi per ciclu
tol_in  = 1e-8;   % criteriu oprire interior
tol_out = 1e-8;   % criteriu oprire exterior: variatie cost intre cicluri
eps_bd  = 1e-9;   % margine minima stricta fata de frontiera

%% Precalculare M si d_f (inainte de bucle)
idx    = setdiff(1:nb, i_slack); 
B_red  = full(B(idx, idx));     
Cg_red = full(Cg(idx, :));         
Pd_red = Pd(idx);                  
Bf_red = full(Bf(:, idx));        

% Factorizare LU a lui B_red
[L_B, U_B, P_B] = lu(B_red);


BinvCg = U_B \ (L_B \ (P_B * Cg_red));   
BinvPd = U_B \ (L_B \ (P_B * Pd_red));   
M      = Bf_red * BinvCg;                 
d_f    = Bf_red * BinvPd;                 

Pd_total = sum(Pd);
e        = ones(ng, 1);

% Punct initial (alocare proportionala cu capacitatea maxima)
Pg = Pgmax * (Pd_total / sum(Pgmax));

% Metoda barierei necesita un punct strict in interiorul limitelor fizice
Pg = max(Pgmin + eps_bd, min(Pgmax - eps_bd, Pg));
Pg = Pg * (Pd_total / sum(Pg));   % refacem bilantul puterilor

% Initializare multiplicator Lagrange: estimare din conditia de stationaritate
nu = -mean(2*c2.*Pg + c1);

hist_cost    = zeros(n_outer, 1);
hist_mu      = zeros(n_outer, 1);
n_outer_used = n_outer;

%% Bucla exterioara: reducere mu 
for outer = 1:n_outer

    % Bucla interioara: pasi Newton la mu fix 
    for inner = 1:n_inner

        Pf = M * Pg - d_f;
        g_obj  = 2*c2.*Pg + c1;                                   % grad al costului de generare
        g_box  = -mu * (1./(Pg - Pgmin) - 1./(Pgmax - Pg));       % grad al barierei pentru Pgmin < Pg < Pgmax
        g_line =  mu * (M' * (1./(Pmax - Pf) - 1./(Pmax + Pf)));  % grad al barierei pentru -Pmax < Pf < Pmax

        % Reziduu stationaritate si fezabilitate
        r_stat = g_obj + g_box + g_line + nu*e;
        r_feas = sum(Pg) - Pd_total;

        % Criteriu oprire interior
        if norm([r_stat; r_feas]) < tol_in
            break;
        end

        % Hessiana functiei
        H_obj  = 2 * diag(c2); % hessiana cost
        H_box  = mu * diag(1./(Pg - Pgmin).^2 + 1./(Pgmax - Pg).^2); % hessiana barierei box
        D_line = diag(1./(Pmax - Pf).^2 + 1./(Pmax + Pf).^2);  % nl x nl
        H_line = mu * (M' * D_line * M); % hessiana barierei linii
        H_mu   = H_obj + H_box + H_line;   % ng x ng, PD garantat in interior

        % Sistem KKT
        KKT = [H_mu, e; e', 0];
        sol = KKT \ [-r_stat; -r_feas];
        dPg = sol(1:ng);
        dnu = sol(ng+1);

        % Line search:
        % Gaseste alpha_max = cel mai mare pas care mentine TOATE
        % barierele strict pozitive.
        % Principiu: daca xi > 0 si dxi < 0, pasul maxim = xi / (-dxi).
        dPf       = M * dPg;
        alpha_max = 1.0;

        % Limita din barierele box
        mask = dPg < 0;
        if any(mask)
            alpha_max = min(alpha_max, ...
                min((Pg(mask) - Pgmin(mask)) ./ (-dPg(mask))));
        end
        mask = dPg > 0;
        if any(mask)
            alpha_max = min(alpha_max, ...
                min((Pgmax(mask) - Pg(mask)) ./ dPg(mask)));
        end

        % Limita din barierele de flux pe linii
        mask = dPf > 0;
        if any(mask)
            alpha_max = min(alpha_max, ...
                min((Pmax(mask) - Pf(mask)) ./ dPf(mask)));
        end
        mask = dPf < 0;
        if any(mask)
            alpha_max = min(alpha_max, ...
                min((Pmax(mask) + Pf(mask)) ./ (-dPf(mask))));
        end

        % Factor de siguranta 0.99 (evita atingerea exacta a frontierei)
        step = min(1.0, 0.99 * alpha_max);

        Pg = Pg + step * dPg;
        nu = nu + step * dnu;

    end % bucla interioara

    hist_cost(outer) = sum(c2.*Pg.^2 + c1.*Pg) + sum(c0);
    hist_mu(outer)   = mu;

    % Criteriu de oprire exterior
    if outer > 2 && abs(hist_cost(outer) - hist_cost(outer-1)) < tol_out
        n_outer_used = outer;
        hist_cost    = hist_cost(1:outer);
        hist_mu      = hist_mu(1:outer);
        break;
    end

    mu = mu * gamma;

end % bucla exterioara

%% Theta si flux final
theta        = zeros(nb, 1);
theta(idx)   = U_B \ (L_B \ (P_B * (Cg_red*Pg - Pd_red)));
Pflow        = Bf * theta;
cost         = sum(c2.*Pg.^2 + c1.*Pg) + sum(c0);

%% Metrici de calitate
res_bilant  = norm(Cg*Pg - B*theta - Pd);
viol_box    = max([Pgmin - Pg; Pg - Pgmax]);
viol_lines  = max(abs(Pflow) - Pmax);

%% Save pentru comparatia finala
sol_ip.Pg        = Pg;
sol_ip.theta     = theta;
sol_ip.cost      = cost;
sol_ip.Pflow     = Pflow;
sol_ip.hist_cost = hist_cost;
sol_ip.hist_mu   = hist_mu;
save('sol_ip.mat', 'sol_ip');
fprintf('\nSalvat in: sol_ip.mat\n');