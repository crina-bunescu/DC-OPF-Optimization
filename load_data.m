clear; clc;

mpc = case30;
baseMVA = mpc.baseMVA;   % baza de normalizare

% Dimensiuni
nb = size(mpc.bus, 1);      % numar noduri  
ng = size(mpc.gen, 1);      % numar generatoare
nl = size(mpc.branch, 1);   % numar linii

% Date noduri
Pd   = mpc.bus(:, 3) / baseMVA;   % consum activ
Vmax = mpc.bus(:, 12);            % limita sup tensiune
Vmin = mpc.bus(:, 13);            % limita inf tensiune

i_slack = find(mpc.bus(:, 2) == 3);   % indexul nodului de referinta (=1)

% Date generatoare (mpc.gen)
Pg0   = mpc.gen(:, 2)  / baseMVA;   % putere initiala
Pgmax = mpc.gen(:, 9)  / baseMVA;   % limita superioara generator
Pgmin = mpc.gen(:, 10) / baseMVA;   % limita inferioara generator

% Nodul la care e conectat fiecare generator
gen_bus = mpc.gen(:, 1);

% Cg(i,j)=1 daca gen j e pe nodul i
Cg = sparse(gen_bus, (1:ng)', 1, nb, ng);

% Coeficienti functie cost 
c2 = mpc.gencost(:, 5) * baseMVA^2;
c1 = mpc.gencost(:, 6) * baseMVA;
c0 = mpc.gencost(:, 7);

% Date linii
Pmax = mpc.branch(:, 6) / baseMVA;    % capacitate termica a fiecarei linii
Pmax(Pmax == 0) = 9999;               % 0 in matpower = nelimitat

% B  = susceptanta nodala
% Bf = flux pe linii
[B, Bf, ~, ~] = makeBdc(baseMVA, mpc.bus, mpc.branch);


save('data.mat', ...
     'mpc', 'baseMVA', 'nb', 'ng', 'nl', ...
     'B', 'Bf', 'Cg', 'gen_bus', ...
     'Pd', 'Pgmin', 'Pgmax', 'Pg0', 'Pmax', ...
     'c2', 'c1', 'c0', 'i_slack', 'Vmax', 'Vmin');
