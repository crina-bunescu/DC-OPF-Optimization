clear; clc; close all;

%% 1. Incarcare date
load('data.mat');
load('sol_cvx.mat');
load('sol_gp.mat');
load('sol_ip.mat');
load('sol_fmincon.mat');

baseMVA_d = double(baseMVA);
Pgmin_d   = double(Pgmin);
Pgmax_d   = double(Pgmax);

% Extragere costuri
cost_cvx = sol_cvx.cost;
cost_gp  = sol_gp.cost;
cost_ip  = sol_ip.cost;
cost_fm  = sol_fmincon.cost;

%% 2. Calcul si tabel

% Calcul Eroare Fezabilitate 
err_cvx = norm(Cg*sol_cvx.Pg - B*sol_cvx.theta - Pd);
err_gp  = norm(Cg*sol_gp.Pg - B*sol_gp.theta - Pd);
err_ip  = norm(Cg*sol_ip.Pg - B*sol_ip.theta - Pd);
err_fm  = norm(Cg*sol_fmincon.Pg - B*sol_fmincon.theta - Pd);

% Calcul Scor R2 fata de referinta CVX
SST = sum((sol_cvx.Pg - mean(sol_cvx.Pg)).^2);
R2_cvx = 1;
R2_gp  = 1 - sum((sol_gp.Pg - sol_cvx.Pg).^2) / SST;
R2_ip  = 1 - sum((sol_ip.Pg - sol_cvx.Pg).^2) / SST;
R2_fm  = 1 - sum((sol_fmincon.Pg - sol_cvx.Pg).^2) / SST;

% Pregatire structuri de date pentru tabel 
Metoda = {'CVX'; 'Grad. Proiectat'; 'Bariera Int.'; 'fmincon'};
Cost_Total_USD = [cost_cvx; cost_gp; cost_ip; cost_fm];
Eroare_Fezabilitate = [err_cvx; err_gp; err_ip; err_fm];
Scor_R2_vs_CVX = [R2_cvx; R2_gp; R2_ip; R2_fm];
Iteratii = {'N/A'; num2str(length(sol_gp.hist_cost)); num2str(length(sol_ip.hist_cost)); num2str(sol_fmincon.iterations)};

% Creare tabel si afisare
T = table(Cost_Total_USD, Eroare_Fezabilitate, Scor_R2_vs_CVX, Iteratii, 'RowNames', Metoda);
fprintf('\n===== REZULTATE COMPARATIVE =====\n\n');
disp(T);

%% 3. FIGURA 1 — Graf retea 
coords = [
    0.00  1.00; 0.20  1.00; 0.05  0.75; 0.20  0.75; 0.35  1.00;
    0.35  0.75; 0.50  1.00; 0.50  0.75; 0.55  0.55; 0.55  0.40;
    0.65  0.55; 0.40  0.55; 0.40  0.40; 0.55  0.65; 0.45  0.65;
    0.35  0.55; 0.45  0.40; 0.30  0.65; 0.25  0.55; 0.40  0.30;
    0.65  0.30; 0.75  0.30; 0.55  0.75; 0.70  0.55; 0.80  0.55;
    0.90  0.40; 0.80  0.40; 0.55  0.85; 0.90  0.30; 0.90  0.20;
];

branch_from = mpc.branch(:,1);
branch_to   = mpc.branch(:,2);
congestion  = abs(sol_cvx.Pflow) ./ Pmax;

figure(1); clf;
set(gcf, 'Name', 'Fig 1 — Retea IEEE 30-bus', 'Position', [50 50 900 650]);

cmap = colormap(hot(256));
colormap(hot);
caxis([0 1]);
hold on;

% Linii
for k = 1:nl
    fr = branch_from(k); to = branch_to(k);
    cong = min(congestion(k), 1.0);
    lw   = 1 + 4 * cong;
    cidx = max(1, round(cong * 255) + 1);
    plot([coords(fr,1) coords(to,1)], [coords(fr,2) coords(to,2)], ...
         '-', 'Color', cmap(cidx, :), 'LineWidth', lw);
end

% Noduri
bus_type = mpc.bus(:,2);
gen_nodes = unique(gen_bus);

for i = 1:nb
    if i == i_slack
        col = [0.13 0.55 0.13]; sz = 120;
    elseif ismember(i, gen_nodes)
        col = [1.0 0.55 0.0];   sz = 100;
    else
        col = [0.2 0.4 0.8];    sz = 60;
    end
    scatter(coords(i,1), coords(i,2), sz, col, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    text(coords(i,1)+0.015, coords(i,2)+0.015, num2str(i), 'FontSize', 7, 'FontWeight', 'bold');
end

h_gen   = scatter(nan, nan, 100, [1.0 0.55 0.0], 'filled', 'MarkerEdgeColor', 'k');
h_load  = scatter(nan, nan, 60,  [0.2 0.4 0.8],  'filled', 'MarkerEdgeColor', 'k');
h_slack = scatter(nan, nan, 120, [0.13 0.55 0.13],'filled', 'MarkerEdgeColor', 'k');

legend([h_gen, h_load, h_slack], {'Generator', 'Sarcina (load)', 'Slack bus'}, ...
       'Location', 'SouthEast', 'FontSize', 10, 'AutoUpdate', 'off');

cb = colorbar('eastoutside');
ylabel(cb, 'Grad congestionare  |P_{flow}| / P_{max}', 'FontSize', 9);
axis off; title('Retea IEEE 30-bus (solutia CVX)', 'FontSize', 13, 'FontWeight', 'bold');
hold off;

%% 4. FIGURA 2 — Alocarea puterilor pe generatoare
figure(2); clf;
set(gcf, 'Name', 'Fig 2 — Alocare Puteri', 'Position', [100 100 900 500]);

data_bar = [sol_cvx.Pg, sol_gp.Pg, sol_ip.Pg, sol_fmincon.Pg] * baseMVA_d; % MW
bar_h = bar(data_bar, 'grouped');
bar_h(1).FaceColor = 'g'; % CVX (Verde)
bar_h(2).FaceColor = 'b'; % GP (Albastru)
bar_h(3).FaceColor = 'r'; % IP (Rosu)
bar_h(4).FaceColor = 'k'; % fmincon (Negru)

hold on;
for i = 1:ng
    x = i + [-0.45 0.45];
    plot(x, [Pgmax_d(i) Pgmax_d(i)] * baseMVA_d, '--r', 'LineWidth', 1.5);
    if Pgmin_d(i) > 0
        plot(x, [Pgmin_d(i) Pgmin_d(i)] * baseMVA_d, '--b', 'LineWidth', 1.5);
    end
end
plot(nan, nan, '--r', 'LineWidth', 1.5, 'DisplayName', 'P_{g,max}');
plot(nan, nan, '--b', 'LineWidth', 1.5, 'DisplayName', 'P_{g,min}');

xticks(1:ng);
xticklabels(arrayfun(@(i) sprintf('G%d (nod %d)', i, gen_bus(i)), 1:ng, 'UniformOutput', false));
ylabel('Putere generata [MW]', 'FontSize', 11);
title('Alocarea Puterilor pe Generatoare (Comparație Metode)', 'FontSize', 13, 'FontWeight', 'bold');
legend([bar_h, plot(nan,nan,'--r'), plot(nan,nan,'--b')], ...
       'CVX', 'Grad. Proiectat', 'Bariera Int.', 'fmincon', 'P_{g,max}', 'P_{g,min}', ...
       'Location', 'NorthEast', 'FontSize', 9);
grid on; box on; hold off;

%% 5. FIGURA 3 — Compararea evolutiilor criteriilor de oprire (Delta Cost)
figure(3); clf;
set(gcf, 'Name', 'Fig 3 — Criterii de oprire', 'Position', [150 150 800 480]);

delta_gp = abs(diff(sol_gp.hist_cost));
delta_ip = abs(diff(sol_ip.hist_cost));

semilogy(1:length(delta_gp), delta_gp, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Grad. Proiectat');
hold on;
semilogy(1:length(delta_ip), delta_ip, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Bariera Int.');

xlabel('Iteratie / Ciclu exterior', 'FontSize', 11);
ylabel('\Delta Cost [$/h]', 'FontSize', 11);
title('Evoluția Criteriului de Oprire (Variația Costului Obiectiv)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'NorthEast', 'FontSize', 10);
grid on; box on; hold off;

%% 6. FIGURA 4 — Compararea convergentei costului obiectiv
figure(4); clf;
set(gcf, 'Name', 'Fig 4 — Convergenta Cost', 'Position', [200 200 800 480]);

plot(sol_gp.hist_cost, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName', 'Grad. Proiectat');
hold on;
plot(sol_ip.hist_cost, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'Bariera Int.');
yline(cost_cvx, 'g--', 'LineWidth', 2, 'DisplayName', 'Referința CVX');
yline(cost_fm, 'k:', 'LineWidth', 2, 'DisplayName', 'Referința fmincon');

xlabel('Iteratie / Ciclu exterior', 'FontSize', 11);
ylabel('Cost Obiectiv [$/h]', 'FontSize', 11);
title('Convergența Costului Obiectiv', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'NorthEast', 'FontSize', 10);
grid on; box on; hold off;