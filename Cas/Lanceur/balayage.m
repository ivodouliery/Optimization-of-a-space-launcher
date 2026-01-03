function x_best = balayage(fc, x_init)
% BALAYAGE Effectue une recherche sur grille pour Theta0 et Theta1
% afin de trouver un bon point de depart minimisant les contraintes.

fprintf('\n=== PHASE 0 : Balayage pour trouver les bons Theta0 / Theta1 ===\n');

best_score = inf;
best_x = x_init;

% Plages de recherche
t0_grid = linspace(0, 0.05, 5);  % Proche de la verticale
t1_grid = linspace(0, 0.5, 10);  % 0 a 30 degres env

for t0_val = t0_grid
    for t1_val = t1_grid
        % On teste avec t2, t3 nominaux (0.05)
        x_test = [t0_val; t1_val; 0.05; 0.05];

        [~, c_val] = fc(x_test);

        % Critère : On cherche a minimiser la violation des contraintes
        score = norm(c_val);

        if score < best_score
            best_score = score;
            best_x = x_test;
        end
    end
end

fprintf('Meilleur point de départ trouvé : Theta0=%.4f, Theta1=%.4f (Score=%.2e)\n', best_x(1), best_x(2), best_score);
x_best = best_x;

end
