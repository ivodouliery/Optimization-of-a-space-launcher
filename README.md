# Optimization-of-a-space-launcher

Optimisation d'un lanceur spatial 3 étages pour mise en orbite à 200 km.

## Utilisation rapide

```matlab
% 1. Initialiser les chemins
setup_project

% 2. Résoudre le problème d'étagement pour Vp donné
cd Cas/PE
PE_Newton      % Solution analytique
PE_SQP         % Solution numérique (validation)

% 3. Simuler et optimiser la trajectoire
cd ../Lanceur
lanceur        % Lance l'optimisation complète
```

## Structure

| Dossier | Description |
|---------|-------------|
| `SQP/` | Optimiseur SQP générique |
| `Cas/PE/` | Problème d'étagement (Newton + SQP) |
| `Cas/Lanceur/` | Simulateur de trajectoire + optimisation |
| `Cas/MHW4D/` | Cas test validation SQP |
| `Cas/Ariane1_Masse/` | Cas test Ariane 1 |

## Paramètres du lanceur

Modifier dans `lanceur.m` :
- `Vp` : vitesse propulsive (m/s)
- `mu` : charge utile (kg)
- `Hc` : altitude cible (m)

## Résultat

Lanceur optimal : **M₀ ≈ 35 tonnes** pour 1000 kg en orbite à 200 km.