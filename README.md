# Modèle de transmission de *Clostridioides difficile*

Ce dépôt contient le code R utilisé dans le cadre de mon mémoire pour développer, calibrer et analyser un modèle compartimental déterministe de transmission de *Clostridioides difficile* entre l’hôpital et la communauté.

## Contenu du dépôt

- `0_packages.R` : chargement des packages nécessaires
- `1_modele.R` : définition du modèle compartimental et des sorties calculées
- `2_calibration.R` : fonctions de calibration du modèle
- `3_scenarios.R` : simulation des scénarios d’intervention
- `4_analyse_prcc.R` : analyse de sensibilité par PRCC
- `plot.R` : fonctions de visualisation des résultats
- `main.R` : script principal permettant de lancer l’ensemble de l’analyse

## Utilisation

Pour lancer l’analyse complète, ouvrir le fichier `main.R` dans RStudio puis exécuter le script.

Les packages nécessaires sont chargés dans le fichier :

```r
source("0_packages.R")
