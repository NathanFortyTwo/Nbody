## Télécom SudParis CSC-5001 : High Performance Systems

Parallelisation de différents algorithmes
(Brute force et Barnes-Hut) pour le calcul de la force gravitationnelle entre N corps.

[Le sujet du projet](https://www-inf.telecom-sudparis.eu/COURS/CSC5001/Supports/Projet/Projet2022/NBody/sujet.php)

Le code source est dans le dossier `src/`
Les fichiers contenant les différentes solutions sont :

- omp_nbody_brute_force.c
- mpi_nbody_brute_force.c
- cuda_nbody_brute_force.cu

Les résultats des temps d'exécution sont dans le dossier `runtimes/`.
Les temps sont calculés pour t=2.
