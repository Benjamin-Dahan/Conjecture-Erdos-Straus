Conjecture d'Erdos-Straus

Ce projet explore la conjecture d'Erdos-Straus et implémente des algorithmes en C et Python pour tester des hypothèses numériques autour de cette conjecture. Un fichier LaTeX est inclus pour expliquer les concepts et les résultats.

Contenu du projet :

conjecture_erdos_straus.pdf : Un document en LaTeX expliquant la conjecture et les résultats associés.

Test pour n = 10^10 (optimisé) : L'implémentation de l'algorithme en C.

Test pour n = 10^10 (Python) : L'implémentation de l'algorithme en Python.

Installation :
Aucune installation spécifique n'est nécessaire. Télécharge simplement les fichiers et lance les scripts Python et C directement.

Comment utiliser :
Pour exécuter l'algorithme en C, compilez et lancez :
gcc "Test pour n = 10^10 (optimisé).c" -o "Test_n_10_10_optimisé"
./Test_n_10_10_optimisé 10000000000

Pour exécuter l'algorithme en Python, lancez simplement :
python "Test pour n = 10^10 (Python).py" 10000000000

Ouvrez le fichier conjecture_erdos_straus.pdf pour une explication détaillée de la conjecture et des résultats obtenus.

Résultats :
L'algorithme Python est une facon plus abordable de comprendre le code à la première lecture, tandis que l'algorithme C a été optimisé pour les grands nombres (n = 10^10 mais encore plus aussi !).
Attention la valeur de départ n = 3362 car 3361 est la seule exception connue aux programmes sur les 10^10 testées.

Auteurs :
Benjamin Dahan
