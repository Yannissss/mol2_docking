# Partie - On danse ensemble ?

Le code présent en l'état est incomplet. Vous trouverez ci-dessous dans les commentaires la manière dont nous aurions procédé pour aboutir à un code fonctionnel si nous n'avions pas été bloqués par la partie évaluation.

## Répartition du travail

Hafsa
- Function `sample`
- Function 'eval_solution'
- Function 'check_hydrogen_bound'

Yannis
- Function 'generate_atom'
- Function 'compute_angle'

## Commentaires

Seule la partie génération des solutions initiales est fonctionnelle. Pour cela, notre code génère des positions aléatoires pour les atomes du ligand tout en s'assurant que leurs positions n'entraînent pas une liaison covalente avec un des atomes du site.

La partie évaluation d'une solution est celle qui nous a posé problème. Seulement son squelette est disponible dans le code.

Concernant les parties non abordées, nous avions quand même réfléchi comment nous nous y prendre pour les implémenter.

Pour la partie sélection des parents, nous comptions tirer aléatoirement des couples parmi les X% meilleures solutions. X étant un hyper-paramètre que nous aurions essayé de choisir empiriquement. Le reste de la population aurait été régénéré aléatoire pour l'itération suivante.

Pour la partie recombinaison des individus, nous aurions créé une solution fille entre prenant aléatoirement des positions des atomes dans le ligand parmi celles dans ses parents.

Pour la partie mutation des individus, nous aurions tiré aléatoirement les positions d'un faible nombre d'atomes de la solution existante.

Concernant la parallélisation, l'évaluation, la mutation et la génération des solutions initiales peuvent être parallélisées. La partie génération des solutions est la plus simple à paralléliser, mais rapporte peu de gain en performances. En revanche, la parallélisation de l'évaluation aurait été une piste intéressante pour obtenir des gains de performances à chaque itération de l'algorithme.