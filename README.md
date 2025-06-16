# Mémoire
Afin de quantifier l'impact du modèle, des taux et de l'interaction entre les deux sur les marges afférentes aux produits non échéancés dans une optique de gestion de risque de taux en ALM, nous avons dû :
- modéliser l'écoulement des produits non échéancés (qui n'ont pas de maturité contractuelle mais on établit une maturité comportementale sur la base d'une étude économétrique),
- modéliser les projections de taux d'intérêt sur différents scénarios (standard, choc parallèle de la courbe,...). Pour cela, nous nous sommes appuyés sur le *dynamic Nelson-Siegel model*, proposant une approche simple (VAR/AR, Diebold & Li (2006)) puis plus sophistiquée (VECM, TVP-VAR (Koop & Korobilis, 2012/2013)).
- réunir ces éléments pour calculer un proxy de marge nette d'intérêt, et en étudier les comportements face aux facteurs explicatifs taux, modèles et interactions.
