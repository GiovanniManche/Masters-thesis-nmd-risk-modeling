"""
Ensemble de fonctions qui servent dans le Notebook principal
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from écoulement.regression import Regressor

def lorenz_curve_and_gini(outstandings: np.ndarray) -> float:
    """
    Fonction qui, à partir des données granulaires d'une quinzaine considérée,
    trace la courbe de Lorenz et calcule le coefficient de Gini.
    """
    # Tri par ordre décroissant car on veut voir à quel point une petite part des comptes détient les encours
    sorted_outstandings = -np.sort(-outstandings)
    cumulative = np.cumsum(sorted_outstandings)
    cumulative = np.insert(cumulative, 0, 0)
    cumulative = cumulative/cumulative[-1]
    x = np.linspace(0,1, len(cumulative))

    # On applique un moins devant le coefficient de Gini car tri décroissant
    gini_coefficient = - (1 - 2* np.trapz(cumulative, x))

    plt.figure(figsize=(10,6))
    plt.plot(x, cumulative, label = "Courbe de Lorenz")
    plt.plot([0,1], [0,1], color = 'red', label = "Egalité parfaite")
    plt.title("Courbe de Lorenz - Analyse graphique de la concentration")
    plt.xlabel("Part cumulée du nombre de compte")
    plt.ylabel("Part cumulée des encours")
    plt.axvline(x=0.05, color = "black", linestyle = "--", label = "5%")
    plt.axvline(x=0.01, color = "darkgreen", linestyle = "--", label = "1%")
    plt.legend()
    plt.text(0.75, 0.5, f"Indice de Gini : {gini_coefficient:2f}",
             bbox=dict(facecolor='red', alpha=0.5))

    plt.show()
    return gini_coefficient

def detect_additive_or_multiplicative(considered_series: pd.Series, period: int = 12) -> str:
    """
    Fonction qui détecte si une série est additive ou multiplicative
    L'intuition est qu'une série avec cycle saisonnier multiplicatif verra 
    l'amplitude, et donc la variance, de ce cycle augmenter avec le niveau moyen de la série. 
    """
    if not isinstance(considered_series.index, pd.DatetimeIndex):
        raise ValueError("L'index de la série doit être un DatetimeIndex.")
        
    # On crée les groupes saisonniers en fonction de la taille du cycle
    cycle = considered_series.index.to_period(f'{int(12/period)}M').strftime('%m') if period == 12 else (
        considered_series.index.to_series().dt.to_period(f'{int(12/period)}M')).dt.month % period
        
    # On regroupe les encours par cycle
    df = pd.DataFrame({'value': considered_series.values, 'cycle': cycle})
    grouped = df.groupby('cycle')['value']

    # On calcule les moyennes et écarts-types par cycle
    means = grouped.mean()
    stds = grouped.std()

    # Régression linéaire des écarts-types sur les moyennes 
    model = Regressor(endog = stds, exog = means, add_constant=True).linear_regression()
    print(model.summary())
    pvalue = model.pvalues[1]
    return 'multiplicative' if pvalue < 0.05 else 'additive'

def histogram(x_axis: list, y_axis: np.ndarray, 
              title: str = "Histogramme") -> None:
    plt.figure(figsize = (10,6))
    plt.bar(x_axis, y_axis)
    plt.title(title)
    plt.ylabel("Valeur")
    plt.tight_layout()
    plt.show()


