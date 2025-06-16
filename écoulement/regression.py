import pandas as pd 
import statsmodels.api as sm
from typing import Union
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller, kpss
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.stats.diagnostic import acorr_ljungbox
from scipy.stats import jarque_bera
from statsmodels.stats.api import het_breuschpagan

class Regressor:

    def __init__(self, endog: pd.Series, exog: Union[pd.Series, pd.DataFrame], 
                add_constant: bool = True) -> None:
        self.X = exog 
        self.X = sm.add_constant(self.X) if add_constant else self.X
        self.y = endog

    def linear_regression(self) :
        self.linear_model = sm.OLS(self.y, self.X).fit()
        return self.linear_model
    
    @staticmethod
    def adf_test(variable: pd.Series, sequential_strategy: bool = True) -> bool:
        """
        Fonction qui affiche les résultats des tests ADF + renvoie True si non stationnaire. 
        """
        if not sequential_strategy:
            pvalue = adfuller(variable)[1]
            print(f"P-value du test de Dickey-Fuller : {pvalue:.4f}")
            return True if pvalue > 0.05 else False
        
        adf_trend = adfuller(variable, regression='ct') 
        print(f"P-value du test de Dickey-Fuller avec tendance et constante : {adf_trend[1]:.4f}")
        if adf_trend[1] < 0.05:
            print("La série est stationnaire autour d'une tendance")
            return False

        adf_drift = adfuller(variable, regression='c') 
        print(f"P-value du test de Dickey-Fuller avec constante : {adf_drift[1]:.4f}")
        if adf_drift[1] < 0.05:
            print("La série est stationnaire autour d'une constante")
            return False

        adf_none = adfuller(variable, regression='n') 
        print(f"P-value du test de Dickey-Fuller sans tendance ni constante : {adf_none[1]:.4f}")
        if adf_none[1] < 0.05:
            print("La série est stationnaire")
            return False
        
        print("La série est non stationnaire")
        return True
        
    @staticmethod
    def kpss_test(variable: pd.Series, sequential_strategy: bool = True) -> bool:
        """
        Fonction qui affiche les résultats des tests KPSS + renvoie True si non stationnaire. 
        """
        if not sequential_strategy:
            statistic, pvalue, _, _ = kpss(variable)
            print(f"P-value du test KPSS : {pvalue:.4f}")
            return True if pvalue < 0.05 else False
        
        statistic_trend, pvalue_trend, _, _ = kpss(variable, regression='ct') 
        print(f"P-value du test KPSS avec tendance et constante : {pvalue_trend:.4f}")
        if pvalue_trend > 0.05:
            print("La série est stationnaire autour d'une tendance")
            return False

        statistic_drift, pvalue_drift, _, _ = kpss(variable, regression='c') 
        print(f"P-value du test KPSS avec constante : {pvalue_drift:.4f}")
        if pvalue_drift > 0.05:
            print("La série est stationnaire autour d'une constante")
            return False
        
        print("La série est non stationnaire")
        return True
    
    @staticmethod
    def checkresiduals(model,ljb_lags: list = [10]) -> None:
        """
        Fonction qui effectue les tests usuels sur les résidus : 
            - test de Ljung Box d'absence d'autocorrélation
            - test de normalité de Jarque-Bera
            - tesr d'homoscédasticité de Breusch-Pagan
        """
        residuals = model.resid
        ljb_test = acorr_ljungbox(residuals, lags = ljb_lags)
        print("============= TEST D'ABSENCE D'AUTOCORRELATION (LJUNG-BOX) =============")
        print(f"P-value du test de Ljung-Box : {ljb_test["lb_pvalue"].iloc[0]}")
        if ljb_test["lb_pvalue"].iloc[0] > 0.05:
            print("Les résidus sont non autocorrélés")
        else:
            print("Attention ! Les résidus présentent de l'autocorrélation")
        
        jb_test = jarque_bera(residuals)
        print("=============TEST DE NORMALITE (JARQUE-BERA) =============")
        print(f"P-value du test de Jarque-Bera : {jb_test.pvalue}")
        print(f"Statistique du test de Jarque-Bera : {jb_test.statistic}")
        if jb_test.pvalue < 0.05:
            print("Les résidus sont non gaussiens")
        else:
            print("Les résidus sont gaussiens")
        
        bp_test = het_breuschpagan(residuals, model.model.exog)
        labels = ["LM stat", "LM p-value", "F-stat", "F-pvalue"]
        print("=============TEST D'HOMOSCEDASTICITE (BREUSCH-PAGAN) =============")
        for l, val in zip(labels, bp_test):
            print(f"{l}: {val}")
        
        plot_acf(residuals, lags = max(ljb_lags), 
                title = "Autocorrélation des résidus")
        plt.show()


    
