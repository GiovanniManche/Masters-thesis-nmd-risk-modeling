import numpy as np
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt


def forecast_taux_client(EUR12M_forecast, r_client_last, r_ref_last, xi, gamma, c, cLT, alpha):
    forecast = []
    r_client_prev = r_client_last
    r_ref_prev = r_ref_last

    for t in range(len(EUR12M_forecast)):
        delta_ref = EUR12M_forecast[t] - r_ref_prev
        res_eq = r_client_prev - cLT - alpha * r_ref_prev
        delta_client = c + xi * delta_ref + gamma * res_eq

        r_client_new = r_client_prev + delta_client
        forecast.append(r_client_new)

        # update
        r_ref_prev = EUR12M_forecast[t]
        r_client_prev = r_client_new

    return np.array(forecast)


def compute_tvmm_dynamic_from_curve(df_curve, replacement_ratio=1/72, duration_months=72):
    """
    Calcule le TVMM dynamique à encours constant en identifiant les colonnes de taux correspondant à +6 ans.

    Parameters:
    - df_curve : DataFrame avec colonnes comme échéances (dates) et lignes comme dates projetées
    - replacement_ratio : fraction de l'encours remplacée chaque mois (ex : 1/72 pour 6 ans)
    - duration_years : durée de placement des strates (ex: 6 ans)

    Returns:
    - pd.Series des TVMM calculés
    """

    # Assurer que la colonne 'Date' est bien au format datetime
    df_curve = df_curve.copy()
    df_curve['Date'] = pd.to_datetime(df_curve['Date'])

    # Convertir les colonnes (hors 'Date') en datetime pour les échéances
    curve_cols = [col for col in df_curve.columns if col != 'Date']
    maturity_dates = pd.to_datetime(curve_cols, dayfirst=True)
    
    tvmm_series = []
    previous_tvmm = None

    for idx, row in df_curve.iterrows():
        t_date = row['Date']
        # Date cible = t + 6 ans
        target_date = t_date + pd.DateOffset(months=duration_months)

        # Trouver la colonne la plus proche de t + 6 ans
        diffs = np.abs(maturity_dates - target_date)
        closest_col_idx = np.argmin(diffs)
        closest_col_name = curve_cols[closest_col_idx]

        taux_replacement = row[closest_col_name]

        if previous_tvmm is None:
            tvmm = taux_replacement  # Initialisation à t0
        else:
            tvmm = (1 - replacement_ratio) * previous_tvmm + replacement_ratio * taux_replacement

        tvmm_series.append(tvmm)
        previous_tvmm = tvmm

    return pd.Series(tvmm_series, index=df_curve['Date'], name="TVMM")

def compute_aggregated_tvmm(tvmm_strates, weights):
    """
    Calcule le TVMM agrégé pondéré pour chaque scénario, selon les poids fournis.
    """
    scenarios = tvmm_strates[0].keys()
    tvmm_total = {}

    for scenario in scenarios:
        weighted_sum = sum(weights[i] * tvmm_strates[i][scenario] for i in range(len(weights)))
        tvmm_total[scenario] = weighted_sum

    return tvmm_total

def compute_mni_dict(tvmm_model_dict, tx_client_forecasts):
    """
    Calcule la MNI agrégée pour un modèle donné (sur tous les scénarios).
    
    Paramètres :
        tvmm_model_dict : dict de séries TVMM (clé = nom scénario)
        tx_client_forecasts : dict de séries ou arrays taux client (même clés)

    Retourne :
        dict des séries de MNI
    """
    mni_model_dict = {
        name: tvmm_model_dict[name] - np.array(tx_client_forecasts[name])
        for name in tvmm_model_dict
    }
    return mni_model_dict

def plot_mni_all_models_by_scenario(mni_models,  tx_client_forecasts, tvmm_average_model):
    scenarios = next(iter(mni_models.values())).keys()  # liste des scénarios

    for scenario in scenarios:
        plt.figure(figsize=(12, 6))

        # Tracer les MNI de chaque modèle
        for model_name, mni_dict in mni_models.items():
            plt.plot(mni_dict[scenario], label=f"MNI - {model_name}", lw=2)

        # Taux client
        if scenario in tx_client_forecasts:
            reference_model = next(iter(mni_models))
            reference_index = mni_models[reference_model][scenario].index
            tx_client_series = pd.Series(tx_client_forecasts[scenario], index=reference_index)
            plt.plot(tx_client_series, label="Taux client (ECM)", linestyle="--", color="black")


        # TVMM moyen (modèle estimé par exemple)
        if scenario in tvmm_average_model:
            plt.plot(tvmm_average_model[scenario], label="TVMM estimé", linestyle="--", color="gray")

        plt.title(f"MNI comparée — Scénario : {scenario}", fontsize=14)
        plt.xlabel("Date")
        plt.ylabel("Taux (%)")
        plt.axhline(0, color="black", linestyle="--", linewidth=1)
        plt.grid(True, linestyle="--", alpha=0.3)
        plt.legend(title="Modèle d'écoulement")
        plt.tight_layout()
        plt.show()
