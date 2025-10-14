import pandas as pd
from itertools import combinations
from typing import Tuple


def prune_isotopes(df: pd.DataFrame, isotope_tol: float = 1.0) -> pd.DataFrame:
    """
    If isotopes within prescribed tolerance exist, retain only that with the smallest m/z.
    """
    df = df.sort_values("m/z").reset_index(drop=True)
    to_drop = []
    for i in range(len(df) - 1):
        j = i + 1
        while j < len(df) and abs(df.loc[i, "m/z"] - df.loc[j, "m/z"]) < isotope_tol:
            to_drop.append(j)
            j += 1
    return df.drop(to_drop).reset_index(drop=True)


def threshold_peaks(
    df: pd.DataFrame, threshold: float = 0.10, m_z_range: Tuple[int, int] = (0, 5000)
) -> pd.DataFrame:
    """
    Discard peaks with intensity below a specified % of the tallest peak.
    """
    max_intensity = df["Intensity"].max()
    df = df[df["Intensity"] > threshold * 0.01 * max_intensity].reset_index()
    df = df[(df["m/z"] >= m_z_range[0]) & (df["m/z"] <= m_z_range[1])]
    return df


def preprocess_data(
    df: pd.DataFrame, threshold: float, m_z_range: Tuple[int, int], isotope_tol=1.0
) -> pd.DataFrame:
    """
    Threshold and prune isotopes of input spectrum.
    """
    df = threshold_peaks(df, threshold=threshold, m_z_range=m_z_range)
    df = prune_isotopes(df, isotope_tol=isotope_tol)
    return df


def compute_peak_differences(
    df: pd.DataFrame, df_mono: pd.DataFrame, mass_tol: float
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute peak differences for input spectrum and assign wherever possible.
    """
    data = []
    for a, b in combinations(df["m/z"], 2):
        diff = abs(a - b)
        match = df_mono[abs(df_mono["m/z"] - diff) < mass_tol]
        if not match.empty:
            row = match.iloc[0]  # just takes first. hmmm
            data.append(
                {
                    "Peak 1": a,
                    "Peak 2": b,
                    "Peak Difference": diff,
                    "Assigned": row["name"],
                    "Assigned Symbol": row["symbol"],
                    "Assigned Mass": row["m/z"],
                }
            )
    df_diffs = pd.DataFrame(data).sort_values("Assigned Mass", ascending=False)
    df_diffs_assigned = df_diffs[df_diffs["Assigned Mass"] != 0].reset_index(drop=True)
    
    df_unmatched = df[
        ~df["m/z"].isin(df_diffs_assigned["Peak 1"])
        & ~df["m/z"].isin(df_diffs_assigned["Peak 2"])
    ].reset_index(drop=True)
    df_unmatched = df_unmatched.drop(columns=["index"])
    
    return df_diffs, df_diffs_assigned, df_unmatched


def all_combinations(weights, target, idx=0, current=None):
    """
    Find all combinations of weights that exactly add up to a target.
    """
    if current is None:
        current = []

    if target == 0:
        yield current
        return
    if idx >= len(weights):
        return

    max_count = target // weights[idx]
    for k in range(max_count + 1):
        yield from all_combinations(
            weights, target - k * weights[idx], idx + 1, current + [k]
        )