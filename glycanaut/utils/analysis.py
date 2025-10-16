import pandas as pd
from itertools import combinations
from typing import Tuple, List
import itertools
from . import mono


def prune_isotopes(df: pd.DataFrame, isotope_tol: float = 1.0) -> pd.DataFrame:
    """
    If isotopes within prescribed tolerance exist, retain only that with the smallest m/z.
    """
    charge_states = range(1, 6)
    tol = 1e-3
    df = df.sort_values("m/z").reset_index(drop=True)
    to_drop = []
    for i in range(len(df) - 1):
        j = i + 1
        diffs = []
        while j < len(df) and abs(df.loc[i, "m/z"] - df.loc[j, "m/z"]) < isotope_tol:
            diffs.append(abs(df.loc[i, "m/z"] - df.loc[j, "m/z"]))
            to_drop.append(j)
            j += 1
        matches = [
            z
            for z in charge_states
            if all(abs((1 / z) - diff) <= tol for diff in diffs)
        ]
        df.loc[i, "Charge State"] = 1 if not matches else matches[0]

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


def generate_polysaccharides(
    df_mono: pd.DataFrame, monosaccharide_list: List[str], length: int = 1
) -> pd.DataFrame:
    """
    Generate n-length combinations of monosaccharides in the list provided.
    """
    df_mono.set_index("Name", inplace=True)
    data = []
    polysaccharides = []
    for r in range(2, length + 1):
        polysaccharides.extend(
            itertools.combinations_with_replacement(monosaccharide_list, r)
        )
    for poly_sacc in polysaccharides:
        item = {
            "Name": "",
            "Symbol": "",
            "m/z": 0,
            "Symbol Description": "Polysaccharide",
            "Ion Type": "",
        }
        for mono_sacc in poly_sacc:
            df_mono_sacc = df_mono.loc[mono_sacc]
            item["Name"] += df_mono_sacc.name + " + "
            item["Symbol"] += df_mono_sacc["Symbol"] + " + "
            item["m/z"] += df_mono_sacc["m/z"]
            item["Ion Type"] += df_mono_sacc["Ion Type"]
        
        item["Name"] = item["Name"][:-3]
        item["Symbol"] = item["Symbol"][:-3]
        item["length"] = len(poly_sacc)
        item["Type"] = "Polysaccharide"
        data.append(item)
    
    return pd.DataFrame(data)


def assign_polysaccharides(
    df_diffs: pd.DataFrame,
    df_mono: pd.DataFrame,
    length: int = 1,
    mass_tol: float = 1.0,
):
    """
    Assign peak differences for polysaccharides of length >= 2.
    """
    if length > 1:
        monosaccharide_list = df_diffs[df_diffs["Assigned Mass"] != 0][
            "Assigned"
        ].unique()
        df_poly = generate_polysaccharides(df_mono, monosaccharide_list, length)
        for idx, row in df_diffs.iterrows():
            if row["Assigned Mass"] != 0:
                continue
            diff = row["Peak Difference"]
            if diff >= df_poly["m/z"].min():
                poly_sacc = df_poly[abs(df_poly["m/z"] - diff) < mass_tol]
                if not poly_sacc.empty:
                    poly_row = poly_sacc.iloc[0]
                    df_diffs.loc[idx, "Assigned"] = poly_row["Name"]
                    df_diffs.loc[idx, "Assigned Symbol"] = poly_row["Symbol"]
                    df_diffs.loc[idx, "Assigned Mass"] = poly_row["m/z"]
                    df_diffs.loc[idx, "Ion Type"] = poly_row["Ion Type"]
                    df_diffs.loc[idx, "Length"] = poly_row["length"]
                    df_diffs.loc[idx, "Type"] = poly_row["Type"]
    df_diffs_assigned = df_diffs[df_diffs["Assigned Mass"] != 0].reset_index(drop=True)
    df_diffs_unassigned = df_diffs[df_diffs["Assigned Mass"] == 0].reset_index(
        drop=True
    )
    return df_diffs_assigned, df_diffs_unassigned


def analyse_spectrum(
    df: pd.DataFrame,
    df_mono: pd.DataFrame,
    mass_tol: float,
    length: int = 1,
    use_mods: bool = False,
    use_b_y: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute peak differences for input spectrum and assign wherever possible.
    """
    data = []

    df_mono = df_mono if not use_b_y else mono.make_b_y_df_mono(df_mono)
    df_mono = df_mono if not use_mods else pd.concat([df_mono, mono.MODS])

    for a, b in combinations(df["m/z"], 2):
        diff = abs(a - b)
        diff_dict = {
            "Peak 1": a,
            "Peak 2": b,
            "Peak Difference": diff,
        }
        if diff >= df_mono["m/z"].min():
            mono_sacc = df_mono[abs(df_mono["m/z"] - diff) < mass_tol]
            if not mono_sacc.empty:
                row = mono_sacc.iloc[0]
                assigned = (
                    row["Name"],
                    row["Symbol"],
                    row["m/z"],
                    row["Ion Type"],
                    row["Type"],
                )
            else:
                assigned = ("No match", "", 0, "", "")
        else:
            assigned = ("Too small", "", 0, "", "")
        diff_dict |= {
            "Assigned Mass": assigned[2],
            "Assigned": assigned[0],
            "Assigned Symbol": assigned[1],
            "Ion Type": assigned[3],
            "Type": assigned[4],
            "Length": 1,
        }
        data.append(diff_dict)

    df_diffs = pd.DataFrame(data)

    if not df_diffs.empty:
        df_diffs = df_diffs.sort_values("Assigned Mass", ascending=False)
        df_diffs_assigned, df_diffs_unassigned = assign_polysaccharides(
            df_diffs, df_mono, length, mass_tol
        )
        df_diffs_unassigned = df_diffs_unassigned.drop(
            columns=[
                "Assigned",
                "Assigned Mass",
                "Assigned",
                "Ion Type",
                "Type",
                "Length",
                "Assigned Symbol",
            ]
        )
        df_unmatched = df[
            ~df["m/z"].isin(df_diffs_assigned["Peak 1"])
            & ~df["m/z"].isin(df_diffs_assigned["Peak 2"])
        ].reset_index(drop=True)
        df_unmatched = df_unmatched.drop(columns=["index"])
    else:
        return None, None, None, None

    return df_diffs, df_diffs_assigned, df_diffs_unassigned, df_unmatched


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