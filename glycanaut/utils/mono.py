import pandas as pd

H20_MASS = 18.010565

MODS = pd.DataFrame(
    [
        {
            "Name": "Acetylation",
            "m/z": 42.0106,
            "Symbol Description": "Modification",
            "Symbol": "Ac",
            "Ion Type": "",
            "type": "Modification",
        },
        {
            "Name": "Methylation",
            "m/z": 14.0157,
            "Symbol Description": "Modification",
            "Symbol": "Me",
            "Ion Type": "",
            "type": "Modification",
        },
        {
            "Name": "Phosphorylation",
            "m/z": 79.9663,
            "Symbol Description": "Modification",
            "Symbol": "Ph",
            "Ion Type": "",
            "type": "Modification",
        },
        {
            "Name": "Sulfonation",
            "m/z": 79.9568,
            "Symbol Description": "Modification",
            "Symbol": "Su",
            "Ion Type": "",
            "type": "Modification",
        },
    ]
)


def make_df_mono(uploaded_file: str = "mono.json") -> pd.DataFrame:
    """
    Read monosaccharide reference file.
    """
    df_mono = pd.read_json(uploaded_file)
    df_mono = df_mono[df_mono["m/z"] > 0]
    df_mono["m/z"] = (df_mono["m/z"] - H20_MASS).round(2)
    df_mono["Ion Type"] = ""
    return df_mono


def make_b_y_df_mono(df_mono: pd.DataFrame) -> pd.DataFrame:
    """
    Split up monosaccharide reference into B and Y ions.
    """
    df_mono_b = df_mono.copy()
    df_mono_y = df_mono.copy()
    df_mono_y["m/z"] = (df_mono_y["m/z"] + H20_MASS).round(2)

    df_mono_y["Symbol"] += " ʸ"
    df_mono_y["Name"] += " (Y)"
    df_mono_y["Ion Type"] = "Y"

    df_mono_b["Symbol"] += " ᵇ"
    df_mono_b["Name"] += " (B)"
    df_mono_b["Ion Type"] = "B"

    df_mono = pd.concat([df_mono_b, df_mono_y], ignore_index=True)

    return df_mono