import pandas as pd

H20_MASS = 18.010565

MODS = pd.DataFrame(
    [
        {
            "name": "Acetylation",
            "m/z": 42.0106,
            "symbol_description": "Modification",
            "symbol": "Ac",
            "ion_type": "",
            "type": "Modification",
        },
        {
            "name": "Methylation",
            "m/z": 14.0157,
            "symbol_description": "Modification",
            "symbol": "Me",
            "ion_type": "",
            "type": "Modification",
        },
        {
            "name": "Phosphorylation",
            "m/z": 79.9663,
            "symbol_description": "Modification",
            "symbol": "Ph",
            "ion_type": "",
            "type": "Modification",
        },
        {
            "name": "Sulfonation",
            "m/z": 79.9568,
            "symbol_description": "Modification",
            "symbol": "Su",
            "ion_type": "",
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
    df_mono["ion_type"] = ""
    return df_mono


def make_b_y_df_mono(df_mono: pd.DataFrame) -> pd.DataFrame:
    """
    Split up monosaccharide reference into B and Y ions.
    """
    df_mono_b = df_mono.copy()
    df_mono_y = df_mono.copy()
    df_mono_y["m/z"] = (df_mono_y["m/z"] + H20_MASS).round(2)

    df_mono_y["symbol"] += " ʸ"
    df_mono_y["name"] += " (Y)"
    df_mono_y["ion_type"] = "Y"

    df_mono_b["symbol"] += " ᵇ"
    df_mono_b["name"] += " (B)"
    df_mono_b["ion_type"] = "B"

    df_mono = pd.concat([df_mono_b, df_mono_y], ignore_index=True)

    return df_mono