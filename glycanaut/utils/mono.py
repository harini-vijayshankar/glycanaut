import pandas as pd 

H20_MASS = 18.010565


def make_df_mono(uploaded_file : str = "mono.json") -> pd.DataFrame:
    """
    Read monosaccharide reference file.
    """
    df_mono = pd.read_json(uploaded_file)
    df_mono = df_mono[df_mono["m/z"] > 0]
    df_mono["m/z"]  = (df_mono["m/z"] - H20_MASS).round(2)
    return df_mono