import pandas as pd 

H20_MASS = 18.010565


def make_df_mono(uploaded_file : str = "mono.json") -> pd.DataFrame:
    """
    Read monosaccharide reference file.
    """
    df_mono = pd.read_json(uploaded_file)
    df_mono = df_mono[df_mono["m/z"] > 0]
    df_mono["m/z"]  = (df_mono["m/z"] - H20_MASS).round(2)

    # df_mono_y = pd.read_json(uploaded_file)
    # df_mono_y = df_mono_y[df_mono_y["m/z"] > 0]
    # df_mono_y["Ion type"] = "Y"
    
    # df_mono_b = df_mono_y.copy()
    # df_mono_b["m/z"]  = (df_mono_b["m/z"] - H20_MASS).round(2)
    # df_mono_y["Ion type"] = "B"
    
    # df_mono = pd.concat([df_mono_b, df_mono_y], ignore_index=True)
    
    return df_mono