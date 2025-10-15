import streamlit as st
import pandas as pd
import numpy as np
from glycanaut.utils import plotting, analysis, mono

# ----------------------- Layout config ----------------------- #

st.set_page_config(
    page_title="Glycanaut",
    page_icon="🍩",    
    layout="wide"
)

st.markdown(
    """
    <style>
        [data-testid="stSidebar"] {
            min-width: 25%;
            max-width: 25%;
            padding-top: 0rem;
        }
        .block-container {
            padding-left: 2rem;
            padding-right: 2rem;
            padding-top: 2rem;
        }
        [data-testid="stSidebarCollapseButton"] {
            display: none
        
    </style>
""",
    unsafe_allow_html=True,
)

# ----------------------- Sidebar ----------------------- #

st.sidebar.header("Upload spectrum file to begin.")

if "expander_open" not in st.session_state:
    st.session_state.expander_open = False

with st.sidebar.expander("(Optional) Configure SNFG Reference", expanded=False):
    uploaded_mono_file = st.file_uploader(
        "Upload monosaccharide reference file", type=["json"]
    )

uploaded_spectrum_file = st.sidebar.file_uploader("Upload spectrum file", type=["xlsx"])
xls_sheet_name = st.sidebar.empty()

with st.sidebar.expander("Parameters", expanded=st.session_state.expander_open):
    with st.form("process_form"):
        m_z_range_c = st.empty()
        with m_z_range_c.container():
            min_m_z = 0 if "m_z_range_min" not in st.session_state else st.session_state.m_z_range_min
            max_m_z = 5000 if "m_z_range_max" not in st.session_state else st.session_state.m_z_range_max
            m_z_range_val = (min_m_z, max_m_z) if "m_z_range" not in st.session_state else st.session_state.m_z_range
            m_z_range = st.slider(
                    "m/z Range",
                    min_value=min_m_z,
                    max_value=max_m_z,
                    value=m_z_range_val,
                )
        threshold = st.slider(
            "Threshold %", min_value=1, max_value=99, value=15, step=1
        )
        mass_tol = st.number_input(
            "Mass tolerance for ID (Da)",
            value=1.0,
            step=0.00001,
            max_value=100.0,
            min_value=0.00001,
        )
        isotope_tol = st.number_input(
            "Isotope tolerance (Da)",
            value=1.2,
            step=0.00001,
            max_value=100.0,
            min_value=0.00001,
        )
        length = st.number_input(
            "Polysaccharide length to match",
            value=1,
            step=1,
            max_value=8,
            min_value=1,
        )
        use_mods = st.checkbox("Use modifications", value=False)
        use_b_y = st.checkbox("Detect B / Y ions", value=False)
        with st.sidebar:
            submit_button = st.form_submit_button("Analyse")

if uploaded_spectrum_file:
    df = None
    if uploaded_spectrum_file.name.endswith(".csv"):
        df = pd.read_csv(uploaded_spectrum_file)
    elif uploaded_spectrum_file.name.endswith(".tsv"):
        df = pd.read_csv(uploaded_spectrum_file, sep="\t")
    elif uploaded_spectrum_file.name.endswith(".xlsx"):
        with xls_sheet_name.container():
            xls = pd.ExcelFile(uploaded_spectrum_file)
            sheet_name = st.selectbox("Select a sheet", xls.sheet_names)
            if sheet_name:
                df = pd.read_excel(xls, sheet_name=sheet_name)[1:]
            else:
                st.stop()
    else:
        st.error("Unsupported file type.")
        st.stop()

    st.session_state.expander_open = True

    df = df.rename(columns={df.columns[0]: "m/z", df.columns[1]: "Intensity"})
    df = df[pd.to_numeric(df.iloc[:, 0], errors="coerce").notna()]

    with m_z_range_c.container():
        try:
            min_m_z = np.round(df["m/z"].min()).astype(int)
            max_m_z = np.round(df["m/z"].max()).astype(int)
            m_z_range = st.slider(
                "m/z Range",
                min_value=min_m_z,
                max_value=max_m_z,
                value=(min_m_z, max_m_z),
            )
            st.session_state.m_z_range_min = min_m_z
            st.session_state.m_z_range_max = max_m_z
            st.session_state.m_z_range = (min_m_z, max_m_z)
        except:
            pass

# ----------------------- Results ----------------------- #

st.markdown('<h1 style="color: #669673;">Glycanaut</h1>', unsafe_allow_html=True)

if submit_button and not uploaded_spectrum_file:
    st.error("No file uploaded!")

if submit_button:
    top, bottom = st.container(), st.container()
    df_mono = (
        mono.make_df_mono()
        if not uploaded_mono_file
        else mono.make_df_mono(uploaded_mono_file)
    )
    df = analysis.preprocess_data(
        df, threshold=threshold, m_z_range=m_z_range, isotope_tol=isotope_tol
    )
    df_diffs, df_diffs_assigned, df_diffs_unassigned, df_unmatched = analysis.compute_peak_differences(
        df, df_mono, mass_tol=mass_tol, length=length, use_mods=use_mods
    )

    if df_diffs is None:
        st.error("No peaks found. Check the spectrum file / widen the m/z range / reduce the threshold and try again.")
        st.stop()

    with top:
        st.subheader("Mass spectrum")
        st.plotly_chart(plotting.plot_mass_spectrum(df), use_container_width=True)

    with bottom:
        with st.expander("SNFG Reference", expanded=False):
            st.dataframe(df_mono)

        with st.expander("Modifications", expanded=False):
            st.dataframe(mono.MODS.drop(columns=["ion_type"]))
            
        st.subheader("MS2 Identification")

        st.markdown("##### Overview of peak differences")
        st.plotly_chart(plotting.plot_peak_diff_histogram(df_diffs, df_diffs_assigned))

        st.markdown("##### Assigned peak differences")
        st.dataframe(df_diffs_assigned)

        with st.expander("Unassigned peak differences", expanded=False):
            st.dataframe(df_diffs_unassigned)

        with st.expander("Unmatched peaks", expanded=False):
            st.dataframe(df_unmatched)

        st.markdown("##### Peak losses")
        st.write("Little green dots represent peaks. The green path highlighted shows the shortest path between the largest and smallest peaks.")
        st.plotly_chart(
            plotting.plot_peak_diff_graph(df_diffs_assigned), use_container_width=True
        )

else:
    st.info("Upload a spectrum file and click Analyse to begin.")