"""Review coverage sensitivity in four historical donor categories."""
from pathlib import Path

import pandas as pd
import streamlit as st

PANEL = Path(__file__).parent / "results" / "scenarios.csv"

st.set_page_config(page_title="HLA-A coverage", page_icon="🧬", layout="wide")
st.title("HLA-A carrier coverage")
st.write("Compare fixed allele groups across four historical U.S. donor registry categories.")
st.caption("The categories describe this registry study and do not represent countries.")
table = pd.read_csv(PANEL)
f = st.select_slider("Departure from random pairing, F", options=sorted(table.F.unique()))
filtered = table[table.F.eq(f)]
focus = st.selectbox("Allele group", filtered.panel.unique(), index=2)
chart = filtered[filtered.panel.eq(focus)].set_index("donor_group")
st.bar_chart(chart, y="carrier_coverage", color="#258579", y_label="Carrier fraction")
display = filtered[["panel", "donor_group", "allele_mass", "carrier_coverage"]].copy()
st.dataframe(display, hide_index=True, column_config={
    "allele_mass": st.column_config.NumberColumn("Allele mass", format="percent"),
    "carrier_coverage": st.column_config.NumberColumn("Carrier coverage", format="percent"),
})
st.download_button("Download selected scenario", display.to_csv(index=False),
                   file_name="hla-coverage.csv", mime="text/csv")
