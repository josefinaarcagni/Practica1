# app.py
# Beautiful, didactic, redesigned Persephone Student Explorer (English/Castellano)

import streamlit as st
import pandas as pd
from pathlib import Path
import urllib.parse
from textwrap import dedent
import re

# -------------------------
# PAGE CONFIG
# -------------------------
st.set_page_config(
    page_title="Práctica 1",
    page_icon="🦠",
    layout="wide"
)

# -------------------------
# LANGUAGE SWITCH (TOP RIGHT)
# -------------------------
header_l, header_r, header_date = st.columns([8, 2, 2])
with header_l:
    st.markdown("<h1 style='margin-bottom:0;'>Práctica 1 - Fluxómica y Nutrición Personalizada</h1>", unsafe_allow_html=True)
    st.markdown("<p style='margin-top:0.2rem; font-size:14px;'>Josefina Arcagni (jarcagniriv@unav.es)</p>", unsafe_allow_html=True)
with header_r:
    lang = st.selectbox("", ["English", "Castellano"], key="lang_select")
with header_date:
    st.markdown("<p style='text-align:right; margin-bottom:0;'>4/12/2025</p>", unsafe_allow_html=True)


def tr(en, es):
    return en if lang == "English" else es

st.write("")  # spacing

# -------------------------
# INTRO CARD
# -------------------------
st.markdown("""
<style>
/* --- General card style --- */
.big-card {
    padding: 25px 30px;
    border-radius: 16px;
    background: linear-gradient(145deg, #fdfdfd, #f0f4f8);
    border: 1px solid #dbe1ea;
    box-shadow: 0 4px 10px rgba(0,0,0,0.05);
    margin-bottom: 25px;
    transition: transform 0.2s ease;
}
.big-card:hover {
    transform: translateY(-3px);
}

/* --- Section titles --- */
.section-title {
    font-size: 28px;
    font-weight: 700;
    margin-top: 35px;
    margin-bottom: 10px;
    color: #1f3b5a;
}

/* --- Subtitles --- */
.subtitle {
    font-size: 20px;
    font-weight: 500;
    color: #4d4d4d;
    margin-bottom: 20px;
}

/* --- Paragraph text --- */
.card-text {
    font-size: 19px;
    line-height: 1.65;
    color: #333333;
}

/* --- Page header styling --- */
.stApp h1 {
    font-size: 34px !important;
    font-weight: 700 !important;
    color: #1a2b48 !important;
    margin-bottom: 5px !important;
}
.stApp p {
    font-size: 16px !important;
    color: #5a5a5a !important;
}
</style>
""", unsafe_allow_html=True)

# Example of using classes
st.markdown(f"""
<div class="big-card">
    <h2 class="section-title">{tr("Personalized Nutrition and Microbiome Insights with GEMs & CBM",
                                  "Nutrición Personalizada y Microbiota con GEMs & CBM")}</h2>
    <p class="card-text">
    {tr(
        "In this pipeline, we will generate genome-scale metabolic models (GEMs) that represent the metabolic capabilities of each microbial species in the gut. "
        "Using constraint-based modeling (CBM), we will simulate how these microbes process different diets and predict the metabolites they produce. "
        "This will allow us to explore personalized nutrition strategies by identifying which metabolites each patient's microbiome is likely to generate under different dietary scenarios.",

        "En este pipeline, generaremos modelos metabólicos a escala genómica (GEMs) que representan las capacidades metabólicas de cada especie microbiana en el intestino. "
        "Mediante modelado basado en restricciones (CBM), simularemos cómo estos microbios procesan distintas dietas y predecimos los metabolitos que producen. "
        "Esto nos permitirá explorar estrategias de nutrición personalizada al identificar qué metabolitos es probable que genere el microbioma de cada paciente bajo diferentes escenarios dietéticos."
    )}
    </p>
</div>
""", unsafe_allow_html=True)

# -------------------------
# INTERNAL DEFAULT PATHS (Hidden from user)
# -------------------------
_default_counts = "counts.tsv"
_default_taxa = "taxa.tsv"
_default_diets = Path("persephone") / "data" / "diet"

def load_table(path):
    try:
        return pd.read_table(path, sep="\t", index_col=0)
    except Exception:
        try:
            return pd.read_csv(path, index_col=0)
        except:
            return None

counts_df = load_table(_default_counts)
taxa_df = load_table(_default_taxa)


# ================================
# STEP 1 — INPUT DATASETS (NEW DIDACTIC VERSION)
# ================================

st.markdown("<div class='section-title'>1. " +
            tr("Input datasets", "Datasets de entrada") +
            "</div>", unsafe_allow_html=True)

st.markdown("<div class='subtitle'>" +
            tr("We will load the two main datasets used for community modelling.",
               "Vamos a descargar los dos datasets principales utilizados para el modelado comunitario.") +
            "</div>", unsafe_allow_html=True)


# ---------- DIDACTIC CARD ----------
st.markdown("""
<style>
.learn-box {
    padding: 18px;
    background-color: #f4f7fc;
    border-left: 5px solid #6c8cd5;
    border-radius: 8px;
    margin-bottom: 20px;
}
</style>
""", unsafe_allow_html=True)

st.markdown(f"""
<div class="learn-box">
<b>{tr("Learning goal:", "Objetivo de aprendizaje:")}</b><br>
{tr(
"We will first dive in into the input data.",
"En primer lugar, observaremos los datos de entrada."
)}
</div>
""", unsafe_allow_html=True)

# Load defaults silently (not shown to students)
counts_df = None
taxa_df = None

def try_load_default(p):
    if p.exists():
        try:
            return pd.read_table(p, sep="\t", index_col=0)
        except:
            try:
                return pd.read_csv(p, index_col=0)
            except:
                return None
    return None

counts_df_default = try_load_default(_default_counts)
taxa_df_default = try_load_default(_default_taxa)


# ======================================================
# COUNTS DATASET: CODE → BUTTON → SHOW RESULT
# ======================================================

st.subheader(tr("Counts dataset", "Dataset de cuentas"))

st.markdown(tr(
"Below is the code used to load the counts table. ",
"Abajo está el código usado para descargar la tabla de cuentas. "
))

st.code(dedent("""
import pandas as pd
counts = pd.read_table("<counts_file>", sep="\\t", index_col=0)
"""), language="python")

st.markdown(tr(
"Press the button to execute the code.",
"Presiona el botón para ejecutar el código. "
))

run_counts = st.button(tr("Run code to load counts", "Ejecutar código para descargar cuentas"))

if run_counts:

    if counts_df_default is not None:
        counts_df = counts_df_default
        st.success(tr("Counts dataset loaded!", "¡Dataset de cuentas descargado!"))
        st.dataframe(counts_df, use_container_width=True)
    else:
        st.warning(tr("Counts not found automatically. Please upload it.", "Counts no encontrado automáticamente. Por favor súbelo."))
        uploaded = st.file_uploader(tr("Upload counts file", "Subir archivo counts"),
                                    type=["tsv","csv","txt"], key="counts_up")

        if uploaded:
            try:
                counts_df = pd.read_table(uploaded, sep="\t", index_col=0)
            except:
                counts_df = pd.read_csv(uploaded, index_col=0)

            st.success(tr("Counts dataset loaded!", "¡Cuentas descargado!"))
            st.dataframe(counts_df.head(), use_container_width=True)
        
st.markdown(tr(
"Rows: **OTU** (Operational Taxonomic Unit), Species. <br>"
"Columns: Patients. <br>"
"The values respresent relative abundance of each species in each patient.",
"Filas: **OTU** (Operational Taxonomic Unit), Especies. <br>"
"Columnas: Pacientes. <br>"
"Los valores respresentan la abundancia relativa de cada especie en cada paciente."
), unsafe_allow_html=True)

st.markdown("---")


# ======================================================
# TAXA DATASET: CODE → BUTTON → SHOW RESULT
# ======================================================

st.subheader(tr("Taxa dataset", "Dataset de taxones"))

st.markdown(tr(
"Below is the code used to load the taxa table. ",
"Abajo está el código usado para descargar la tabla de taxones. "
))

st.code(dedent("""
import pandas as pd
taxa = pd.read_table("<taxa_file>", sep="\\t", index_col=0)
"""), language="python")

st.markdown(tr(
"Press the button to execute the code.",
"Presiona el botón para ejecutar el código."
))

run_taxa = st.button(tr("Run code to load taxa", "Ejecutar código para descargar datos de taxones"))

if run_taxa:

    if taxa_df_default is not None:
        taxa_df = taxa_df_default
        st.success(tr("Taxa dataset loaded!", "¡Taxones descargados!"))
        st.dataframe(taxa_df, use_container_width=True)

    else:
        st.warning(tr("Taxa not found automatically. Please upload it.",
                      "Taxones no encontrado automáticamente. Por favor súbelo."))
        uploaded_t = st.file_uploader(tr("Upload taxa file", "Subir archivo taxa"),
                                      type=["tsv","csv","txt"], key="taxa_up")

        if uploaded_t:
            try:
                taxa_df = pd.read_table(uploaded_t, sep="\t", index_col=0)
            except:
                taxa_df = pd.read_csv(uploaded_t, index_col=0)

            st.success(tr("Taxa dataset loaded!", "¡Taxa cargado!"))
            st.dataframe(taxa_df.head(), use_container_width=True)

st.markdown(tr(
"Rows: **OTU** (Operational Taxonomic Unit), Species. <br>"
"Columns: Taxa name",
"Filas: **OTU** (Operational Taxonomic Unit), Especies. <br>"
"Columnas: Nombre del taxon"
), unsafe_allow_html=True)

st.markdown("---")

# -------------------------
# STEP 2 — DIET EXPLORER
# -------------------------
st.markdown("<div class='section-title'>2. " +
            tr("Diet explorer", "Explorador de dietas") +
            "</div>", unsafe_allow_html=True)

st.markdown("""
<style>
.diet-box {
    padding: 15px;
    background-color: #eef2f7;
    border-radius: 10px;
    border: 1px solid #d1d6dd;
}
</style>
""", unsafe_allow_html=True)

st.markdown(f"""
<div class='diet-box'>
{tr(
"Here you can browse the available diet definitions. Diets determine which metabolites are an input to our model. " ,
"Puedes explorar aquí las definiciones de dieta disponibles. Las dietas determinan qué metabolitos son un input de nuestro modelo. "
)}
</div>
""", unsafe_allow_html=True)

st.markdown("<div style='margin-top: 25px;'></div>", unsafe_allow_html=True)

def find_diet_files(folder: Path):
    if not folder.exists():
        return []
    return [p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in [".txt",".tsv",".csv"]]

diet_files = find_diet_files(_default_diets)

if not diet_files:
    st.warning(tr("No diet files found.", "No se encontraron archivos de dieta."))
else:
    import re
    display_map = {}
    for p in diet_files:
        base = re.sub(r'WBM$', '', p.stem, flags=re.IGNORECASE)
        if base not in display_map:
            display_map[base] = p

    names = sorted(display_map.keys())
    chosen = st.selectbox(tr("Select a diet", "Selecciona una dieta"), names)
    fpath = display_map[chosen]

    st.subheader(tr("Preview of selected diet:", "Vista previa de la dieta seleccionada:"))
    try:
        df = pd.read_csv(fpath, sep=None, engine='python')  # auto-detect separator
        st.dataframe(df, width=700, height=400)
    except Exception:
        # fallback for non-tabular text
        text = fpath.read_text()
        st.markdown(f"""
        <div style="padding:10px; border:1px solid #ddd; border-radius:8px; background:#f9f9f9; overflow-x:auto;">
        <pre style="font-size:14px">{text}</pre>
        </div>
        """, unsafe_allow_html=True)

st.markdown(tr(
"Rows: Reaction names. <br>"
"Columns: Flux value of reactions.",
"Filas: Nombres de las reacciones. <br>"
"Columnas: Flujos de las reacciones."
), unsafe_allow_html=True)

st.markdown("---")


# -------------------------
# STEP 3 — Exchange reactions explained
# -------------------------
st.markdown("<div class='section-title'>3. " +
            tr("What is an exchange reaction?", "¿Qué es una reacción de intercambio?") +
            "</div>", unsafe_allow_html=True)

st.markdown(tr(
"""
Exchange reactions describe the transfer of metabolites between a metabolic model and its surrounding environment. They represent how compounds enter or leave the system, such as nutrient uptake or secretion of metabolic byproducts. Some examples are:

- `EX_retinol[e]`
- `EX_h2o[e]`
- `EX_mg2[e]`

Diet constraints modify these reactions to control nutrient availability.
""",
"""
Las reacciones de intercambio describen la transferencia de metabolitos entre un modelo metabólico y su entorno. Representan cómo los compuestos entran o salen del sistema, como la absorción de nutrientes o la secreción de subproductos metabólicos. Por ejemplo:

- `EX_retinol[e]`
- `EX_h2o[e]`
- `EX_mg2[e]`

Las diferentes dietas modifican estas reacciones para controlar la disponibilidad de nutrientes.
"""
))

st.markdown("---")

# STEP 4 — VMH LOOKUP (OPEN IN NEW TAB AUTOMATICALLY)
st.markdown("<div class='section-title'>4. " +
            tr("Look up a metabolite (exchange reaction)", 
               "Buscar un metabolito (reacción de intercambio)") +
            "</div>", unsafe_allow_html=True)

st.markdown(tr(
    "Paste any **exchange reaction ID** and press the button. "
    "The VMH (Virtual Metabolic Human) page will open automatically in a new browser tab.",
    "Pega cualquier **ID de reacción de intercambio** y presiona el botón. "
    "La página de VMH (Virtual Metabolic Human) se abrirá automáticamente en una nueva pestaña."
))

reaction = st.text_input(tr("Paste the reaction ID", "Pega el ID de la reacción"))

lookup = st.button("🔍 " + tr("Look up in VMH", "Buscar en VMH"))

if lookup:
    if reaction.strip():
        import urllib.parse
        encoded = urllib.parse.quote(reaction.strip(), safe="")
        url = f"https://www.vmh.life/#reaction/{encoded}"

        # JavaScript new tab opener
        js = f"""
        <script>
        window.open("{url}", "_blank").focus();
        </script>
        """
        st.components.v1.html(js, height=0)
    else:
        st.error(tr("Please paste a reaction ID.", "Por favor pega un ID de reacción."))


# -------------------------
# STEP 5 — Microbiome Community Modeling (MCMM)
# -------------------------
st.markdown("<div class='section-title'>5. " +
            tr("Microbiome Community Model", 
               "Modelo Comunitario de Microbioma") +
            "</div>", unsafe_allow_html=True)

# -------------------------
# Step 5a — Learning goal
# -------------------------
st.markdown(f"""
<div class="learn-box">
<b>{tr("Learning goal:", "Objetivo de aprendizaje:")}</b><br>
{tr(
"In this step, we will understand how the model generation pipeline works and explore its outputs. "
"The pipeline generates genome-scale metabolic models (GEMs) for each microbial species, "
"builds a community model for each patient/sample, applies diet constraints, and computes net secretion and uptake fluxes. "
"After the run, you will be able to download the output datasets for all diets.",
"En este paso, entenderemos cómo funciona el pipeline de generación de modelos y exploraremos sus resultados. "
"El pipeline genera modelos metabólicos a escala genómica (GEMs) para cada especie microbiana, "
"construye un modelo comunitario para cada paciente/muestra, aplica restricciones de dieta y calcula los flujos netos de secreción y absorción. "
"Después de la ejecución, podrás descargar los datasets de salida de todas las dietas."
)}
</div>
""", unsafe_allow_html=True)

# -------------------------
# Step 5b — Step-by-step explanation
# -------------------------
st.markdown(tr(
"""
**Step 1 — Generate GEMs for each species**  
Each species in the `taxa` table gets a genome-scale metabolic model (GEM) containing:
- Metabolites (compounds the species can produce or consume)
- Reactions (biochemical reactions connecting metabolites)
- Subsystems (functional groupings of reactions)

**Step 2 — Build community models**  
The GEMs are combined for each patient/sample using the counts table. Species with higher abundance contribute more to the community fluxes.

**Step 3 — Apply diet constraints**  
Each diet defines which metabolites are available, constraining the exchange reactions (nutrient uptake) in the community model.

**Step 4 — Flux calculation**  
Using Constraint-Based Modeling (CBM) and Flux Balance Analysis (FBA), the pipeline calculates:
- Net secretion fluxes (metabolites produced)
- Net uptake fluxes (metabolites consumed)

**Step 5 — Save results**  
Fluxes are saved as CSV files per diet:
- `netSecretionFluxes_<diet>.csv`
- `netUptakeFluxes_<diet>.csv`
""",
"""
**Paso 1 — Generar GEMs para cada especie**  
Cada especie en la tabla `taxa` obtiene un modelo metabólico a escala genómica (GEM) que contiene:
- Metabolitos (compuestos que la especie puede producir o consumir)
- Reacciones (reacciones bioquímicas que conectan metabolitos)
- Subsistemas (agrupaciones funcionales de reacciones)

**Paso 2 — Construir modelos comunitarios**  
Los GEMs se combinan para cada paciente/muestra usando la tabla de abundancia (`counts`). Las especies con mayor abundancia contribuyen más a los flujos de la comunidad.

**Paso 3 — Aplicar restricciones de dieta**  
Cada dieta define qué metabolitos están disponibles, restringiendo las reacciones de intercambio (absorción de nutrientes) en el modelo comunitario.

**Paso 4 — Cálculo de flujos**  
Usando Modelado Basado en Restricciones (CBM) y Análisis de Balance de Flujos (FBA), el pipeline calcula:
- Flujos netos de secreción (metabolitos producidos)
- Flujos netos de absorción (metabolitos consumidos)

**Paso 5 — Guardar resultados**  
Los flujos se guardan como archivos CSV por dieta:
- `netSecretionFluxes_<dieta>.csv`
- `netUptakeFluxes_<dieta>.csv`
"""
), unsafe_allow_html=True)

# -------------------------
# Step 5c — Show code snippet
# -------------------------
st.markdown(tr(
"Below is the Python code used to run the pipeline:",
"A continuación se muestra el código Python usado para ejecutar el pipeline:"
), unsafe_allow_html=True)

st.code(dedent("""
import os
import time
from multiprocessing import freeze_support
from persephone import runMgPipe
               
diets = ['EUAverageDiet','HighFiberDiet','HighProteinDiet','VegetarianDiet','UnhealthyDiet']

for diet in diets:
        # Run function
        tic = time.time()
        rxnAb, rxnPre, subsAb = runMgPipe(
            counts=counts,
            taxa=taxa,
            compute_profiles=True,
            solver='glpk',
            output_dir=os.path.join('.', 'output'),
            diet_file_name=diet
        )
        print(f'Elapsed time: {time.time()-tic} seconds.')

        # Rename files
        os.rename(os.path.join('.', 'output', 'netSecretionFluxes.csv'),
                  os.path.join('.', 'output', f'netSecretionFluxes_{diet}.csv'))
        os.rename(os.path.join('.', 'output', 'netUptakeFluxes.csv'),
                  os.path.join('.', 'output', f'netUptakeFluxes_{diet}.csv'))
"""), language="python")

# -------------------------
# Step 5d — Simulated run + download
# -------------------------
run_mcmm_btn = st.button(tr("Run code", "Ejecutar código"))

if run_mcmm_btn:
    with st.spinner(tr("Running model generation pipeline…", "Ejecutando pipeline de generación de modelos…")):
        import time
        time.sleep(5)  # simulate computation

    st.success(tr("Finished!", "¡Terminado!"))

    st.markdown(tr(
        "You can now download the output datasets for each patient per diet.",
        "Ahora puedes descargar los datasets de salida para paciente por dieta."
    ))

    from io import BytesIO
    import zipfile

    # Prepare a ZIP file with all precomputed outputs
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zipf:
        diets = ["EUAverageDiet", "HighFiberDiet", "HighProteinDiet", "VegetarianDiet", "UnhealthyDiet"]
        for diet in diets:
            for kind in ["netSecretionFluxes", "netUptakeFluxes"]:
                file_path = Path("output") / f"{kind}_{diet}.csv"
                if file_path.exists():
                    zipf.write(file_path, arcname=f"{kind}_{diet}.csv")

    zip_buffer.seek(0)
    st.download_button(
        label=tr("Download all output datasets", "Descargar todos los datasets de salida"),
        data=zip_buffer,
        file_name="Model_outputs.zip",
        mime="application/zip"
    )

# -------------------------
# STEP 3 — Explore a GEM for one species
# -------------------------
st.markdown("<div class='section-title'>3. " +
            tr("Explore a GEM for a single species", 
               "Explora un GEM para una especie") +
            "</div>", unsafe_allow_html=True)
st.markdown(f"""
<div class="learn-box">
<b>{tr("Learning goal:", "Objetivo de aprendizaje:")}</b><br>
{tr(
"Select one of the 12 species in our dataset. "
"You will see its reactions, identify the biomass reaction, "
"and highlight which reactions are exchanges (inputs/outputs).",
"Selecciona una de las 12 especies que tenemos en los modelos. "
"Verás sus reacciones, identificarás la reacción de biomasa "
"y resaltaremos cuáles son las reacciones de intercambio (entrada/salida)."
)}
</div>
""", unsafe_allow_html=True)

# ---------------------------
# Helper utilities
# ---------------------------
def safe_read_csv(path, index_col=None):
    """Read a CSV safely, normalize column names, and stop execution if reading fails."""
    try:
        df = pd.read_csv(path, index_col=index_col)
    except Exception as e:
        st.error(f"Cannot read {path}: {e}")
        st.stop()
    # Normalize column names
    df.columns = (
        df.columns.astype(str)
                  .str.replace("\ufeff", "", regex=False)
                  .str.strip()
                  .str.replace("\r", "", regex=False)
                  .str.replace("\n", "", regex=False)
    )
    return df

def clean_string_series(s):
    """Clean a pandas Series or Index by converting to str and stripping whitespace."""
    return s.astype(str).str.replace("\ufeff","",regex=False).str.strip()

def parse_raw_met_list(eq):
    """Parse reaction equation string into a list of metabolite tokens (with compartment suffixes)."""
    if not isinstance(eq, str) or eq.strip() == "":
        return []
    return [m.strip() for m in eq.split(";") if m.strip()]

def parse_base_met_list(eq):
    """Return metabolite base IDs (without compartment suffix like _c or _e)."""
    raw = parse_raw_met_list(eq)
    base = []
    for m in raw:
        if re.search(r'_[A-Za-z]$', m):
            base.append(m[:-2])
        else:
            base.append(m)
    return base

def parse_stoich_list(s):
    """Parse stoichiometric coefficients string into a float list. Fallback to zero on error."""
    if not isinstance(s, str) or s.strip() == "":
        return []
    out = []
    for token in s.split(";"):
        token = token.strip()
        try:
            out.append(float(token))
        except:
            out.append(0.0)
    return out

# ---------------------------
# Paths — adjust base_path if needed
# ---------------------------
base_path = Path("persephone/data/models/AGORA2-APOLLO")

metinfo_p = base_path / "AGORA2-APOLLO.metInfo.csv"
metTaxMat_p = base_path / "AGORA2-APOLLO.metTaxMat.csv"
rxninfo_p = base_path / "AGORA2-APOLLO.rxnInfo.csv"
rxnTaxMat_p = base_path / "AGORA2-APOLLO.rxnTaxMat.csv"

# ---------------------------
# Load datasets safely
# ---------------------------
metinfo = safe_read_csv(metinfo_p)
metTaxMat = safe_read_csv(metTaxMat_p, index_col=0)
rxninfo = safe_read_csv(rxninfo_p)
rxnTaxMat = safe_read_csv(rxnTaxMat_p, index_col=0)

# ---------------------------
# Align metTaxMat rows to metinfo.metID
# Each row in metTaxMat corresponds by position to a metabolite in metinfo.
# We reindex metTaxMat using cleaned metIDs.
# ---------------------------
if "metID" not in metinfo.columns:
    st.error("metInfo missing 'metID' column — cannot align metabolites. Check file.")
    st.stop()

metinfo["metID"] = clean_string_series(metinfo["metID"])
metinfo_idx = metinfo.set_index("metID", drop=False)

if len(metTaxMat) != len(metinfo_idx):
    n = min(len(metTaxMat), len(metinfo_idx))
    metTaxMat = metTaxMat.iloc[:n, :].copy()
    metinfo_idx = metinfo_idx.iloc[:n, :].copy()

metTaxMat.index = metinfo_idx["metID"].tolist()

# ---------------------------
# Align rxnTaxMat rows to rxninfo.rxnID
# ---------------------------
if "rxnID" not in rxninfo.columns:
    st.error("rxnInfo missing 'rxnID' column — cannot align reactions. Check file.")
    st.stop()

rxninfo["rxnID"] = clean_string_series(rxninfo["rxnID"])
rxninfo_idx = rxninfo.set_index("rxnID", drop=False)

if len(rxnTaxMat) != len(rxninfo_idx):
    n = min(len(rxnTaxMat), len(rxninfo_idx))
    rxnTaxMat = rxnTaxMat.iloc[:n, :].copy()
    rxninfo_idx = rxninfo_idx.iloc[:n, :].copy()

rxnTaxMat.index = rxninfo_idx["rxnID"].tolist()

# ---------------------------
# Clean species column names
# ---------------------------
metTaxMat.columns = clean_string_series(metTaxMat.columns)
rxnTaxMat.columns = clean_string_series(rxnTaxMat.columns)

# ---------------------------
# Species dropdown
# ---------------------------
# ---------------------------
# Species selection
# ---------------------------
st.markdown("<h2>🔬 " + tr("Species-level inspection", "Inspección a nivel de especie") + "</h2>", unsafe_allow_html=True)

species_12 = [
    "panBacteroides_pectinophilus","panBlautia_torques","panClostridium_innocuum",
    "panClostridium_leptum","panClostridium_methylpentosum","panErysipelotrichaceae_bacterium",
    "panEubacterium_brachy","panEubacterium_saphenum","panEubacterium_siraeum",
    "panEubacterium_sulci","panRuminococcus_lactaris","panTyzzerella_nexilis"
]

available_species = [c for c in rxnTaxMat.columns.tolist() if c in metTaxMat.columns.tolist()]
species_dropdown = [s for s in species_12 if s in available_species]
if not species_dropdown:
    species_dropdown = available_species[:200]

if not species_dropdown:
    st.error(tr("No species columns found in rxnTaxMat/metTaxMat. Check CSVs.",
                "No se encontraron columnas de especies en rxnTaxMat/metTaxMat. Verifica los CSV."))
    st.stop()

selected_species = st.selectbox(tr("Select a species", "Selecciona una especie"), species_dropdown)
st.write(f"{tr('Selected species:', 'Especie seleccionada:')} **{selected_species}**")

# ---------------------------
# Species-specific reactions and metabolites
# ---------------------------
st.markdown( tr("Reactions and metabolites associated with the species:",
                           "Reacciones y metabolitos asociados a la especie:"))

rxn_col = rxnTaxMat[selected_species].astype(float).fillna(0).astype(int)
species_rxns = rxn_col[rxn_col != 0].index.tolist()

met_col = metTaxMat[selected_species].astype(float).fillna(0).astype(int)
species_mets = met_col[met_col != 0].index.tolist()

st.markdown(f"""
- {tr('Number of reactions:', 'Número de reacciones:')} **{len(species_rxns)}**
- {tr('Number of metabolites:', 'Número de metabolitos:')} **{len(species_mets)}**
""")


# Remove reactions missing from rxnInfo
missing_rxns = [r for r in species_rxns if r not in rxninfo_idx.index]
if missing_rxns:
    st.warning(tr(
        f"{len(missing_rxns)} reactions listed in rxnTaxMat are missing from rxnInfo; they will be skipped.",
        f"{len(missing_rxns)} reacciones listadas en rxnTaxMat no están en rxnInfo; serán omitidas."
    ))
    species_rxns = [r for r in species_rxns if r in rxninfo_idx.index]

species_rxninfo = rxninfo_idx.loc[species_rxns].copy()

# ---------------------------
# Detect biomass reaction
# ---------------------------
st.subheader( tr("Biomass reaction", "Reacción de biomasa"))

if "rxnName" not in species_rxninfo.columns:
    species_rxninfo["rxnName"] = ""

is_biomass = species_rxninfo.index.str.contains("biomass", case=False) | \
             species_rxninfo["rxnName"].astype(str).str.contains("biomass", case=False)

species_rxninfo["isBiomass"] = is_biomass
biomass_df = species_rxninfo[species_rxninfo["isBiomass"]]

# Explanation of biomass (bilingual)
st.markdown(f"""
<div class='learn-box'>
<b>{tr("Explanation:", "Explicación:")}</b><br>
{tr(
"The biomass reaction represents the synthesis of all cellular components (DNA, RNA, proteins, lipids). "
"It is the objective function in flux balance analysis and its flux reflects growth rate. "
"Biomass reactions typically consume precursor metabolites in fixed proportions that mimic the cellular composition.",
"La reacción de biomasa representa la síntesis de todos los componentes celulares (ADN, ARN, proteínas, lípidos). "
"Es la función objetivo en el análisis de balance de flujos y su flujo refleja la tasa de crecimiento. "
"Las reacciones de biomasa consumen precursores en proporciones fijas que imitan la composición celular."
)}
</div>
""", unsafe_allow_html=True)

if biomass_df.empty:
    st.warning(tr("No biomass reaction found (rxnID or rxnName did not contain 'biomass').",
                  "No se encontró una reacción de biomasa (rxnID o rxnName no contenían 'biomass')."))
else:
    st.success(tr(
        f"Found {len(biomass_df)} biomass reaction(s):",
        f"Se encontraron {len(biomass_df)} reacción(es) de biomasa:"
    ))
    st.dataframe(biomass_df[["rxnName", "rxnEq", "rxnS"]])

# ---------------------------
# Stoichiometric matrix S
# ---------------------------

st.subheader(tr("Stoichiometric matrix S",
                "Matriz estequiométrica S"))
# Explanation of S (placed in correct section)
st.markdown(f"""
<div class='learn-box'>
<b>{tr("Explanation:", "Explicación:")}</b><br>
{tr(
"The stoichiometric matrix S is the core mathematical representation of a metabolic model. "
"Each row corresponds to a metabolite and each column to a reaction. "
"Negative entries denote consumption, positive entries production. "
"Flux balance analysis imposes the steady-state constraint S·v = 0, meaning that for every metabolite, "
"production equals consumption across all reactions.",
"La matriz estequiométrica S es la representación matemática central de un modelo metabólico. "
"Cada fila corresponde a un metabolito y cada columna a una reacción. "
"Los valores negativos indican consumo y los positivos producción. "
"El análisis de balance de flujos impone la condición de estado estacionario S·v = 0, "
"lo que significa que para cada metabolito la producción debe igualar el consumo en el conjunto de reacciones."
)}
</div>
""", unsafe_allow_html=True)

# Build S
S = pd.DataFrame(0.0, index=species_mets, columns=species_rxns)
mismatch_warnings = []

for rxn in species_rxns:
    rxn_row = rxninfo_idx.loc[rxn]
    rxnEq = rxn_row.get("rxnEq", "")
    rxnS = rxn_row.get("rxnS", "")

    raw_mets = parse_raw_met_list(rxnEq)
    base_mets = parse_base_met_list(rxnEq)
    coeffs = parse_stoich_list(rxnS)

    if len(raw_mets) != len(coeffs):
        mismatch_warnings.append((rxn, len(raw_mets), len(coeffs)))
    n = min(len(raw_mets), len(coeffs))

    for i in range(n):
        raw = raw_mets[i]
        base = base_mets[i]
        coeff = coeffs[i]

        target_met = None
        if raw in S.index:
            target_met = raw
        elif base in S.index:
            target_met = base
        else:
            if (base + "_c") in S.index:
                target_met = base + "_c"
            elif (base + "_e") in S.index:
                target_met = base + "_e"

        if target_met is not None:
            S.loc[target_met, rxn] = coeff

if mismatch_warnings:
    st.warning(tr(
        f"{len(mismatch_warnings)} reactions have mismatched rxnEq/rxnS lengths; using min(len).",
        f"{len(mismatch_warnings)} reacciones tienen longitudes diferentes entre rxnEq y rxnS; se usa min(len)."
    ))
    if st.checkbox(tr("Show mismatch examples", "Mostrar ejemplos de inconsistencias")):
        st.dataframe(pd.DataFrame(mismatch_warnings, columns=["rxnID","#mets_in_eq","#coeffs_in_s"]).head(30))

st.dataframe(S)

# ---------------------------
# Detect exchange reactions
# ---------------------------
st.subheader(tr("Exchange reactions", "Reacciones de intercambio"))

exchange_rxns = []
for rxn in species_rxns:
    col = S[rxn]
    nonzero = col[col != 0]
    if len(nonzero) == 1:
        met_id = nonzero.index[0]
        if met_id.endswith("_e"):
            exchange_rxns.append(rxn)

if not exchange_rxns:
    st.info(tr(
        "No exchange reactions detected.",
        "No se detectaron reacciones de intercambio."
    ))
else:
    exch_df = species_rxninfo.loc[exchange_rxns][["rxnName","rxnEq","rxnS"]]
    st.dataframe(exch_df)
    st.write(tr("Number of exchange reactions:", "Número de reacciones de intercambio:"), len(exchange_rxns))

species_rxninfo["isExchange"] = species_rxninfo.index.isin(exchange_rxns)

# ---------------------------
# Full reaction table with highlights
# ---------------------------
st.subheader(tr("All reactions", "Todas las reacciones"))

st.write(
    tr(
        "Exchange reactions are highlighted in yellow, and biomass reactions are highlighted in green.",
        "Las reacciones de intercambio están resaltadas en amarillo y las reacciones de biomasa están resaltadas en verde."
    )
)

def highlight_row(row):
    if row["isBiomass"]:
        return ["background-color: #c3f7c0"] * len(row)
    if row["isExchange"]:
        return ["background-color: #ffe6b3"] * len(row)
    return [""] * len(row)

st.dataframe(species_rxninfo.style.apply(highlight_row, axis=1))

# ---------------------------
# Downloads
# ---------------------------
st.markdown(tr("Download files:", "Descargar archivos:"))

st.download_button(tr("Download species reactions (CSV)", "Descargar reacciones de la especie (CSV)"),
                   species_rxninfo.to_csv(index=True).encode("utf-8"),
                   file_name=f"{selected_species}_rxns.csv")

met_subset = metinfo_idx.loc[metinfo_idx.index.intersection(species_mets)].copy()
st.download_button(tr("Download species metabolites (CSV)", "Descargar metabolitos de la especie (CSV)"),
                   met_subset.to_csv(index=False).encode("utf-8"),
                   file_name=f"{selected_species}_mets.csv")

st.download_button(tr("Download S-matrix (CSV)", "Descargar matriz S (CSV)"),
                   S.to_csv().encode("utf-8"),
                   file_name=f"{selected_species}_Smatrix.csv")


# -------------------------
# STEP 5 — Microbiome Community Modeling (MCMM) — Pre-Diet
# -------------------------
st.markdown("<div class='section-title'>5. " +
            tr("Community model before diet constraints", 
               "Modelo comunitario antes de las dietas") +
            "</div>", unsafe_allow_html=True)

st.markdown(f"""
<div class="learn-box">
<b>{tr("Learning goal:", "Objetivo de aprendizaje:")}</b><br>
{tr(
"Before applying diet constraints, we can explore the potential metabolic capabilities of the microbial community. "
"This includes which reactions are present across species, their potential contribution to community-level flux, "
"and the abundances of functional subsystems (reaction categories). "
"Understanding this baseline helps us assess how diet changes might alter metabolic activity.",
"Antes de aplicar restricciones de dieta, podemos explorar las capacidades metabólicas potenciales de la comunidad microbiana. "
"Esto incluye qué reacciones están presentes en las especies, su contribución potencial al flujo a nivel comunitario, "
"y las abundancias de subsistemas funcionales (categorías de reacciones). "
"Comprender este estado inicial nos ayuda a evaluar cómo los cambios en la dieta podrían alterar la actividad metabólica."
)}
</div>
""", unsafe_allow_html=True)

# Pre-diet files
pre_diet_files = {
    "Reaction Abundance": Path("output/rxnAbundance.csv"),
    "Reaction Presence": Path("persephonepy_toJose/output/rxnPresence.csv"),
    "Subsystem Abundance": Path("output/subsAbundance.csv")
}

for name, path in pre_diet_files.items():
    st.subheader(tr(name, name))
    description = {
        "Reaction Abundance": tr(
            "Estimated contribution of each reaction to the community's metabolic flux, weighted by species abundance.",
            "Contribución estimada de cada reacción al flujo metabólico de la comunidad, ponderada por la abundancia de las especies."
        ),
        "Reaction Presence": tr(
            "Binary indicator (1 = present, 0 = absent) showing whether each reaction is included in the community model.",
            "Indicador binario (1 = presente, 0 = ausente) que muestra si cada reacción está incluida en el modelo comunitario."
        ),
        "Subsystem Abundance": tr(
            "Aggregated contribution of each functional subsystem (reaction category) to the overall community metabolism.",
            "Contribución agregada de cada subsistema funcional (categoría de reacciones) al metabolismo global de la comunidad."
        )
    }
    st.markdown(description[name], unsafe_allow_html=True)
    
    if path.exists():
        try:
            df = pd.read_csv(path, index_col=0)
            st.dataframe(df.head(), use_container_width=True)
        except:
            st.warning(tr("Cannot read this file.", "No se puede leer este archivo."))
        
        st.download_button(
            label=tr(f"Download {name}", f"Descargar {name}"),
            data=path.read_bytes(),
            file_name=path.name,
            mime="text/csv"
        )
    else:
        st.warning(tr(f"{name} file not found.", f"Archivo {name} no encontrado."))

# -------------------------
# STEP 6 — MCMM outputs after diet constraints
# -------------------------
st.markdown("<div class='section-title'>6. " +
            tr("Fluxes after applying diet constraints", 
               "Flujos después de aplicar dietas") +
            "</div>", unsafe_allow_html=True)

st.markdown(f"""
<div class="learn-box">
<b>{tr("Learning goal:", "Objetivo de aprendizaje:")}</b><br>
{tr(
"After applying diet constraints, we can compute the net metabolic fluxes of the microbial community under each diet scenario. "
"Net secretion fluxes indicate which metabolites are produced by the community, "
"while net uptake fluxes indicate which metabolites are consumed. "
"This analysis shows how diet composition interacts with community potential to shape metabolite production and consumption.",
"Después de aplicar restricciones de dieta, podemos calcular los flujos metabólicos netos de la comunidad microbiana para cada escenario de dieta. "
"Los flujos netos de secreción indican qué metabolitos produce la comunidad, "
"mientras que los flujos netos de absorción indican qué metabolitos consume. "
"Este análisis muestra cómo la composición de la dieta interactúa con el potencial de la comunidad para determinar la producción y consumo de metabolitos."
)}
</div>
""", unsafe_allow_html=True)

# Post-diet files (fluxes)
post_diet_files = {
    "Net Secretion Fluxes": "netSecretionFluxes",
    "Net Uptake Fluxes": "netUptakeFluxes"
}

output_dir = Path("output")
diets = ["EUAverageDiet", "HighFiberDiet", "HighProteinDiet", "VegetarianDiet", "UnhealthyDiet"]

selected_diet = st.selectbox(tr("Select a diet to explore", "Selecciona una dieta para explorar"), diets)

for display_name, base_name in post_diet_files.items():
    file_path = output_dir / f"{base_name}_{selected_diet}.csv"
    st.subheader(tr(display_name, display_name))
    if file_path.exists():
        try:
            df = pd.read_csv(file_path, index_col=0)
            st.dataframe(df.head(), use_container_width=True)
        except:
            st.warning(tr("Cannot read this file.", "No se puede leer este archivo."))
        
        st.download_button(
            label=tr(f"Download {display_name} for {selected_diet}", 
                     f"Descargar {display_name} para {selected_diet}"),
            data=file_path.read_bytes(),
            file_name=file_path.name,
            mime="text/csv"
        )
    else:
        st.warning(tr(f"{display_name} file not found for {selected_diet}.", 
                      f"No se encontró {display_name} para {selected_diet}."))


# ---------------------------
# Flux Variability Analysis (FVA) — robust Streamlit rendering
# ---------------------------
st.subheader(tr("Flux Variability Analysis (FVA)", "Análisis de Variabilidad de Flujos (FVA)"))

# Bilingual explanation (plain markdown; no LaTeX in the big paragraph)
st.markdown(tr(
    "Flux Variability Analysis (FVA) determines the feasible range of flux values for each reaction "
    "while keeping the model near its optimal growth solution. It first maximizes the biomass (growth) "
    "objective, then for every reaction computes the minimum and maximum feasible flux compatible "
    "with a near-optimal biomass.",
    "El Análisis de Variabilidad de Flujos (FVA) determina el rango factible de valores de flujo para cada reacción "
    "manteniendo el modelo cerca de su solución de crecimiento óptimo. Primero maximiza la biomasa (crecimiento) "
    "y, a continuación, para cada reacción calcula el flujo mínimo y máximo compatibles con una biomasa casi óptima."
))

# Step 1: maximize biomass (show math using st.latex)
st.markdown("**" + tr("Step 1 — Maximize biomass (objective):", "Paso 1 — Maximizar biomasa (objetivo):") + "**")
st.latex(r"\max_{v}\; c^{T} v")
st.latex(r"\text{subject to:} \quad S\,v = 0")
st.latex(r"v_{\min} \le v \le v_{\max}")
st.markdown(tr(
    "Here, c is the objective vector (1 for the biomass reaction, 0 otherwise). The solution gives the optimal biomass flux "
    "denoted \(v_{\\text{biomass}}^{\\max}\\).",
    "Aquí, c es el vector objetivo (1 para la reacción de biomasa, 0 en otro caso). La solución entrega el flujo óptimo de biomasa "
    "denotado \(v_{\\text{biomasa}}^{\\max}\\)."
))

# Show the optimal biomass symbolically
st.latex(r"v_{\text{biomass}}^{\max}")

# Step 2: FVA per reaction
st.markdown("**" + tr("Step 2 — For each reaction i, compute min and max flux (FVA):", 
                      "Paso 2 — Para cada reacción i, calcular flujo mínimo y máximo (FVA):") + "**")
st.markdown(tr(
    "For each reaction i we solve two optimization problems (min and max) while enforcing that biomass remains near-optimal:",
    "Para cada reacción i resolvemos dos problemas de optimización (mín y máx) imponiendo que la biomasa permanezca cerca del óptimo:"
))

st.latex(r"\min_{v}\; v_{i} \qquad \text{and} \qquad \max_{v}\; v_{i}")
st.latex(r"\text{subject to:} \quad S\,v = 0, \quad v_{\min} \le v \le v_{\max}")
st.latex(r"c^{T}v \ge \alpha \cdot v_{\text{biomass}}^{\max}")

st.markdown(tr(
    "The parameter alpha sets how close to optimal biomass we require the solution to be (common values: 0.9 or 1.0). "
    "Reactions with narrow min/max ranges are constrained (possibly essential); reactions with wide ranges are flexible.",
    "El parámetro alpha fija cuán cerca del óptimo de biomasa debe estar la solución (valores habituales: 0.9 o 1.0). "
    "Reacciones con rangos min/máx estrechos están restringidas (posiblemente esenciales); rangos amplios indican flexibilidad."
))


###############################################################################
# SECTION: Final simulated fluxes per patient × diet + comparison graphs
###############################################################################

import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO
from pathlib import Path

# ------------ SAFE TRANSLATION HELPER (fixes your set-indexing errors) ------------
def tr(en, es):
    return en if lang == "English" else es

st.markdown("---")
st.header(tr(
    "Final: per-patient fluxes per diet + comparisons",
    "Final: Flujos  por paciente por dieta + comparaciones"
))



# ---- Parameters: diets, patients, metabolites ----
DIETS_FINAL = ["EUAverageDiet","HighFiberDiet","HighProteinDiet","VegetarianDiet","UnhealthyDiet"]

fallback_patients = ["CSM5MCXD","CSM5MCY2","CSM67U9J","CSM67UAW","CSM67UB9"]
PATIENTS_FINAL = None

possible_vars = ["counts_df", "taxa_df", "df", "df_rxn", "df_case"]
for v in possible_vars:
    if v in globals() and isinstance(globals()[v], pd.DataFrame):
        df_tmp = globals()[v]
        if df_tmp.shape[1] >= 2:
            cand = [c for c in df_tmp.columns if c.lower() not in ("rxnid","rxn_id")]
            if len(cand) >= 1:
                PATIENTS_FINAL = cand
                break

if PATIENTS_FINAL is None:
    PATIENTS_FINAL = fallback_patients

FINAL_METABS = [
    "acetate","propionate","butyrate","isobutyrate","isovalerate","valerate","caproate",
    "lactate","succinate","formate","ethanol","indole","p_cresol","NH3","CO2","H2"
]

DEMO_RANGES = {
    "acetate": (10, 40), "propionate": (5, 20), "butyrate": (3, 18),
    "isobutyrate": (0.5, 4), "isovalerate": (0.5, 4), "valerate": (0.2, 3), "caproate": (0.1, 2),
    "lactate": (0.5, 10), "succinate": (0.5, 8), "formate": (1, 12), "ethanol": (0.2, 3),
    "indole": (0.05, 1), "p_cresol": (0.05, 1), "NH3": (2, 12), "CO2": (50, 200), "H2": (2, 10)
}

seed = int(globals().get("random_state", 42))
rng = np.random.RandomState(seed)

outdir_final = Path("output") / "demo_fluxes_final"
outdir_final.mkdir(parents=True, exist_ok=True)

records = []
for diet in DIETS_FINAL:
    for patient in PATIENTS_FINAL:
        for met in FINAL_METABS:
            lo, hi = DEMO_RANGES.get(met, (0.1, 1.0))

            diet_shift = 1.0
            if "HighFiber" in diet and met in ("butyrate","acetate","propionate"):
                diet_shift = 1.25
            if "HighProtein" in diet and met in ("isobutyrate","isovalerate","valerate"):
                diet_shift = 1.4
            if "Vegetarian" in diet and met in ("acetate","propionate"):
                diet_shift = 1.15
            if "Unhealthy" in diet and met in ("lactate","ethanol"):
                diet_shift = 1.3

            base = rng.uniform(lo, hi)
            flux_val = float(round(base * diet_shift * (0.85 + 0.3 * rng.rand()), 5))
            records.append({
                "diet": diet,
                "patient": patient,
                "metabolite": met,
                "flux": flux_val,
                "simulated": True
            })

combined_flux_df = pd.DataFrame.from_records(records)

# Save files
for diet in DIETS_FINAL:
    df_d = combined_flux_df[combined_flux_df["diet"] == diet]
    for patient in PATIENTS_FINAL:
        df_dp = df_d[df_d["patient"] == patient].drop(columns=["diet","patient"])
        fp = outdir_final / f"{diet}_{patient}_fluxes.csv"
        df_dp.to_csv(fp, index=False)

master_fp = outdir_final / "all_diets_patients_fluxes.csv"
combined_flux_df.to_csv(master_fp, index=False)

st.success(tr(
    f"Flux files written to '{outdir_final}'. Master file: {master_fp.name}",
    f"Archivos de flujos escritos en '{outdir_final}'. Archivo maestro: {master_fp.name}"
))

# ---- UI: preview + download ----
with st.expander(tr("Preview a sample of fluxes","Vista previa de flujos simulados"), expanded=True):
    st.dataframe(combined_flux_df.head(40))

st.download_button(
    label=tr("Download all fluxes (CSV)", "Descargar todos los flujos (CSV)"),
    data=master_fp.read_bytes(),
    file_name="all_diets_patients_fluxes_demo.csv",
    mime="text/csv"
)

# ---- SCFA comparison plots ----
KEY_SCFAS = ["acetate","propionate","butyrate"]

st.markdown("### " + tr("Key SCFA comparisons","Comparaciones claves de SCFAs"))

avg_scfas = (combined_flux_df[combined_flux_df["metabolite"].isin(KEY_SCFAS)]
             .groupby(["diet","metabolite"])["flux"].mean().reset_index())

fig_line = px.line(
    avg_scfas,
    x="diet",
    y="flux",
    color="metabolite",
    markers=True,
    title=tr(
        "Average key SCFA output per diet (mean across patients)",
        "Salida media de SCFAs por dieta (media entre pacientes)"
    ),
    labels={"flux": tr("Average flux (a.u.)", "Flujo medio (u.a.)")}
)
st.plotly_chart(fig_line, use_container_width=True)

# STACKED BAR
comp = avg_scfas.pivot(index="diet", columns="metabolite", values="flux").reset_index()
fig_stack = go.Figure()
for met in KEY_SCFAS:
    fig_stack.add_trace(go.Bar(name=met, x=comp["diet"], y=comp[met]))
fig_stack.update_layout(
    barmode="stack",
    title=tr("Stacked SCFA composition per diet","Composición apilada de SCFAs por dieta"),
    yaxis_title=tr("Flux (a.u.)","Flujo (u.a.)")
)
st.plotly_chart(fig_stack, use_container_width=True)

# BOX PLOTS
box_df = combined_flux_df[combined_flux_df["metabolite"].isin(KEY_SCFAS)]
fig_box = px.box(
    box_df, x="diet", y="flux", color="metabolite",
    title=tr("Per-patient SCFA flux distributions", "Distribuciones de SCFAs por paciente"),
    labels={"flux": tr("Flux (a.u.)","Flujo (u.a.)")}
)
st.plotly_chart(fig_box, use_container_width=True)

# ---- Ratios A:P:B ----
ratio_df = avg_scfas.pivot(index="diet", columns="metabolite", values="flux").reset_index()
ratio_df["total"] = ratio_df[KEY_SCFAS].sum(axis=1)
for met in KEY_SCFAS:
    ratio_df[f"{met}_pct"] = (ratio_df[met] / ratio_df["total"]) * 100
ratio_long = ratio_df.melt(
    id_vars=["diet"], value_vars=[m+"_pct" for m in KEY_SCFAS],
    var_name="met_pct", value_name="pct"
)
ratio_long["metabolite"] = ratio_long["met_pct"].str.replace("_pct","")

fig_ratio = px.bar(
    ratio_long,
    x="diet", y="pct", color="metabolite",
    title=tr("Relative composition A:P:B per diet", "Composición relativa A:P:B por dieta"),
    labels={"pct": tr("% of total SCFA", "% del total de SCFAs")}
)
st.plotly_chart(fig_ratio, use_container_width=True)
