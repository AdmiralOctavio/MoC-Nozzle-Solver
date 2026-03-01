import streamlit as st
import engine.Parameters as P
import main
import matplotlib.pyplot as plt
import utils.Output as Output
import hashlib
from engine.Fuels import PROPELLANTS, PROPELLANT_OPTIONS

st.set_page_config(page_title="MOC Nozzle Designer", layout="wide")

# ── Theme state ───────────────────────────────────────────────────────────────
if "theme" not in st.session_state:
    st.session_state.theme = "lcars"

theme = st.session_state.theme

# ── Theme toggle button CSS — pinned into top-right corner ────────────────────
if theme == "lcars":
    st.markdown("""
    <style>
    div[data-testid="stButton"]:has(button[data-testid="baseButton-secondary"]#theme-btn) button,
    div[data-testid="stButton"] button[key="theme-btn"] {
        position: fixed !important;
        top: 0px !important;
        right: 0px !important;
        width: 98px !important;
        height: 74px !important;
        background-color: #b566ff !important;
        color: #000000 !important;
        border: none !important;
        border-radius: 0 32px 0 0 !important;
        font-family: 'Antonio', sans-serif !important;
        font-size: 0.78rem !important;
        font-weight: 700 !important;
        letter-spacing: 0.15em !important;
        text-transform: uppercase !important;
        z-index: 1000 !important;
        margin: 0 !important;
        padding: 0 !important;
        line-height: 1.3 !important;
    }
    div[data-testid="stButton"] button[key="theme-btn"]:hover {
        background-color: #cc88ff !important;
        box-shadow: none !important;
    }
    /* Nudge the button column to zero padding so it sits flush */
    div[data-testid="stHorizontalBlock"] > div:last-child {
        padding: 0 !important;
    }
    </style>
    """, unsafe_allow_html=True)
else:
    st.markdown("""
    <style>
    div[data-testid="stButton"] button[key="theme-btn"] {
        position: fixed !important;
        top: 0px !important;
        right: 0px !important;
        width: 98px !important;
        height: 56px !important;
        background-color: #0a0500 !important;
        color: #ff6a00 !important;
        border: 1px solid #ff6a00 !important;
        border-radius: 0 !important;
        font-family: 'Share Tech Mono', monospace !important;
        font-size: 0.72rem !important;
        font-weight: 400 !important;
        letter-spacing: 0.12em !important;
        text-transform: uppercase !important;
        z-index: 1000 !important;
        margin: 0 !important;
        padding: 0 !important;
    }
    div[data-testid="stButton"] button[key="theme-btn"]:hover {
        background-color: #ff6a00 !important;
        color: #0a0500 !important;
        box-shadow: 0 0 12px #ff6a00 !important;
    }
    </style>
    """, unsafe_allow_html=True)

# Render the toggle button (will be repositioned by CSS above)
if st.button(
    "▓  CRT\nMODE" if theme == "lcars" else "◈  LCARS\nMODE",
    key="theme-btn"
):
    st.session_state.theme = "crt" if theme == "lcars" else "lcars"
    st.rerun()


# ══════════════════════════════════════════════════════════════════════════════
# LCARS THEME
# ══════════════════════════════════════════════════════════════════════════════
if theme == "lcars":

    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Antonio:wght@400;700&display=swap');

    :root {
        --lcars-black:   #000000;
        --lcars-bg:      #0d0a08;
        --lcars-orange:  #ff9900;
        --lcars-gold:    #ffcc66;
        --lcars-peach:   #ffccaa;
        --lcars-red:     #cc4444;
        --lcars-purple:  #b566ff;
        --lcars-blue:    #4488ff;
        --lcars-teal:    #44aaaa;
        --lcars-text:    #ffffff;
        --lcars-label:   #ffcc66;
    }

    html, body,
    [data-testid="stAppViewContainer"],
    [data-testid="stApp"] {
        background-color: var(--lcars-bg) !important;
        color: var(--lcars-text) !important;
        font-family: 'Antonio', sans-serif !important;
    }

    [data-testid="stHeader"], header[data-testid="stHeader"],
    [data-testid="stSidebarCollapsedControl"],
    [data-testid="collapsedControl"],
    button[kind="headerNoPadding"],
    [data-testid="stSidebarNavCollapseIcon"],
    [data-testid="stDecoration"] { display: none !important; }

    [data-testid="stMain"], .main .block-container {
        background-color: var(--lcars-bg) !important;
        padding: 0 1.5rem 2rem 1.5rem !important;
        max-width: 100% !important;
    }

    [data-testid="stSidebar"] {
        background-color: var(--lcars-bg) !important;
        border-right: none !important;
        padding-top: 0 !important;
    }
    [data-testid="stSidebar"] > div:first-child { padding-top: 0 !important; }
    [data-testid="stSidebar"] * {
        color: var(--lcars-text) !important;
        font-family: 'Antonio', sans-serif !important;
    }
    [data-testid="stSidebar"] h1,
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3 {
        color: var(--lcars-black) !important;
        background: var(--lcars-purple) !important;
        padding: 4px 16px 4px 12px !important;
        border-radius: 0 20px 20px 0 !important;
        letter-spacing: 0.12em;
        font-size: 0.82rem !important;
        text-transform: uppercase;
        margin: 14px 0 6px 0 !important;
        display: inline-block;
    }

    input, [data-baseweb="input"] input {
        background-color: #111122 !important;
        color: var(--lcars-gold) !important;
        border: 2px solid var(--lcars-purple) !important;
        border-radius: 0 20px 20px 0 !important;
        font-family: 'Antonio', sans-serif !important;
        font-size: 1rem !important;
        caret-color: var(--lcars-orange);
    }
    input:focus {
        border-color: var(--lcars-orange) !important;
        box-shadow: 0 0 8px var(--lcars-orange) !important;
        outline: none !important;
    }
    [data-testid="stSidebar"] label,
    [data-testid="stSidebar"] p {
        color: var(--lcars-gold) !important;
        font-size: 0.78rem !important;
        text-transform: uppercase;
        letter-spacing: 0.08em;
    }

    button[data-testid="stNumberInputStepDown"],
    button[data-testid="stNumberInputStepUp"],
    [data-baseweb="spinner"] button {
        background-color: var(--lcars-purple) !important;
        color: var(--lcars-black) !important;
        border: none !important;
        border-radius: 4px !important;
        font-family: 'Antonio', sans-serif !important;
        font-weight: 700;
    }
    button[data-testid="stNumberInputStepDown"]:hover,
    button[data-testid="stNumberInputStepUp"]:hover,
    [data-baseweb="spinner"] button:hover { background-color: var(--lcars-orange) !important; }
    button[data-testid="stNumberInputStepDown"] svg,
    button[data-testid="stNumberInputStepUp"] svg,
    [data-baseweb="spinner"] button svg { fill: var(--lcars-black) !important; }

    [data-testid="stSlider"] [data-baseweb="slider"] [role="slider"] {
        background-color: var(--lcars-orange) !important;
        border-color: var(--lcars-orange) !important;
    }
    [data-testid="stToggle"] span { background-color: #333355 !important; }
    [data-testid="stToggle"] input:checked + span { background-color: var(--lcars-teal) !important; }

    /* All buttons except the theme toggle */
    [data-testid="stButton"] button:not([key="theme-btn"]) {
        background-color: var(--lcars-orange) !important;
        color: var(--lcars-black) !important;
        border: none !important;
        border-radius: 0 30px 30px 0 !important;
        font-family: 'Antonio', sans-serif !important;
        font-size: 1rem !important;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.15em;
        padding: 0.6rem 1.4rem !important;
        transition: background-color 0.15s, box-shadow 0.15s;
    }
    [data-testid="stButton"] button:not([key="theme-btn"]):hover {
        background-color: var(--lcars-gold) !important;
        box-shadow: 0 0 16px var(--lcars-orange);
    }

    [data-testid="stVerticalBlockBorderWrapper"] > div {
        background-color: #100d0a !important;
        border: none !important;
        border-left: 8px solid var(--lcars-purple) !important;
        border-radius: 0 !important;
        padding: 1rem 1.2rem !important;
        margin-bottom: 6px;
    }

    h2 {
        color: var(--lcars-black) !important;
        background: var(--lcars-teal) !important;
        font-family: 'Antonio', sans-serif !important;
        font-size: 0.88rem !important;
        letter-spacing: 0.15em;
        text-transform: uppercase;
        padding: 5px 18px 5px 12px !important;
        border-radius: 0 20px 20px 0 !important;
        display: inline-block;
        margin-bottom: 10px !important;
    }
    h3 {
        color: var(--lcars-gold) !important;
        font-family: 'Antonio', sans-serif !important;
        letter-spacing: 0.1em;
        text-transform: uppercase;
    }

    [data-testid="stAlert"] {
        border-radius: 0 20px 20px 0 !important;
        font-family: 'Antonio', sans-serif !important;
        letter-spacing: 0.06em;
    }
    [data-testid="stPyplotRootElement"] {
        border-left: 6px solid var(--lcars-teal);
        border-bottom: 2px solid var(--lcars-purple);
    }
    hr { border-color: #1a1a33 !important; }
    [data-testid="stSpinner"] p {
        color: var(--lcars-gold) !important;
        font-family: 'Antonio', sans-serif !important;
        letter-spacing: 0.1em;
    }
    ::-webkit-scrollbar { width: 6px; background: var(--lcars-bg); }
    ::-webkit-scrollbar-thumb { background: var(--lcars-purple); border-radius: 20px; }
    ::-webkit-scrollbar-thumb:hover { background: var(--lcars-orange); }
    [data-testid="stAppDeployButton"] button {
        background: var(--lcars-purple) !important;
        color: var(--lcars-black) !important;
        border: none !important;
        border-radius: 20px !important;
        font-family: 'Antonio', sans-serif !important;
    }
    </style>
    """, unsafe_allow_html=True)

    # ── Colour helpers ────────────────────────────────────────────────────────
    LCARS_PALETTE = ["#ff9900", "#ffcc66", "#b566ff", "#44aaaa", "#ffccaa", "#4488ff", "#cc6644"]

    def _seeded_color(seed_str: str, offset: int = 0) -> str:
        h = int(hashlib.md5(seed_str.encode()).hexdigest(), 16)
        return LCARS_PALETTE[(h + offset) % len(LCARS_PALETTE)]

    def _contrast_pair(seed_str: str):
        h = int(hashlib.md5(seed_str.encode()).hexdigest(), 16)
        i1 = h % len(LCARS_PALETTE)
        i2 = (h // len(LCARS_PALETTE) + 1) % len(LCARS_PALETTE)
        if i2 == i1:
            i2 = (i2 + 1) % len(LCARS_PALETTE)
        return LCARS_PALETTE[i1], LCARS_PALETTE[i2]

    def lcars_metric(label: str, value: str, color: str = None):
        if color is None:
            color = _seeded_color(label)
        st.markdown(f"""
        <div style="background:#0d0a08; border-left:6px solid {color};
                    border-bottom:2px solid {color}55;
                    padding:10px 14px 8px 14px; margin-bottom:6px;">
            <div style="font-family:'Antonio',sans-serif; font-size:0.63rem;
                        letter-spacing:0.2em; text-transform:uppercase;
                        color:{color}; margin-bottom:4px;">{label}</div>
            <div style="font-family:'Antonio',sans-serif; font-size:1.4rem;
                        font-weight:700; color:#ffffff; letter-spacing:0.04em;">{value}</div>
        </div>
        """, unsafe_allow_html=True)

    def lcars_metric_status(label: str, value: str, status_ok: bool = True):
        color      = "#44cc88" if status_ok else "#cc4444"
        badge_text = "NOMINAL" if status_ok else "WARNING"
        st.markdown(f"""
        <div style="background:#0d0a08; border-left:6px solid {color};
                    border-bottom:2px solid {color}55;
                    padding:10px 14px 8px 14px; margin-bottom:6px;">
            <div style="font-family:'Antonio',sans-serif; font-size:0.63rem;
                        letter-spacing:0.2em; text-transform:uppercase; color:{color};
                        margin-bottom:4px; display:flex; align-items:center; gap:8px;">
                {label}
                <span style="background:{color}; color:#000; padding:1px 8px;
                             border-radius:10px; font-size:0.58rem; font-weight:700;
                             letter-spacing:0.15em;">{badge_text}</span>
            </div>
            <div style="font-family:'Antonio',sans-serif; font-size:1.4rem;
                        font-weight:700; color:#ffffff; letter-spacing:0.04em;">{value}</div>
        </div>
        """, unsafe_allow_html=True)

    def lcars_divider(label: str, primary: str = "#ff9900"):
        c1, c2 = _contrast_pair(label + "A")
        c3, _  = _contrast_pair(label + "B")
        st.markdown(f"""
        <div style="display:flex; align-items:stretch; height:32px;
                    margin:20px 0 10px 0; gap:0;">
            <div style="width:48px; background:{primary};
                        border-radius:16px 0 0 16px; flex-shrink:0;"></div>
            <div style="width:6px; background:#0d0a08;"></div>
            <div style="flex:1; background:{primary}; display:flex;
                        align-items:center; padding:0 18px;">
                <span style="font-family:'Antonio',sans-serif; font-size:0.88rem;
                             font-weight:700; letter-spacing:0.22em; color:#000;
                             text-transform:uppercase; white-space:nowrap;">{label}</span>
            </div>
            <div style="width:6px; background:#0d0a08;"></div>
            <div style="width:55px; background:{c1};"></div>
            <div style="width:6px; background:#0d0a08;"></div>
            <div style="width:35px; background:{c2};"></div>
            <div style="width:6px; background:#0d0a08;"></div>
            <div style="width:70px; background:{c3}; border-radius:0 16px 16px 0;"></div>
        </div>
        """, unsafe_allow_html=True)

    # ── Top bar (purple segment omitted — replaced by fixed theme button) ─────
    st.markdown("""
    <div style="display:flex; align-items:stretch; height:64px; margin-bottom:6px; gap:0;">
        <div style="width:64px; background:#ff9900; border-radius:32px 0 0 0; flex-shrink:0;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="flex:1; background:#ff9900; display:flex; align-items:center;
                    padding:0 20px; gap:20px; white-space:nowrap; overflow:hidden;">
            <span style="font-family:'Antonio',sans-serif; font-size:1.4rem; font-weight:700;
                         letter-spacing:0.18em; color:#000; text-transform:uppercase;">
                MOC NOZZLE DESIGNER
            </span>
            <span style="font-size:0.68rem; color:#0d0a08; opacity:0.65; letter-spacing:0.1em;
                         overflow:hidden; text-overflow:ellipsis;">
                AXISYMMETRIC · METHOD OF CHARACTERISTICS · STARDATE 2401.7
            </span>
        </div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="width:55px; background:#44aaaa; flex-shrink:0;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="width:40px; background:#ffcc66; flex-shrink:0;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="width:98px; background:#b566ff; border-radius:0 32px 0 0; flex-shrink:0;"></div>
    </div>
    <div style="display:flex; height:10px; margin-bottom:16px; gap:6px;">
        <div style="width:64px; background:#b566ff; border-radius:0 0 0 6px;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="flex:3; background:#ff9900;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="flex:1; background:#44aaaa;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="width:40px; background:#ffccaa;"></div>
        <div style="width:8px; background:#0d0a08;"></div>
        <div style="width:98px; background:#ffcc66; border-radius:0 0 6px 0;"></div>
    </div>
    """, unsafe_allow_html=True)

    # ── Sidebar ───────────────────────────────────────────────────────────────
    with st.sidebar:
        st.markdown("""
        <div style="background:#b566ff; border-radius:28px 0 0 0; height:56px;
                    display:flex; align-items:center; padding:0 16px; margin-bottom:0;">
            <span style="font-family:'Antonio',sans-serif; font-size:1.1rem; font-weight:700;
                         letter-spacing:0.2em; color:#000000; text-transform:uppercase;">
                PARAMETERS
            </span>
        </div>
        <div style="display:flex; height:10px; margin-bottom:12px; gap:4px;">
            <div style="flex:2; background:#ff9900;"></div>
            <div style="width:4px; background:#0d0a08;"></div>
            <div style="flex:1; background:#44aaaa;"></div>
            <div style="width:4px; background:#0d0a08;"></div>
            <div style="flex:1; background:#ffcc66;"></div>
        </div>
        """, unsafe_allow_html=True)

        st.header("01 — PERFORMANCE")
        P.Selected_Propellant = st.selectbox(
            "Propellant",
            options=[k for k, _ in PROPELLANT_OPTIONS],
            format_func=lambda k: PROPELLANTS[k].name,
        )
        pc_mPa                 = st.number_input("Chamber Pressure (Bar)", 10, 200, 34)
        P.Ambient_P            = st.number_input("Ambient Pressure (Pa)", 0, 101325, 101325)
        of_val                 = st.number_input("O/F Ratio", value=5.13)
        thrust_val             = st.number_input("Target Thrust (N)", value=600)
        P.Combustion_Efficiency = st.number_input("Combustive Efficiency (%)", 0, 100, 85) / 100
        P.Refinement           = st.number_input("Refinement Factor", value=100)
        P.Shorten_Percentage   = st.slider("Total Nozzle Length (%)", 50, 100, 75, 5) / 100

        st.header("02 — CHAMBER")
        P.L_combustion      = st.number_input("Total Chamber Length", value=83.02)
        P.Contraction_ratio = st.number_input("Contraction Ratio", value=16)
        P.Chamber_Slope     = st.number_input("Convergent Chamber Slope", value=45)
        P.R1                = st.number_input("Radius 1", value=10)
        P.R2                = st.number_input("Radius 2", value=50)

        st.header("03 — OUTPUT")
        P.Graph2d = st.toggle("2D Characteristic Grid", value=True)
        P.Graph3d = st.toggle("3D Characteristic Grid", value=False)
        P.Stl     = st.toggle("Generate STL", value=False)
        P.Dxf     = st.toggle("Generate DXF", value=True)

        st.markdown("<div style='height:12px'></div>", unsafe_allow_html=True)

    # ── Run button ────────────────────────────────────────────────────────────
    if st.sidebar.button("▶  INITIALIZE SOLVER", use_container_width=True):
        P.P_combustion        = pc_mPa * 1e5
        P.Oxidiser_Fuel_Ratio = of_val
        P.Thrust              = thrust_val
        P.update_engine_data(P.P_combustion, P.Oxidiser_Fuel_Ratio)

        with st.spinner("PROCESSING — PLEASE STAND BY..."):
            opt_results = main.run(gui_mode=True)
            data = Output.outputTable(opt_results["rt"], opt_results["mdot"], opt_results["mach"], theme)

        lcars_divider("GEOMETRY + DIMENSIONS", "#ff9900")
        with st.container(border=True):
            c1, c2, c3 = st.columns(3)
            with c1:
                lcars_metric("Nozzle Length",   f"{data['wall_x'][-1]:.2f} mm")
                lcars_metric("Total Length",     f"{data['total_length']:.2f} mm")
            with c2:
                lcars_metric("Throat Radius",    f"{data['y_min']:.2f} mm")
                lcars_metric("Throat Diameter",  f"{data['y_min']*2:.2f} mm")
            with c3:
                lcars_metric("Exit Radius",      f"{data['exit_radius']:.2f} mm")
            lcars_metric_status("Exit Angle", f"{data['Exit_Angle']:.2f}°",
                                status_ok=data["Exit_Angle"] <= 6)

        lcars_divider("PERFORMANCE + COMBUSTION", "#b566ff")
        with st.container(border=True):
            p1, p2, p3 = st.columns(3)
            with p1:
                lcars_metric("Chamber Temperature",     f"{data['T_comb']:.1f} K")
                lcars_metric("Gamma",                   f"{data['g']:.3f}")
            with p2:
                lcars_metric("SL Optimal Pressure Ratio",  f"{data['P_comb']/101325:.2f}")
                lcars_metric("SL Optimal Exit Mach",       f"{data['M_opt']:.2f}")
                lcars_metric("SL Optimal Expansion Ratio", f"{data['AR_opt']:.2f}")
            with p3:
                lcars_metric("Predicted Exit Mach",     f"{data['M_exit_predicted']:.2f}")

        lcars_divider("THRUST + FLOW STABILITY", "#44aaaa")
        with st.container(border=True):
            t1, t2 = st.columns(2)
            with t1:
                lcars_metric("Total Predicted Thrust", f"{data['Thrust_total']:.1f} N")
            with t2:
                flow_ok = data["P_exit"] >= 0.4 * P.Ambient_P
                lcars_metric_status("Exit Pressure", f"{data['P_exit']/1000:.1f} kPa",
                                    status_ok=flow_ok)
                if not flow_ok:
                    st.error("⚠  FLOW SEPARATION RISK — RECOMMEND REDESIGN")
                else:
                    st.success("FLOW STABLE — NOMINAL OPERATING CONDITIONS")

        lcars_divider("EFFICIENCY + MASS FLOW", "#ffcc66")
        with st.container(border=True):
            e1, e2 = st.columns(2)
            with e1:
                lcars_metric("Specific Impulse (Design)", f"{data['Isp_design']:.1f} s")
                lcars_metric("CEA Ideal Isp",             f"{data['Isp_cea']:.1f} s")
            with e2:
                lcars_metric("Mass Flow Rate", f"{data['mdot']:.3f} kg/s")

        if data["fig"]:
            lcars_divider("METHOD OF CHARACTERISTICS GRID", "#4488ff")
            st.pyplot(data["fig"])

        st.markdown("""
        <div style="margin-top:28px; display:flex; height:18px; gap:6px;">
            <div style="width:64px; background:#b566ff; border-radius:0 0 0 12px;"></div>
            <div style="width:8px; background:#0d0a08;"></div>
            <div style="flex:3; background:#ff9900;"></div>
            <div style="width:8px; background:#0d0a08;"></div>
            <div style="width:40px; background:#ffccaa;"></div>
            <div style="width:8px; background:#0d0a08;"></div>
            <div style="flex:1; background:#44aaaa;"></div>
            <div style="width:8px; background:#0d0a08;"></div>
            <div style="width:98px; background:#ffcc66; border-radius:0 0 12px 0;"></div>
        </div>
        """, unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# CRT THEME
# ══════════════════════════════════════════════════════════════════════════════
elif theme == "crt":

    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Share+Tech+Mono&display=swap');

    :root {
        --crt-bg:      #0a0500;
        --crt-surface: #0f0800;
        --crt-border:  #7a3000;
        --crt-amber:   #ff6a00;
        --crt-glow:    #ffaa44;
        --crt-dim:     #7a3000;
        --crt-green:   #88ff00;
    }

    html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"] {
        background-color: var(--crt-bg) !important;
        color: var(--crt-amber) !important;
        font-family: 'Share Tech Mono', monospace !important;
    }
    [data-testid="stMain"], .main .block-container {
        background-color: var(--crt-bg) !important;
        padding-top: 1.5rem;
    }
    [data-testid="stSidebar"] {
        background-color: var(--crt-surface) !important;
        border-right: 1px solid var(--crt-border) !important;
    }
    [data-testid="stSidebar"] * {
        color: var(--crt-amber) !important;
        font-family: 'Share Tech Mono', monospace !important;
    }
    [data-testid="stSidebar"] h1,
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3 {
        color: var(--crt-glow) !important;
        border-bottom: 1px solid var(--crt-border);
        padding-bottom: 4px;
        letter-spacing: 0.1em;
    }
    h1 {
        color: var(--crt-glow) !important;
        font-family: 'Share Tech Mono', monospace !important;
        letter-spacing: 0.12em;
        text-transform: uppercase;
        border-bottom: 1px solid var(--crt-border);
        padding-bottom: 0.4rem;
    }
    h2, h3 {
        color: var(--crt-amber) !important;
        font-family: 'Share Tech Mono', monospace !important;
        letter-spacing: 0.08em;
        text-transform: uppercase;
    }
    [data-testid="stVerticalBlockBorderWrapper"] > div {
        background-color: var(--crt-surface) !important;
        border: 1px solid var(--crt-border) !important;
        border-radius: 0 !important;
        padding: 1rem !important;
    }
    [data-testid="stMetric"] {
        background-color: var(--crt-bg) !important;
        border: 1px solid var(--crt-border) !important;
        border-radius: 0 !important;
        padding: 0.5rem 0.75rem !important;
    }
    [data-testid="stMetricLabel"] {
        color: var(--crt-dim) !important;
        font-size: 0.7rem !important;
        text-transform: uppercase;
        letter-spacing: 0.1em;
    }
    [data-testid="stMetricValue"] {
        color: var(--crt-glow) !important;
        font-size: 1.4rem !important;
    }
    [data-testid="stMetricDelta"] { color: var(--crt-green) !important; }
    input, textarea, [data-baseweb="input"] input {
        background-color: var(--crt-bg) !important;
        color: var(--crt-amber) !important;
        border: 1px solid var(--crt-border) !important;
        border-radius: 0 !important;
        font-family: 'Share Tech Mono', monospace !important;
        caret-color: var(--crt-amber);
    }
    input:focus {
        border-color: var(--crt-amber) !important;
        box-shadow: 0 0 6px var(--crt-amber) !important;
        outline: none !important;
    }
    [data-testid="stSlider"] [data-baseweb="slider"] div[role="slider"] {
        background-color: var(--crt-amber) !important;
        border-color: var(--crt-amber) !important;
    }
    [data-testid="stToggle"] span { background-color: var(--crt-border) !important; }
    [data-testid="stToggle"] input:checked + span { background-color: var(--crt-amber) !important; }

    /* All buttons except theme toggle */
    [data-testid="stButton"] button:not([key="theme-btn"]) {
        background-color: var(--crt-bg) !important;
        color: var(--crt-amber) !important;
        border: 1px solid var(--crt-amber) !important;
        border-radius: 0 !important;
        font-family: 'Share Tech Mono', monospace !important;
        text-transform: uppercase;
        letter-spacing: 0.1em;
        transition: box-shadow 0.15s ease, background-color 0.15s ease;
    }
    [data-testid="stButton"] button:not([key="theme-btn"]):hover {
        background-color: var(--crt-amber) !important;
        color: var(--crt-bg) !important;
        box-shadow: 0 0 12px var(--crt-amber);
    }
    [data-testid="stSpinner"] { color: var(--crt-amber) !important; }
    hr { border-color: var(--crt-border) !important; }
    [data-testid="stPyplotRootElement"] { border: 1px solid var(--crt-border); }
    ::-webkit-scrollbar { width: 6px; background: var(--crt-bg); }
    ::-webkit-scrollbar-thumb { background: var(--crt-border); border-radius: 0; }
    ::-webkit-scrollbar-thumb:hover { background: var(--crt-amber); }
    [data-testid="stSidebarCollapsedControl"],
    [data-testid="collapsedControl"],
    button[kind="headerNoPadding"],
    [data-testid="stSidebarNavCollapseIcon"] { display: none !important; }
    [data-testid="stHeader"], header[data-testid="stHeader"] {
        background-color: var(--crt-bg) !important;
        border-bottom: 1px solid var(--crt-border) !important;
    }
    button[data-testid="stNumberInputStepDown"],
    button[data-testid="stNumberInputStepUp"] {
        background-color: var(--crt-bg) !important;
        color: var(--crt-amber) !important;
        border: 1px solid var(--crt-border) !important;
        border-radius: 0 !important;
    }
    button[data-testid="stNumberInputStepDown"]:hover,
    button[data-testid="stNumberInputStepUp"]:hover {
        background-color: var(--crt-amber) !important;
        color: var(--crt-bg) !important;
    }
    [data-testid="stNumberInputStepDown"] svg,
    [data-testid="stNumberInputStepUp"] svg { fill: var(--crt-amber) !important; }
    [data-testid="stAppDeployButton"] button {
        background: transparent !important;
        color: var(--crt-amber) !important;
        border: 1px solid var(--crt-border) !important;
        border-radius: 0 !important;
        font-family: 'Share Tech Mono', monospace !important;
    }
    [data-testid="stAppViewContainer"]::after {
        content: "";
        pointer-events: none;
        position: fixed;
        inset: 0;
        background: repeating-linear-gradient(
            to bottom,
            transparent 0px, transparent 3px,
            rgba(0,0,0,0.18) 3px, rgba(0,0,0,0.18) 4px
        );
        z-index: 9999;
    }
    </style>
    """, unsafe_allow_html=True)

    st.title("// MOC NOZZLE DESIGNER")
    st.caption("METHOD OF CHARACTERISTICS  ·  AXISYMMETRIC NOZZLE SOLVER  ·  T-001V")

    with st.sidebar:
        st.header("01 // PERFORMANCE")
        P.Selected_Propellant = st.selectbox(
            "Propellant",
            options=[k for k, _ in PROPELLANT_OPTIONS],
            format_func=lambda k: PROPELLANTS[k].name,
        )
        pc_mPa                  = st.number_input("Chamber Pressure (Bar)", 10, 200, 34)
        P.Ambient_P             = st.number_input("Ambient Pressure (Pa)", 0, 101325, 101325)
        of_val                  = st.number_input("O/F Ratio", value=5.13)
        thrust_val              = st.number_input("Target Thrust (N)", value=600)
        P.Combustion_Efficiency = st.number_input("Combustive Efficiency (%)", 0, 100, 85) / 100
        P.Refinement            = st.number_input("Refinement Factor", value=100)
        P.Shorten_Percentage    = st.slider("Total Nozzle Length", 50, 100, 75, 5) / 100

        st.header("02 // CHAMBER")
        P.L_combustion      = st.number_input("Total Chamber Length", value=83.02)
        P.Contraction_ratio = st.number_input("Contraction Ratio", value=16)
        P.Chamber_Slope     = st.number_input("Convergent Chamber Slope", value=45)
        P.R1                = st.number_input("Radius 1", value=10)
        P.R2                = st.number_input("Radius 2", value=50)

        st.header("03 // OUTPUT")
        P.Graph2d = st.toggle("Show 2D Characteristic Grid", value=True)
        P.Graph3d = st.toggle("Show 3D Characteristic Grid", value=False)
        P.Stl     = st.toggle("Generate STL", value=False)
        P.Dxf     = st.toggle("Generate DXF", value=True)

    if st.sidebar.button("▶  RUN DESIGN OPTIMIZATION", use_container_width=True):
        P.P_combustion        = pc_mPa * 1e5
        P.Oxidiser_Fuel_Ratio = of_val
        P.Thrust              = thrust_val
        P.update_engine_data(P.P_combustion, P.Oxidiser_Fuel_Ratio)

        with st.spinner("CONVERGING ON GEOMETRY..."):
            opt_results = main.run(gui_mode=True)
            data = Output.outputTable(opt_results["rt"], opt_results["mdot"], opt_results["mach"], theme)

        st.header("FINAL NOZZLE DESIGN SPECIFICATIONS")

        with st.container(border=True):
            st.subheader("GEOMETRY + DIMENSIONS")
            c1, c2, c3 = st.columns(3)
            c1.metric("Nozzle Length",   f"{data['wall_x'][-1]:.2f} mm")
            c1.metric("Total Length",    f"{data['total_length']:.2f} mm")
            c2.metric("Throat Radius",   f"{data['y_min']:.2f} mm")
            c2.metric("Throat Diameter", f"{data['y_min']*2:.2f} mm")
            c3.metric("Exit Radius",     f"{data['exit_radius']:.2f} mm")
            angle_color = "normal" if data["Exit_Angle"] <= 6 else "inverse"
            st.metric("Exit Angle", f"{data['Exit_Angle']:.2f}°", delta_color=angle_color)

        with st.container(border=True):
            st.subheader("PERFORMANCE + COMBUSTION")
            p1, p2, p3 = st.columns(3)
            p1.metric("Chamber Temperature",    f"{data['T_comb']:.1f} K")
            p1.metric("Gamma",                  f"{data['g']:.3f}")
            p2.metric("Optimal Pressure Ratio", f"{data['P_comb']/101325:.2f}")
            p2.metric("Optimal Exit Mach",      f"{data['M_opt']:.2f}")
            p2.metric("Optimal Expansion Ratio",f"{data['AR_opt']:.2f}")
            p3.metric("Predicted Exit Mach",    f"{data['M_exit_predicted']:.2f}")

        with st.container(border=True):
            st.subheader("THRUST + FLOW STABILITY")
            t1, t2 = st.columns(2)
            with t1:
                st.metric("Predicted Thrust", f"{data['Thrust_total']:.1f} N")
            with t2:
                st.metric("Exit Pressure", f"{data['P_exit']/1000:.1f} kPa")
                if data["P_exit"] < 0.4 * P.Ambient_P:
                    st.error("⚠  FLOW SEPARATION RISK HIGH")
                else:
                    st.success("✅  FLOW STABLE")

        with st.container(border=True):
            st.subheader("EFFICIENCY + MASS FLOW")
            e1, e2 = st.columns(2)
            e1.metric("Specific Impulse (Design)", f"{data['Isp_design']:.1f} s")
            e1.write(f"**CEA Ideal Isp:** {data['Isp_cea']:.1f} s")
            e2.metric("Mass Flow Rate", f"{data['mdot']:.3f} kg/s")

        if data["fig"]:
            st.divider()
            st.subheader("METHOD OF CHARACTERISTICS GRID")
            st.pyplot(data["fig"])