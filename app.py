import requests
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

st.set_page_config(page_title="DNA → Protein", page_icon="🧬")


# ---------- 1) DNA to protein ----------
def clean_dna(text):
    # Remove FASTA header lines (starting with >), spaces and line breaks
    lines = [l for l in text.splitlines() if not l.startswith(">")]
    return "".join("".join(lines).split()).upper()


def find_longest_orf(dna_text):
    dna = Seq(dna_text)
    best = ""
    for strand in (dna, dna.reverse_complement()):
        for frame in range(3):
            end = len(strand) - ((len(strand) - frame) % 3)
            protein = str(strand[frame:end].translate())
            for piece in protein.split("*"):
                if "M" in piece:
                    orf = piece[piece.index("M"):]
                    if len(orf) > len(best):
                        best = orf
    return best


# ---------- 2) 3D structure ----------
@st.cache_data(show_spinner=False)
def predict_structure(protein):
    r = requests.post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/",
        data=protein,
        timeout=120,
    )
    r.raise_for_status()
    return r.text


def mean_confidence(pdb):
    vals = [
        float(l[60:66])
        for l in pdb.splitlines()
        if l.startswith("ATOM") and l[12:16].strip() == "CA"
    ]
    avg = sum(vals) / len(vals)
    # The service returns 0-1, so convert to a 0-100 scale
    return avg * 100 if avg <= 1 else avg


# ---------- 3) Function ----------
@st.cache_data(show_spinner=False)
def blast_search(protein, top=5):
    handle = NCBIWWW.qblast("blastp", "swissprot", protein, hitlist_size=top)
    record = NCBIXML.read(handle)
    results = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        results.append(
            (round(100 * hsp.identities / hsp.align_length, 1), aln.title)
        )
    return results


# ---------- Interface ----------
st.title("🧬 DNA → Protein")
st.caption(
    "All results are scientific predictions, not confirmed facts. "
    "Use cDNA or a bacterial gene (no introns)."
)

dna_input = st.text_area("Paste your DNA sequence here", height=150)

if st.button("Show answers", type="primary"):
    st.session_state["show"] = True

if st.session_state.get("show") and dna_input.strip():
    dna = clean_dna(dna_input)

    if not set(dna) <= set("ACGTN"):
        st.error("The sequence contains characters other than A / C / G / T.")
        st.stop()

    protein = find_longest_orf(dna)
    if len(protein) < 50:
        st.error("No suitable protein was found in this sequence.")
        st.stop()

    tab1, tab2, tab3 = st.tabs(["Protein", "3D Structure", "Function"])

    with tab1:
        st.code(protein)
        st.write(f"Length: **{len(protein)}** amino acids")

    with tab2:
        if len(protein) > 400:
            st.warning(
                "The protein is longer than 400 amino acids, "
                "which the free prediction service does not accept."
            )
        else:
            try:
                with st.spinner("Predicting the 3D structure..."):
                    pdb = predict_structure(protein)
                st.metric("Average confidence (out of 100)", f"{mean_confidence(pdb):.1f}")
                view = py3Dmol.view(width=330, height=400)
                view.addModel(pdb, "pdb")
                view.setStyle({"cartoon": {"color": "spectrum"}})
                view.zoomTo()
                components.html(view._make_html(), height=420)
            except Exception as e:
                st.error(f"The structure prediction service is unavailable right now. Try again later. ({e})")

    with tab3:
        try:
            with st.spinner(
                "Searching for similar proteins... this can take a minute or more, keep this page open"
            ):
                hits = blast_search(protein)
            if not hits:
                st.info("No sufficient match found. The function is unknown.")
            for identity, title in hits:
                st.write(f"**{identity}%** — {title[:100]}")
            st.caption(
                "Above about 40% identity: the function is likely similar. "
                "Below 25%: weak similarity."
            )
        except Exception as e:
            st.error(f"The search failed. Please try again. ({e})")
