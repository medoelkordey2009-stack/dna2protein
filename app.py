import requests
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

st.set_page_config(page_title="DNA → Protein", page_icon="🧬")


# ---------- 1) من DNA لبروتين ----------
def clean_dna(text):
    # نشيل سطر الـ FASTA (اللي بيبدأ بـ >) والمسافات والأسطر
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


# ---------- 2) الشكل ثلاثي الأبعاد ----------
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
    return sum(vals) / len(vals)


# ---------- 3) الوظيفة ----------
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


# ---------- الواجهة ----------
st.title("🧬 DNA → Protein")
st.caption(
    "كل النتايج توقعات علمية، مش حقايق مؤكدة. "
    "استخدم cDNA أو جين بكتيري (من غير introns)."
)

dna_input = st.text_area(
    "الصق تسلسل DNA هنا (بعد اللصق اضغط برّه الخانة وهيبدأ التحليل لوحده)",
    height=150,
)

if dna_input.strip():
    dna = clean_dna(dna_input)

    if not set(dna) <= set("ACGTN"):
        st.error("فيه حروف مش A / C / G / T في التسلسل.")
        st.stop()

    protein = find_longest_orf(dna)
    if len(protein) < 50:
        st.error("ملقيتش بروتين مناسب في التسلسل ده.")
        st.stop()

    tab1, tab2, tab3 = st.tabs(["البروتين", "الشكل 3D", "الوظيفة"])

    with tab1:
        st.code(protein)
        st.write(f"الطول: **{len(protein)}** حمض أميني")

    with tab2:
        if len(protein) > 400:
            st.warning("البروتين أطول من 400 حمض أميني، والخدمة المجانية مش هتقبله.")
        else:
            try:
                with st.spinner("بتوقع الشكل..."):
                    pdb = predict_structure(protein)
                st.metric("متوسط الثقة (من 100)", f"{mean_confidence(pdb):.1f}")
                view = py3Dmol.view(width=330, height=400)
                view.addModel(pdb, "pdb")
                view.setStyle({"cartoon": {"color": "spectrum"}})
                view.zoomTo()
                components.html(view._make_html(), height=420)
            except Exception as e:
                st.error(f"خدمة توقع الشكل مش شغالة دلوقتي. جرب بعد شوية. ({e})")

    with tab3:
        try:
            with st.spinner(
                "بدوّر على بروتينات مشابهة... ممكن ياخد دقيقة أو أكتر، سيب الصفحة مفتوحة"
            ):
                hits = blast_search(protein)
            if not hits:
                st.info("مفيش تطابق كافي. الوظيفة غير معروفة.")
            for identity, title in hits:
                st.write(f"**{identity}%** — {title[:100]}")
            st.caption(
                "فوق 40% تقريباً: الوظيفة غالباً مشابهة. أقل من 25%: تشابه ضعيف."
            )
        except Exception as e:
            st.error(f"البحث فشل. جرب تاني. ({e})")
