import streamlit as st
from Bio import SeqIO, Entrez
from io import StringIO

Entrez.email = "ochipash@gmail.com"

@st.cache_data
def load_words(mode='clean', min_length=3):
    word_file = 'naughty_words.txt' if mode == 'dirty' else 'english_words.txt'
    with open(word_file, encoding='utf-8') as f:
        words = set(w.strip().lower() for w in f if len(w.strip()) >= min_length)
    return words

def fetch_fasta(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        data = handle.read()
        handle.close()
        return data
    except Exception:
        return None

def fetch_proteome_fasta(assembly_accession):
    try:
        links = Entrez.elink(dbfrom="assembly", db="protein", id=assembly_accession, linkname="assembly_protein")
        record = Entrez.read(links)
        if not record[0]["LinkSetDb"]:
            return None
        protein_ids = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
        handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()
        return fasta
    except Exception:
        return None
        
def fetch_sequence_from_accession(accession):
    accession = accession.strip()
    if accession.startswith("GCF_") or accession.startswith("GCA_"):
        return fetch_proteome_fasta(accession)
    else:
        return fetch_fasta(accession)


def extract_sequence(fasta_str):
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    return "".join(str(r.seq).lower() for r in records)

def find_words(seq, word_set):
    return sorted({word for word in word_set if word in seq})

st.set_page_config(page_title="WORDOME", layout="centered")

st.title("🔬 WORDOME")
st.markdown("Find hidden English (or naughty) words in protein sequences.")

with st.sidebar:
    st.header("⚙️ Options")
    min_length = st.slider("Minimum word length", 3, 12, 3)
    mode = st.checkbox("😈 Enable Dirty Mode")

if mode:
    st.markdown("""
        <style>
        body {
            background-color: #1e1e1e;
            color: #ff5e79;
        }
        </style>
    """, unsafe_allow_html=True)

word_set = load_words(mode='dirty' if mode else 'clean', min_length=min_length)

upload = st.file_uploader("Upload a FASTA file", type="fasta")
accession = st.text_input("Or paste an NCBI protein accession (e.g. WP_000000001)")

sequence = ""

if upload is not None:
    try:
        raw = upload.read().decode("utf-8")
        sequence = extract_sequence(raw)
    except Exception:
        st.error("Failed to parse uploaded FASTA file.")
elif accession:
    fetched = fetch_sequence_from_accession(accession.strip())
    if fetched:
        sequence = extract_sequence(fetched)
    else:
        st.error("Couldn't fetch that accession from NCBI.")

if sequence:
    st.success("Sequence loaded. Searching for words...")
    results = find_words(sequence, word_set)
    if results:
        st.write(f"Found {len(results)} words:")
        st.text("\n".join(results))
    else:
        st.warning("No words found in this sequence.")

st.markdown("---")
st.caption("Made with 💙 by your local microbe whisperer.")
