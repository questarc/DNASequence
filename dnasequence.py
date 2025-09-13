import streamlit as st
from Bio.Seq import Seq
try:
    from Bio.SeqUtils import GC
except ImportError:
    from Bio.SeqUtils import gc_fraction
    def GC(seq):  # Fallback for newer Biopython versions
        return gc_fraction(seq) * 100

from Bio.Restriction import RestrictionBatch, Analysis

# Define a set of common restriction enzymes
common_enzymes = RestrictionBatch(['EcoRI', 'BamHI', 'HindIII', 'NotI', 'XhoI'])

# Sample DNA sequences for testing
samples = {
    "Sample 1: Basic Sequence with Single ORF": "ATGCGTACGTTAGCTTAG",
    "Sample 2: Sequence with Multiple ORFs": "ATGCGTACGTTAGCTATGCCCGGGTAATGA",
    "Sample 3: Sequence with Restriction Sites": "GAATTCATGCGTCCCGGGATCCGCTAGCTAA",
    "Sample 4: Longer Comprehensive Sequence": "ATGCGTACGTTAGCTCCCGGGATGAATTCCTCGAGGCTAGCTAAATGCCCGGGTAATGACGCGTTAG"
}

# Custom CSS for professional styling
st.markdown("""
    <style>
    .main-header {
        color: #1E3A8A; /* Dark blue for main headers */
        font-size: 32px;
        font-weight: bold;
        margin-bottom: 10px;
    }
    .sub-header {
        color: #3B82F6; /* Lighter blue for subheaders */
        font-size: 24px;
        font-weight: bold;
        margin-top: 20px;
        margin-bottom: 10px;
    }
    .info-box {
        background-color: #F0F7FF; /* Light blue background for info */
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 20px;
    }
    .stButton>button {
        background-color: #2563EB; /* Blue button */
        color: white;
        font-weight: bold;
        padding: 10px 20px;
        border-radius: 5px;
    }
    .stButton>button:hover {
        background-color: #1E40AF; /* Darker blue on hover */
    }
    .stTextArea textarea {
        font-family: 'Courier New', Courier, monospace;
        font-size: 16px;
    }
    </style>
""", unsafe_allow_html=True)

# App title
st.markdown('<div class="main-header">DNA Sequence Analyzer for Students</div>', unsafe_allow_html=True)

# Introduction
st.markdown("""
<div class="info-box">
Welcome to the DNA Sequence Analyzer! This tool helps 10th graders explore DNA sequences. Enter a sequence using A, C, G, T, or try our samples. Click 'Analyze' to discover features like length, GC content, transcription, translation, open reading frames (ORFs), and restriction enzyme sites. Learn how DNA works in a fun, interactive way!
</div>
""", unsafe_allow_html=True)

# Sidebar for sample DNA sequences
st.sidebar.markdown('<div class="sub-header">Sample DNA Sequences</div>', unsafe_allow_html=True)
st.sidebar.markdown("Select a sample to load it into the input box.")
selected_sample = st.sidebar.selectbox("Choose a sample:", list(samples.keys()))
if st.sidebar.button("Load Sample"):
    st.session_state.dna_input = samples[selected_sample]

# Input for DNA sequence
st.markdown('<div class="sub-header">Enter Your DNA Sequence</div>', unsafe_allow_html=True)
st.markdown("*DNA uses bases: A (Adenine), C (Cytosine), G (Guanine), T (Thymine). Use only these letters.*")
dna_input = st.text_area("", value=st.session_state.get('dna_input', ''), height=150, placeholder="e.g., ATGCGTACGTTAGCTTAG")
analyze_button = st.button("Analyze DNA Sequence")

# Process DNA sequence on button click
if analyze_button:
    if dna_input:
        try:
            # Create Seq object
            dna_seq = Seq(dna_input.strip().upper())
            
            # Validate sequence
            if not all(base in 'ACGT' for base in dna_seq):
                st.error("Invalid DNA sequence! Please use only A, C, G, T.")
            else:
                st.success("Sequence validated successfully!")
                
                # Basic Sequence Information
                with st.expander("Basic Sequence Information"):
                    st.markdown("""
                    - **Length**: Total number of bases (in base pairs, bp).
                    - **GC Content**: Percentage of G and C bases (affects DNA stability).
                    - **Complement**: Matching strand (A↔T, C↔G).
                    - **Reverse Complement**: Complement read backwards (used in PCR).
                    - **Transcription (mRNA)**: DNA to RNA (T becomes U).
                    - **Translation (Protein)**: mRNA to protein using codons (3 bases).
                    """)
                    st.write(f"**Length**: {len(dna_seq)} bp")
                    st.write(f"**GC Content**: {GC(dna_seq):.2f}%")
                    st.write(f"**Complement**: {dna_seq.complement()}")
                    st.write(f"**Reverse Complement**: {dna_seq.reverse_complement()}")
                    st.write(f"**Transcription (mRNA)**: {dna_seq.transcribe()}")
                    st.write(f"**Translation (Protein)**: {dna_seq.translate(to_stop=True)}")
                
                # Open Reading Frames (ORFs)
                with st.expander("Open Reading Frames (ORFs)"):
                    st.markdown("""
                    ORFs are potential protein-coding regions starting with ATG (start codon) and ending with TAA, TAG, or TGA (stop codons). We check all 3 reading frames (starting at base 1, 2, or 3). Only ORFs ≥30 bp are shown.
                    """)
                    def find_orfs(seq, min_length=30):
                        orfs = []
                        for frame in range(3):
                            for i in range(frame, len(seq), 3):
                                if seq[i:i+3] == 'ATG':
                                    for j in range(i+3, len(seq), 3):
                                        codon = seq[j:j+3]
                                        if len(codon) != 3:
                                            break
                                        if codon in ['TAA', 'TAG', 'TGA']:
                                            orf = seq[i:j+3]
                                            if len(orf) >= min_length:
                                                orfs.append({
                                                    'frame': frame + 1,
                                                    'start': i + 1,
                                                    'end': j + 3,
                                                    'length': len(orf),
                                                    'sequence': orf,
                                                    'protein': orf.translate(to_stop=True)
                                                })
                                            break
                        return orfs
                    
                    orfs = find_orfs(dna_seq)
                    if orfs:
                        for orf in orfs:
                            st.markdown(f'<div class="sub-header">ORF in Frame {orf["frame"]}</div>', unsafe_allow_html=True)
                            st.write(f"**Position**: {orf['start']} - {orf['end']}")
                            st.write(f"**Length**: {orf['length']} bp")
                            st.write(f"**DNA Sequence**: {orf['sequence']}")
                            st.write(f"**Protein Sequence**: {orf['protein']}")
                    else:
                        st.info("No ORFs found ≥30 bp. Try a longer sequence!")
                
                # Restriction Enzyme Sites
                with st.expander("Restriction Enzyme Sites"):
                    st.markdown("""
                    Restriction enzymes cut DNA at specific sequences. For example, EcoRI cuts at GAATTC. We check for common enzymes, showing where they cut (positions are 1-based).
                    """)
                    analysis = Analysis(common_enzymes, dna_seq, linear=True)
                    sites = analysis.full()
                    
                    if any(sites.values()):
                        for enzyme, positions in sites.items():
                            if positions:
                                st.write(f"**{enzyme}**: Positions {', '.join(map(str, positions))}")
                    else:
                        st.info("No restriction sites found. Try a sequence with sites like GAATTC (EcoRI)!")
                
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
    else:
        st.warning("Please enter a DNA sequence to analyze.")
