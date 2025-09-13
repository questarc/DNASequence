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

st.title("DNA Sequence Analyzer for Students")

st.markdown("""
### Welcome to the DNA Analyzer!
This app is designed for 10th graders to learn about DNA sequences. Enter a DNA sequence (using A, C, G, T) and click 'Analyze' to see its features. 
You'll learn about length, GC content, complements, transcription, translation, open reading frames (ORFs), and restriction sites.
Use the samples below to get started!
""")

# Sidebar for samples
st.sidebar.header("Sample DNA Sequences")
st.sidebar.markdown("Select a sample to load it into the input area.")
selected_sample = st.sidebar.selectbox("Choose a sample:", list(samples.keys()))
if st.sidebar.button("Load Sample"):
    st.session_state.dna_input = samples[selected_sample]

# Input for DNA sequence
dna_input = st.text_area("Enter your DNA sequence (A, C, G, T only):", value=st.session_state.get('dna_input', ''), height=200)
st.markdown("*Tip: DNA is made of bases A (Adenine), C (Cytosine), G (Guanine), and T (Thymine). Sequences must only use these letters.*")

# Button to analyze
if st.button("Analyze DNA Sequence"):
    if dna_input:
        try:
            # Create Seq object
            dna_seq = Seq(dna_input.strip().upper())
            
            # Validate sequence
            if not all(base in 'ACGT' for base in dna_seq):
                st.error("Invalid DNA sequence. Please use only A, C, G, T.")
            else:
                st.success("Sequence validated successfully!")
                
                # Basic features with explanations
                with st.expander("Basic Sequence Information (Click to expand)"):
                    st.markdown("""
                    - **Length**: The total number of bases in the DNA sequence (measured in base pairs, bp).
                    - **GC Content**: The percentage of G and C bases. Higher GC content means the DNA is more stable.
                    - **Complement**: The matching strand (A pairs with T, C with G).
                    - **Reverse Complement**: The complement read backwards, useful for PCR primers.
                    - **Transcription (mRNA)**: DNA to RNA (T becomes U).
                    - **Translation (Protein)**: mRNA to protein using codons (groups of 3 bases).
                    """)
                    st.write(f"**Length:** {len(dna_seq)} bp")
                    st.write(f"**GC Content:** {GC(dna_seq):.2f}%")
                    st.write(f"**Complement:** {dna_seq.complement()}")
                    st.write(f"**Reverse Complement:** {dna_seq.reverse_complement()}")
                    st.write(f"**Transcription (mRNA):** {dna_seq.transcribe()}")
                    st.write(f"**Translation (Protein):** {dna_seq.translate(to_stop=True)}")  # Translate until first stop
                
                # Find Open Reading Frames (ORFs)
                with st.expander("Open Reading Frames (ORFs) (Click to expand)"):
                    st.markdown("""
                    ORFs are potential protein-coding regions starting with ATG (start codon) and ending with a stop codon (TAA, TAG, TGA).
                    We search in all 3 reading frames (starting at position 1, 2, or 3). Only ORFs longer than 30 bp are shown.
                    """)
                    st.write("Searching for ORFs longer than 30 bp in all frames...")
                    
                    def find_orfs(seq, min_length=30):
                        orfs = []
                        for frame in range(3):
                            for i in range(frame, len(seq), 3):
                                if seq[i:i+3] == 'ATG':  # Start codon
                                    for j in range(i+3, len(seq), 3):
                                        codon = seq[j:j+3]
                                        if len(codon) != 3:
                                            break
                                        if codon in ['TAA', 'TAG', 'TGA']:  # Stop codons
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
                            st.subheader(f"ORF in Frame {orf['frame']}")
                            st.write(f"Position: {orf['start']} - {orf['end']}")
                            st.write(f"Length: {orf['length']} bp")
                            st.write(f"DNA Sequence: {orf['sequence']}")
                            st.write(f"Protein Sequence: {orf['protein']}")
                    else:
                        st.info("No ORFs found longer than 30 bp. Try a longer sequence!")
                
                # Restriction Sites
                with st.expander("Restriction Enzyme Sites (Click to expand)"):
                    st.markdown("""
                    Restriction enzymes cut DNA at specific sequences (sites). We check for common ones like EcoRI (cuts at GAATTC).
                    Positions show where the cut starts (1-based indexing).
                    """)
                    analysis = Analysis(common_enzymes, dna_seq, linear=True)
                    sites = analysis.full()
                    
                    if any(sites.values()):
                        for enzyme, positions in sites.items():
                            if positions:
                                st.write(f"**{enzyme}:** Positions {', '.join(map(str, positions))}")
                    else:
                        st.info("No restriction sites found for common enzymes. Try sequences with known sites!")
                
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
    else:
        st.warning("Please enter a DNA sequence to analyze.")
