import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Restriction import RestrictionBatch, Analysis

# Define a set of common restriction enzymes
common_enzymes = RestrictionBatch(['EcoRI', 'BamHI', 'HindIII', 'NotI', 'XhoI'])

st.title("DNA Sequence Feature Analyzer")

# Input for DNA sequence
dna_input = st.text_area("Enter DNA sequence (A, C, G, T only):", height=200)

if dna_input:
    try:
        # Create Seq object
        dna_seq = Seq(dna_input.strip().upper())
        
        # Validate sequence
        if not all(base in 'ACGT' for base in dna_seq):
            st.error("Invalid DNA sequence. Please use only A, C, G, T.")
        else:
            st.success("Sequence validated successfully!")
            
            # Basic features
            st.header("Basic Sequence Information")
            st.write(f"**Length:** {len(dna_seq)} bp")
            st.write(f"**GC Content:** {GC(dna_seq):.2f}%")
            st.write(f"**Complement:** {dna_seq.complement()}")
            st.write(f"**Reverse Complement:** {dna_seq.reverse_complement()}")
            st.write(f"**Transcription (mRNA):** {dna_seq.transcribe()}")
            st.write(f"**Translation (Protein):** {dna_seq.translate(to_stop=True)}")  # Translate until first stop
            
            # Find Open Reading Frames (ORFs)
            st.header("Open Reading Frames (ORFs)")
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
                st.info("No ORFs found longer than 30 bp.")
            
            # Restriction Sites
            st.header("Restriction Enzyme Sites")
            analysis = Analysis(common_enzymes, dna_seq, linear=True)
            sites = analysis.full()
            
            if any(sites.values()):
                for enzyme, positions in sites.items():
                    if positions:
                        st.write(f"**{enzyme}:** Positions {', '.join(map(str, positions))}")
            else:
                st.info("No restriction sites found for common enzymes.")
            
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
