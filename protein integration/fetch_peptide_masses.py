import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def fetch_protein_sequence(uniprot_id):
    """
    Fetch the protein sequence from UniProt using the accession number.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        sequence = ''.join(fasta_data.split('\n')[1:])
        return sequence
    else:
        raise ValueError(f"Unable to fetch data for UniProt ID: {uniprot_id}")

def calculate_molecular_mass(peptide_sequence):
    """
    Calculate the molecular mass of a given peptide sequence.
    """
    analysed_seq = ProteinAnalysis(peptide_sequence)
    return analysed_seq.molecular_weight()

def get_peptide_mass(uniprot_id, peptide_sequence):
    """
    Fetch the full protein sequence using UniProt ID and calculate the 
    molecular mass of the given peptide sequence.
    """
    protein_sequence = fetch_protein_sequence(uniprot_id)
    if peptide_sequence in protein_sequence:
        peptide_mass = calculate_molecular_mass(peptide_sequence)
        return peptide_mass
    else:
        raise ValueError(f"The peptide sequence '{peptide_sequence}' is not found in the protein sequence of UniProt ID: {uniprot_id}")

# Example usage:
uniprot_id = "P01308"  # Example UniProt ID
peptide_sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKA"  # Example peptide sequence

try:
    peptide_mass = get_peptide_mass(uniprot_id, peptide_sequence)
    print(f"The molecular mass of the peptide '{peptide_sequence}' is {peptide_mass:.2f} Da")
except ValueError as e:
    print(e)
