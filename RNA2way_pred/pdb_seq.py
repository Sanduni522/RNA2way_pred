from glob import glob
from biopandas.pdb import PandasPdb
import sys
import pandas as pd

def sequence(pdb_file_path):
    ppdb = PandasPdb().read_pdb(pdb_file_path)
    atom_df = ppdb.df['ATOM']
    ter_df = ppdb.df['OTHERS'][ppdb.df['OTHERS']['record_name'] == 'TER']

    c2_atoms = atom_df[atom_df['atom_name'] == 'C2'].copy()

    c2_atoms['is_ter'] = False
    ter_df['is_ter'] = True
    combined_df = pd.concat([c2_atoms, ter_df], ignore_index=True)
    combined_df.sort_values(by=['line_idx'], inplace=True)

    last_residue_number = None
    last_chain_id = None
    sequences = []
    current_sequence = []

    for idx, row in combined_df.iterrows():
        if row['is_ter']:
            if current_sequence:
                sequences.append(''.join(current_sequence))
                current_sequence = []
            continue

        chain_id = row['chain_id']
        residue_number = row['residue_number']

        if chain_id != last_chain_id or (
                last_residue_number is not None and (residue_number - last_residue_number > 1)):
            if current_sequence:
                sequences.append(''.join(current_sequence))
                current_sequence = []
            last_chain_id = chain_id

        current_sequence.append(row['residue_name'])
        last_residue_number = residue_number

    if current_sequence:
        sequences.append(''.join(current_sequence))

    rna_seq = '_'.join(sequences)
    return rna_seq

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print(sequence(sys.argv[1]))
