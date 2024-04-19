import pandas as pd
import glob
import os
from pathlib import Path
import numpy as np
import sys
import rmsd

def pdb_selres():
    df = pd.read_csv("pdb_structures.csv")
    for i, row in df.iterrows():
        pdb = Path(row['pdb_name']).stem
        os.system(f"pdb_tidy pdb/{pdb}/{pdb}.pdb | pdb_selres -4:10,17,18 > pdb/{pdb}/{pdb}.pdb")


def main(a,path,pdb_file, x1, x2, x3, x4):
    df = pd.read_csv(f"{path}/bin/pdb_structures.csv")

    RMSD = []
    pdbs = df['pdb_name']
    all_coords = []
    for pdb in pdbs:
        all_coords.append(rmsd.get_coordinates_pdb(f"{path}/test/resources/{pdb_file[:-4]}/pdb/{pdb}/{pdb}.pdb", [x1, x2, x3, x4]))

    all_coords_native = [rmsd.get_coordinates_pdb(a, [x1, x2, x3, x4])]

    for i in range(0, len(all_coords)):
        rms = rmsd.compute_rmsd_from_coords(np.array(all_coords[i]), np.array(all_coords_native[0]))

        print(rms)
        RMSD.append(rms)
    df['rms'] = RMSD

    df.to_csv(f"{path}/test/resources/{pdb_file[:-4]}/pdb_structures_rms.csv", index=False)


if __name__ == '__main__':
    main(sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
