import re
import pandas as pd
import os
import glob
import freesasa
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import load_model
import math
from math import log
from keras.callbacks import EarlyStopping
import pdb_seq
import subprocess
from xlwt import Workbook
import build_pdb_set
import xlrd
import openpyxl
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy import stats

def prep_fasta_secstruct(pdb_file,path):
    os.chdir(f'{path}/test/resources/')
    subprocess.call(['mkdir', f'{pdb_file[:-4]}'])
    subprocess.call(['mv', pdb_file, f'{pdb_file[:-4]}/'])
    os.chdir(f'{path}/test/resources/{pdb_file[:-4]}')
    motif = pdb_seq.sequence(pdb_file)
    seq_low = motif.lower()
    print(seq_low)
    seq1 = seq_low.split('_')[0]
    seq2 = seq_low.split('_')[1]
    ss1 = '(' + '.' * (len(seq1) - 2) + '('
    ss2 = ')' + '.' * (len(seq2) - 2) + ')'
    fasta = open(f'{seq_low}.fasta', 'w')
    secstruct = open(f'{seq_low}.secstruct', 'w')
    ## additional gc pairs are added for better stability
    text1 = f""">>{seq_low}.pdb
gg{seq1}cc gg{seq2}cc"""
    fasta.write(text1)
    fasta.close()

    text2 = f"""(({ss1}(( )){ss2}))
gg{seq1}cc gg{seq2}cc"""
    secstruct.write(text2)
    secstruct.close()

def renumber_pdb(pdb_file,name): ## name would be the motif (ex: cccg_cccg). The name must be simple letters
    motif = pdb_seq.sequence(pdb_file)
    len_motif1 = len(motif.split('_')[0])
    end1 = len_motif1 + 3 - 1
    start2 = end1 + 5
    len_motif2 = len(motif.split('_')[1])
    end2 = start2 + len_motif2 - 1
    subprocess.call(['renumber_pdb_in_place.py', pdb_file, f'3-{end1}', f'{start2}-{end2}'])
    subprocess.call(['rna_denovo.macosclangrelease', '-fasta', f'{name}.fasta', '-secstruct_file', f'{name}.secstruct',
                     '-minimize_rna', 'true', '-s', pdb_file])
    subprocess.call(['extract_pdbs.macosclangrelease', '-in:file:silent', 'default.out'])
    subprocess.call(['mv', 'S_000001.pdb', f'{name}.pdb']) ## the default pdb name would be S_000001.pdb when extracted from the above code
    subprocess.call(['rm', '-rf', 'default.out'])
    print('renumber complete')

def farfar(name,pdb_file):
    subprocess.call(['rna_denovo.macosclangrelease', '-fasta', f'{name}.fasta', '-secstruct_file', f'{name}.secstruct',
                     '-minimize_rna', 'true', '-nstruct', '3', '-native', f'{name}.pdb', '-exclude_native_fragments',
                     'true'])
    #                 'true', '-fragment_homology_rmsd', '2.0'])
    subprocess.call(['mkdir', 'pdb'])
    subprocess.call(['extract_pdbs.macosclangrelease', '-in:file:silent', 'default.out'])
    subprocess.call('mv S_*.pdb pdb/', shell=True)

def run_dssr(pdb_file,path):
    os.chdir(f'{path}/test/resources/{pdb_file[:-4]}/pdb')
    filenames = sorted(glob.glob('*.pdb'))
    for model_pdb in filenames:
        subprocess.call(['mkdir', f'{model_pdb[:-4]}'])
        subprocess.call(['mv', f'{model_pdb}', f'{model_pdb[:-4]}/'])
        subprocess.call(['x3dna-dssr', f'-i={model_pdb[:-4]}/{model_pdb}', f'-o={model_pdb[:-4]}/{model_pdb[:-4]}.out'])
        subprocess.call(['mv', 'dssr-torsions.txt', f'{model_pdb[:-4]}/'])
        subprocess.call('rm -rf dssr-*', shell=True)

def read_dssr_file(pdb_file,path):
    os.chdir(f'{path}/test/resources/{pdb_file[:-4]}/')
    wb1 = Workbook()
    var_holder = {}
    motif1 = pdb_seq.sequence(pdb_file)
    print(motif1)
    motif11_1 = motif1.split('_')[0]
    motif11_2 = motif1.split('_')[1]
    n = len(motif1) + 7
    k = n - 5
    x1 = len(motif11_1)

    ## This half of the code is for the first strand in the motif
    for i in range(x1):
        var_holder['sheet' + str(i)] = wb1.add_sheet(f'Sheet {i + 3}', cell_overwrite_ok=True)

        x = var_holder['sheet' + str(i)]

        column_headers = ['nt', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'e-z', 'chi', 'phase-angle',
                          'sugar-type', 'ssZp', 'Dp', 'splay', 'eta', 'theta', 'eta_1', 'theta_1', 'eta_2',
                          'theta_2', 'v0', 'v1', 'v2', 'v3', 'v4', 'tm', 'P', 'Puckering', 'bin', 'cluster',
                          'suitness', 'name', 'motif', 'nt-num']
        for j in range(34):
            x.write(0, j, column_headers[j])
            names = sorted(glob.glob("pdb/S_*"))
            out = 'dssr-torsions.txt'
            row = 1
            for name in names:
                name1 = name.split('/')[1]
                lines = open(f"{name}/{out}", 'r').read().splitlines()
                start_index = -1
                start_index1 = -1
                start_index2 = -1
                start_index3 = -1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'splay' in line:  # Checks for the header line
                        start_index = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'theta' in line:  # Checks for the header line
                        start_index1 = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'Puckering' in line:  # Checks for the header line
                        start_index2 = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'suiteness' in line:
                        start_index3 = index + 1

                txt1 = lines[start_index + i + 2]
                txt2 = lines[start_index1 + i + 2]
                txt3 = lines[start_index2 + i + 2]
                txt4 = lines[start_index3 + i + 2]

                nt = txt1[9:10]
                nt_num = (i + 3)
                alpha = txt1[26:32]
                beta = txt1[34:40]
                gamma = txt1[43:48]
                delta = txt1[52:56]
                epsilon = txt1[58:65]
                zeta = txt1[67:73]
                e_z = txt1[75:83]
                chi = txt1[85:98]
                phase = txt1[101:116]
                sugar = txt1[117:127]
                ssZp = txt1[130:135]
                Dp = txt1[138:143]
                splay = txt1[145:151]
                eta = txt2[27:33]
                theta = txt2[34:41]
                eta_1 = txt2[42:49]
                theta_1 = txt2[50:57]
                eta_2 = txt2[58:65]
                theta_2 = txt2[66:73]
                v0 = txt3[28:33]
                v1 = txt3[35:41]
                v2 = txt3[44:49]
                v3 = txt3[51:57]
                v4 = txt3[59:64]
                tm = txt3[67:72]
                P = txt3[75:80]
                Puck = txt3[81:90]
                bin = txt4[25:31]
                clus = txt4[34:38]
                suit = txt4[43:53]

                x.write(row, 31, name1)
                x.write(row, 32, motif1)
                x.write(row, 0, nt)
                x.write(row, 33, nt_num)
                x.write(row, 1, alpha)
                x.write(row, 2, beta)
                x.write(row, 3, gamma)
                x.write(row, 4, delta)
                x.write(row, 5, epsilon)
                x.write(row, 6, zeta)
                x.write(row, 7, e_z)
                x.write(row, 8, chi)
                x.write(row, 9, phase)
                x.write(row, 10, sugar)
                x.write(row, 11, ssZp)
                x.write(row, 12, Dp)
                x.write(row, 13, splay)
                x.write(row, 14, eta)
                x.write(row, 15, theta)
                x.write(row, 16, eta_1)
                x.write(row, 17, theta_1)
                x.write(row, 18, eta_2)
                x.write(row, 19, theta_2)
                x.write(row, 20, v0)
                x.write(row, 21, v1)
                x.write(row, 22, v2)
                x.write(row, 23, v3)
                x.write(row, 24, v4)
                x.write(row, 25, tm)
                x.write(row, 26, P)
                x.write(row, 27, Puck)
                x.write(row, 28, bin)
                x.write(row, 29, clus)
                x.write(row, 30, suit)

                row += 1

    ## For the second strand of the motif.
    for i in range((x1 + 6), (n - 2)):
        var_holder['sheet' + str(i)] = wb1.add_sheet(f'Sheet {i + 1}', cell_overwrite_ok=True)

        x = var_holder['sheet' + str(i)]

        column_headers = ['nt', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'e-z', 'chi', 'phase-angle',
                          'sugar-type', 'ssZp', 'Dp', 'splay', 'eta', 'theta', 'eta_1', 'theta_1', 'eta_2',
                          'theta_2', 'v0', 'v1', 'v2', 'v3', 'v4', 'tm', 'P', 'Puckering', 'bin', 'cluster',
                          'suitness', 'name', 'motif', 'nt-num']

        for j in range(34):
            x.write(0, j, column_headers[j])
            names = sorted(glob.glob("pdb/S_*"))
            out = 'dssr-torsions.txt'
            row = 1
            for name in names:
                name1 = name.split('/')[1]
                lines = open(f"{name}/{out}", 'r').read().splitlines()
                start_index = -1
                start_index1 = -1
                start_index2 = -1
                start_index3 = -1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'splay' in line:  # Checks for the header line
                        start_index = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'theta' in line:  # Checks for the header line
                        start_index1 = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'Puckering' in line:  # Checks for the header line
                        start_index2 = index + 1

                for index, line in enumerate(lines):
                    if 'nt' in line and 'suiteness' in line:
                        start_index3 = index + 1
                      
                txt1 = lines[start_index + i]
                txt2 = lines[start_index1 + i]
                txt3 = lines[start_index2 + i]
                txt4 = lines[start_index3 + i]

                nt = txt1[9:10]
                nt_num = i + 1
                alpha = txt1[26:32]
                beta = txt1[34:40]
                gamma = txt1[43:48]
                delta = txt1[52:56]
                epsilon = txt1[58:65]
                zeta = txt1[67:73]
                e_z = txt1[75:83]
                chi = txt1[85:98]
                phase = txt1[101:116]
                sugar = txt1[117:127]
                ssZp = txt1[130:135]
                Dp = txt1[138:143]
                splay = txt1[145:151]
                eta = txt2[27:33]
                theta = txt2[34:41]
                eta_1 = txt2[42:49]
                theta_1 = txt2[50:57]
                eta_2 = txt2[58:65]
                theta_2 = txt2[66:73]
                v0 = txt3[28:33]
                v1 = txt3[35:41]
                v2 = txt3[44:49]
                v3 = txt3[51:57]
                v4 = txt3[59:64]
                tm = txt3[67:72]
                P = txt3[75:80]
                Puck = txt3[81:90]
                bin = txt4[25:31]
                clus = txt4[34:38]
                suit = txt4[43:53]

                x.write(row, 31, name1)
                x.write(row, 32, motif1)
                x.write(row, 0, nt)
                x.write(row, 33, nt_num)
                x.write(row, 1, alpha)
                x.write(row, 2, beta)
                x.write(row, 3, gamma)
                x.write(row, 4, delta)
                x.write(row, 5, epsilon)
                x.write(row, 6, zeta)
                x.write(row, 7, e_z)
                x.write(row, 8, chi)
                x.write(row, 9, phase)
                x.write(row, 10, sugar)
                x.write(row, 11, ssZp)
                x.write(row, 12, Dp)
                x.write(row, 13, splay)
                x.write(row, 14, eta)
                x.write(row, 15, theta)
                x.write(row, 16, eta_1)
                x.write(row, 17, theta_1)
                x.write(row, 18, eta_2)
                x.write(row, 19, theta_2)
                x.write(row, 20, v0)
                x.write(row, 21, v1)
                x.write(row, 22, v2)
                x.write(row, 23, v3)
                x.write(row, 24, v4)
                x.write(row, 25, tm)
                x.write(row, 26, P)
                x.write(row, 27, Puck)
                x.write(row, 28, bin)
                x.write(row, 29, clus)
                x.write(row, 30, suit)

                row += 1

    wb1.save(f'dssr_data.xls')

def cal_hbond(pdb_file,path):
    os.chdir(f'{path}/test/resources/{pdb_file[:-4]}/pdb/')
    subprocess.call(['mkdir', '../hbonds'])
    filenames = sorted(glob.glob('*/*.pdb'))
    for m in filenames:
        model_pdb = m.split('/')[0]
        subprocess.call(['x3dna-dssr', f'-i={model_pdb}/{model_pdb}.pdb', '--get-hbond', f'-o={model_pdb}/{model_pdb}_FARFAR-hbonds.txt'])
        subprocess.call(['mv', f'{model_pdb}/{model_pdb}_FARFAR-hbonds.txt', '../hbonds/'])

def convert_txt_to_csv(pdb_file,path):
    os.chdir(f'{path}/test/resources/{pdb_file[:-4]}/hbonds')
    filenames = sorted(glob.glob('*_FARFAR-hbonds.txt'))
    for model_pdb in filenames:
        read_file = pd.read_csv(f'{model_pdb}', skiprows=2, delimiter="\s+",
                            names=["Position_1", "Position_2", "Hbond_num", "Type", "Distance", "Hbond_atoms", "Atom_1",
                                   "Atom_2"])
        read_file.to_csv(f"{model_pdb[:-4]}.csv", index=False)

def rmsd(pdb_file,path):
    filenames = sorted(glob.glob(f'{path}/test/resources/{pdb_file[:-4]}/pdb/*/*.pdb'))
    for m in filenames:
        model_pdb = m.split('/')[-2]
        os.chdir(f'{path}/test/resources/{pdb_file[:-4]}/pdb/{model_pdb}/')
        motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}')
        motif1 = motif.split('_')[0]
        motif2 = motif.split('_')[1]
        x1 = 3
        x2 = x1 + len(motif1) - 1
        x3 = x2 + 5
        x4 = x3 + len(motif2) - 1
        build_pdb_set.main(m, path, pdb_file, x1, x2, x3, x4)

def add_RMSD(pdb_file,path):
    df_score = pd.read_csv(f'{path}/test/resources/{pdb_file[:-4]}/pdb_structures_rms.csv')
    rmsd = df_score['rms']

    xls = pd.ExcelFile(f'dssr_data.xls')
    for sheet_name in xls.sheet_names:
        df_final = pd.read_excel(f'dssr_data.xls',sheet_name=sheet_name)
        df_final['RMSD'] = rmsd

        df_final.to_excel(f'{path}/test/resources/{pdb_file[:-4]}/dssr_data_{sheet_name}_plus_rmsd.xlsx', sheet_name=sheet_name, index=False)
        print(df_final)

def multiplyList(myList, char):
    # Multiply elements one by one
    result = char
    fin_list = []
    for x in myList:
        fin_list.append(f'{result}{x}')
    return fin_list

def SASA(pdb_file,path):
    motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}')
    motif1 = motif.split('_')[0].lower()
    motif2 = motif.split('_')[1].lower()

    a_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'a']
    a_pos2 = [(pos + len(motif1) + 7) for pos, char in enumerate(motif2) if char == 'a']
    c_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'c']
    c_pos2 = [(pos + len(motif1) + 7) for pos, char in enumerate(motif2) if char == 'c']

    a_pos = a_pos1 + a_pos2
    c_pos = c_pos1 + c_pos2
    pos = a_pos + c_pos

    for l in pos:
        df = pd.read_excel(f'{path}/test/resources/{pdb_file[:-4]}/dssr_data_Sheet {l}_plus_rmsd.xlsx')  # name of the file that has list of models

        struct = df['name']
        nt_num = df['nt-num']
        index = len(df.index)

        path1 = f'{path}/test/resources/{pdb_file[:-4]}/pdb'  # the models are inside this folder

        sasa_n1 = []
        sasa_n3 = []

        sasa_final = []

        for j in range(index):
            ppdb = PandasPdb()
            stru = struct[j]
            ppdb.read_pdb(f'{path1}/{stru}/{stru}.pdb')
            ATOM = ppdb.df['ATOM']
            resname = ATOM['residue_name']
            atom = ATOM['atom_name']
            resi_number = ATOM['residue_number']
            length = len(ATOM.index)

            for i in range(length):
                if (resi_number[i] == nt_num[j]):
                    if atom[i] == 'N1':
                        if (resname[i] == 'A' or resname[i] == 'C'):
                            structure = freesasa.Structure(f'{path1}/{stru}/{stru}.pdb')
                            result = freesasa.calc(structure, freesasa.Parameters(
                                {'algorithm': freesasa.LeeRichards, 'probe-radius': 2.0}))
                            selection = freesasa.selectArea((f'n1,(name N1) and (resn A) and (resi {resi_number[i]}) ',
                                                             f'n3, (name N3) and (resn C) and (resi {resi_number[i]})'),
                                                            structure, result)
                            sel1 = selection['n1']
                            sel2 = selection['n3']
                            sasa = max(sel1, sel2)
                            print(struct[j], i, sel1, sel2, sasa)
                            sasa_n1.append(sel1)
                            sasa_n3.append(sel2)
                            sasa_final.append(sasa)
                        else:
                            sasa_n1.append('0')
                            sasa_n3.append('0')
                            sasa_final.append('0')
                            print(struct[j], i, 0, 0, 0)

        df['SASA_2'] = sasa_final

        df.to_csv(f'{path}/test/resources/{pdb_file[:-4]}/structural_parameters_{l}.csv', index=False)

def hbond(pdb_file,path):
    motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}')
    var_holder = {}
    p1 = []
    p2 = []
    p3 = []
    p4 = []
    p5 = []
    p6 = []
    p7 = []
    p8 = []
    p9 = []
    p10 = []
    p11 = []
    df_h = pd.DataFrame(
        columns=['Name', 'Atom_1', 'Nuc_1', 'Atom_2', 'Nuc_2', 'H-bond length', 'Type', 'N or other', 'Hbond_atoms',
                 'Hbond_strength', 'Angle'])
    motif1 = motif.split('_')[0]
    motif2 = motif.split('_')[1]
    a1 = motif1.count('A')
    a2 = motif2.count('A')
    c1 = motif1.count('C')
    c2 = motif2.count('C')
    a_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'A']
    a_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'A']
    c_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'C']
    c_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'C']
    pos1 = a_pos1 + a_pos2 + c_pos1 + c_pos2
    a_pos1_1 = []
    a_pos2_1 = []
    c_pos1_1 = []
    c_pos2_1 = []
    if a_pos1 != [] or a_pos2 != [] or c_pos1 != [] or c_pos2 != []:
        a_pos1_1 = multiplyList(a_pos1, 'A')
        a_pos2_1 = multiplyList(a_pos2, 'A')
        c_pos1_1 = multiplyList(c_pos1, 'C')
        c_pos2_1 = multiplyList(c_pos2, 'C')
    a_pos = a_pos1_1 + a_pos2_1
    c_pos = c_pos1_1 + c_pos2_1
    pos = a_pos + c_pos

    for n, l in zip(pos1, pos):
        print(n, l)
        var_holder['n' + str(n)] = l
        if l[0] == "A":
            var_holder['a' + str(n)] = "N1"
        else:
            var_holder['a' + str(n)] = "N3"

        n1 = var_holder['n' + str(n)]
        a1 = var_holder['a' + str(n)]
        name = glob.glob(f'{path}/test/resources/{pdb_file[:-4]}/pdb/*')
        for fn in name:
            fn = fn.split('/')[-1]
            df_fn = pd.read_csv(f'{path}/test/resources/{pdb_file[:-4]}/hbonds/{fn}_FARFAR-hbonds.csv')
            atom_1 = df_fn['Atom_1']
            atom_2 = df_fn['Atom_2']
            dist = df_fn['Distance']
            type = df_fn['Type']
            atom_h = df_fn['Hbond_atoms']
            index_fn = len(df_fn.index)

            for j in range(index_fn):
                if (type[j] == 'p'):
                    ppdb = PandasPdb()
                    ppdb.read_pdb(f'{path}/test/resources/{pdb_file[:-4]}/pdb/{fn}/{fn}.pdb')
                    ATOM = ppdb.df['ATOM']
                    atom = ATOM['atom_name']
                    resi_number = ATOM['residue_number']
                    length = len(ATOM.index)
                    pos_1 = atom_1[j].split('@')[1]
                    ps_1 = atom_1[j].split('@')[0]
                    pos_2 = atom_2[j].split('@')[1]
                    num = int(pos_2[1:len(pos_2)])
                    num_1 = int(pos_1[1:len(pos_1)])

                    ps_2 = atom_2[j].split('@')[0]
                    if (pos_1 == f'{n1}') or (pos_2 == f'{n1}'):
                        x_coord_1 = int(ATOM['x_coord'][atom == ps_1][resi_number == num_1])
                        y_coord_1 = int(ATOM['y_coord'][atom == ps_1][resi_number == num_1])
                        z_coord_1 = int(ATOM['z_coord'][atom == ps_1][resi_number == num_1])
                        x_coord_2 = int(ATOM['x_coord'][atom == ps_2][resi_number == num])
                        y_coord_2 = int(ATOM['y_coord'][atom == ps_2][resi_number == num])
                        z_coord_2 = int(ATOM['z_coord'][atom == ps_2][resi_number == num])
                        a = [x_coord_1, y_coord_1, z_coord_1]
                        b = [x_coord_2, y_coord_2, z_coord_2]
                        mod_a = np.sqrt((x_coord_1 ** 2) + (y_coord_1 ** 2) + (z_coord_1 ** 2))
                        mod_b = np.sqrt((x_coord_2 ** 2) + (y_coord_2 ** 2) + (z_coord_2 ** 2))
                        angle_radian = np.arccos(np.dot(a, b) / (mod_a * mod_b))
                        angle_degrees = math.degrees(angle_radian)

                        x = dist[j]
                        p1.append(x)
                        p2.append(fn)
                        p3.append(atom_1[j])
                        p4.append(atom_2[j])
                        p5.append(type[j])
                        p9.append(angle_degrees)
                        p10.append(pos_1)
                        p11.append(pos_2)
                        print(fn, dist[j])
                        if (ps_1 == f'{a1}') or (ps_2 == f'{a1}'):
                            p6.append('N-included')
                        else:
                            p6.append('Other')
                        p7.append(atom_h[j])
                        if (atom_h[j] == 'O:O'):
                            strength = (2.2 / dist[j]) * 21
                            p8.append(strength)
                        elif (atom_h[j] == 'N:N'):
                            strength = (2.2 / dist[j]) * 13
                            p8.append(strength)
                        elif (atom_h[j] == 'N:O'):
                            strength = (2.2 / dist[j]) * 8
                            p8.append(strength)
                        else:
                            strength = (2.2 / dist[j]) * 8
                            p8.append(strength)

                else:
                    continue

    df_h['Name'] = p2
    df_h['Atom_1'] = p3
    df_h['Nuc_1'] = p10
    df_h['Atom_2'] = p4
    df_h['Nuc_2'] = p11
    df_h['H-bond length'] = p1
    df_h['Type'] = p5
    df_h['N or other'] = p6
    df_h['Hbond_atoms'] = p7
    df_h['Hbond_strength'] = p8
    df_h['Angle'] = p9
    df_h.to_csv(f'{path}/test/resources/{pdb_file[:-4]}/hbond.csv', index=False)

def group_hbond(pdb_file,path):
    motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}').lower()
    df = pd.read_csv(f'{path}/test/resources/{pdb_file[:-4]}/hbond.csv')
    nuc_1 = df['Nuc_1']
    nuc_2 = df['Nuc_2']
    name = df['Name']
    hstrength = df['Hbond_strength']
    hlength = df['H-bond length']
    N = df['N or other']
    angle = df['Angle']
    index = df.index
    var_holder = {}
    a = motif.count('a')
    c = motif.count('c')
    j = a + c

    motif1 = motif.split('_')[0]
    motif2 = motif.split('_')[1]

    a_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'a']
    a_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'a']
    c_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'c']
    c_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'c']
    a_pos1_1 = []
    a_pos2_1 = []
    c_pos1_1 = []
    c_pos2_1 = []
    if a_pos1 != [] or a_pos2 != [] or c_pos1 != [] or c_pos2 != []:
        a_pos1_1 = multiplyList(a_pos1, 'A')
        a_pos2_1 = multiplyList(a_pos2, 'A')
        c_pos1_1 = multiplyList(c_pos1, 'C')
        c_pos2_1 = multiplyList(c_pos2, 'C')
    a_pos = a_pos1_1 + a_pos2_1
    c_pos = c_pos1_1 + c_pos2_1
    pos = a_pos + c_pos

    for n in range(j):
        var_holder['df_AC' + str(n)] = pd.DataFrame()
        var_holder['AC' + str(n) + '1'] = []
        var_holder['AC' + str(n) + '2'] = []
        var_holder['AC' + str(n) + '3'] = []
        var_holder['AC' + str(n) + '4'] = []
        var_holder['AC' + str(n) + '5'] = []
        var_holder['AC' + str(n) + '6'] = []
        var_holder['AC' + str(n) + '7'] = []

        df_AC = var_holder['df_AC' + str(n)]
        AC_1 = var_holder['AC' + str(n) + '1']
        AC_2 = var_holder['AC' + str(n) + '2']
        AC_3 = var_holder['AC' + str(n) + '3']
        AC_4 = var_holder['AC' + str(n) + '4']
        AC_5 = var_holder['AC' + str(n) + '5']
        AC_6 = var_holder['AC' + str(n) + '6']
        AC_7 = var_holder['AC' + str(n) + '7']

        nuc = pos[n]

        for i in range(len(index)):
            if (nuc_1[i] == f'{nuc}') or (nuc_2[i] == f'{nuc}'):
                AC_1.append(name[i])
                AC_2.append(nuc_1[i])
                AC_3.append(nuc_2[i])
                AC_4.append(hlength[i])
                AC_5.append(hstrength[i])
                AC_6.append(N[i])
                AC_7.append(angle[i])
        df_AC['Name'] = AC_1
        df_AC['Nuc_1'] = AC_2
        df_AC['Nuc_2'] = AC_3
        df_AC['H-bond length'] = AC_4
        df_AC['H-bond strength'] = AC_5
        df_AC['N or other'] = AC_6
        df_AC['Angle'] = AC_7

        df1_AC = df_AC.groupby(by=["Name", "N or other"]).agg({
            'H-bond length': 'mean',
            'H-bond strength': 'mean',
            'Angle': 'mean'
        })
        df2_AC = df_AC.groupby(by=["Name"]).agg({
            'H-bond length': 'mean',
            'H-bond strength': 'mean',
            'Angle': 'mean'
        })

        df1_AC.to_csv(f"{path}/test/resources/{pdb_file[:-4]}/N_{nuc}.csv")
        df2_AC.to_csv(f"{path}/test/resources/{pdb_file[:-4]}/hbond_mean_{nuc}.csv")

def H_bond(row, df):
    name_1 = row['name']
    for i in range(len(df.index)):
        name1 = df['Name'][i]
        print(name_1, name1)
        if name1 == name_1:
            result = df['H-bond strength'][i]
            return result

def add_hbond(pdb_file,path):
    motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}')
    motif1 = motif.split('_')[0]
    motif2 = motif.split('_')[1]

    a_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'A']
    a_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'A']

    c_pos1 = [(pos + 3) for pos, char in enumerate(motif1) if char == 'C']
    c_pos2 = [(pos + 7 + len(motif1)) for pos, char in enumerate(motif2) if char == 'C']

    pos1 = a_pos1 + a_pos2 + c_pos1 + c_pos2
    a_pos1_1 = []
    a_pos2_1 = []
    c_pos1_1 = []
    c_pos2_1 = []
    if a_pos1 != [] or a_pos2 != [] or c_pos1 != [] or c_pos2 != []:
        a_pos1_1 = multiplyList(a_pos1, 'A')
        a_pos2_1 = multiplyList(a_pos2, 'A')
        c_pos1_1 = multiplyList(c_pos1, 'C')
        c_pos2_1 = multiplyList(c_pos2, 'C')
    a_pos = a_pos1_1 + a_pos2_1
    c_pos = c_pos1_1 + c_pos2_1
    pos = a_pos + c_pos

    for l, k in zip(pos, pos1):
        sheet_num = k
        nt = l
        df = pd.read_csv(f"{path}/test/resources/{pdb_file[:-4]}/structural_parameters_{sheet_num}.csv")
        df_index = df.index

        df_N = pd.read_csv(f"{path}/test/resources/{pdb_file[:-4]}/N_{nt}.csv")
        df_mean = pd.read_csv(f"{path}/test/resources/{pdb_file[:-4]}/hbond_mean_{nt}.csv")

        df['H-bond strength_for_N'] = df.apply(lambda row: H_bond(row, df_N), axis=1)
        df['H-bond strength_mean'] = df.apply(lambda row: H_bond(row, df_mean), axis=1)

        df.fillna(0, inplace=True)
        df.to_csv(f"{path}/test/resources/{pdb_file[:-4]}/structural_parameters_{sheet_num}.csv", index=False)

def get_structures_less_than_RMSD(pdb_file,path):
    path1 = f'{path}/test/resources/{pdb_file[:-4]}/structural_parameters_*.csv'
    filename = sorted(glob.glob(path1))

    for n, fn in enumerate(filename):
        df_final = pd.read_csv(fn)
        df = df_final[df_final['RMSD'] < 2.5]
        df.to_csv(f'{path}/test/resources/{pdb_file[:-4]}/less_than_rmsd_{n + 1}.csv', index=False)
        print("DONE")

def get_average(pdb_file,path):
    motif = pdb_seq.sequence(f'{path}/test/resources/{pdb_file[:-4]}/{pdb_file}')
    path1 = f'{path}/test/resources/{pdb_file[:-4]}/less_than_rmsd_*.csv'
    filename = sorted(glob.glob(path1))
    Csv = pd.DataFrame()

    nt_dict = []
    alpha_dict = []
    beta_dict = []
    chi_dict = []
    gamma_dict = []
    delta_dict = []
    epsilon_dict = []
    zeta_dict = []
    e_z_dict = []
    ssZp_dict = []
    Dp_dict = []
    splay_dict = []
    eta_dict = []
    theta_dict = []
    eta_1_dict = []
    theta_1_dict = []
    eta_2_dict = []
    theta_2_dict = []
    v0_dict = []
    v1_dict = []
    v2_dict = []
    v3_dict = []
    v4_dict = []
    tm_dict = []
    P_dict = []
    suitness_dict = []
    SASA_2_dict = []
    nt_num_dict = []
    motif_dict = []
    H_bond_N_dict = []
    H_bond_all_dict = []

    for fn in filename:
        df = pd.read_csv(fn)
        motif1 = motif.split('_')[0]
        motif2 = motif.split('_')[1]
        nt = df['nt'].iloc[0]
        alpha = df['alpha'].mean(axis=0)
        beta = df['beta'].mean(axis=0)
        chi = df['chi'].str[1:7].astype(float).mean(axis=0)
        gamma = df['gamma'].mean(axis=0)
        delta = df['delta'].mean(axis=0)
        epsilon = df['epsilon'].mean(axis=0)
        zeta = df['zeta'].mean(axis=0)
        e_z = df['e-z'].str[:4].astype(float).mean(axis=0)
        ssZp = df['ssZp'].mean(axis=0)
        Dp = df['Dp'].mean(axis=0)
        splay = df['splay'].mean(axis=0)
        eta = df['eta'].mean(axis=0)
        theta = df['theta'].mean(axis=0)
        eta_1 = df['eta_1'].mean(axis=0)
        theta_1 = df['theta_1'].mean(axis=0)
        eta_2 = df['eta_2'].mean(axis=0)
        theta_2 = df['theta_2'].mean(axis=0)
        v0 = df['v0'].mean(axis=0)
        v1 = df['v1'].mean(axis=0)
        v2 = df['v2'].mean(axis=0)
        v3 = df['v3'].mean(axis=0)
        v4 = df['v4'].mean(axis=0)
        tm = df['tm'].mean(axis=0)
        P = df['P'].mean(axis=0)
        suitness = df['suitness'].mean(axis=0)
        SASA_2 = df['SASA_2'].mean(axis=0)
        nt_num = df['nt-num'].iloc[0]
        motif = df['motif'].iloc[0]
        H_bond_N = df['H-bond strength_for_N'].mean(axis=0)
        H_bond_all = df['H-bond strength_mean'].mean(axis=0)

        nt_dict.append(nt)
        alpha_dict.append(alpha)
        beta_dict.append(beta)
        chi_dict.append(chi)
        gamma_dict.append(gamma)
        delta_dict.append(delta)
        epsilon_dict.append(epsilon)
        zeta_dict.append(zeta)
        e_z_dict.append(e_z)
        ssZp_dict.append(ssZp)
        Dp_dict.append(Dp)
        splay_dict.append(splay)
        eta_dict.append(eta)
        theta_dict.append(theta)
        eta_1_dict.append(eta_1)
        theta_1_dict.append(theta_1)
        eta_2_dict.append(eta_2)
        theta_2_dict.append(theta_2)
        v0_dict.append(v0)
        v1_dict.append(v1)
        v2_dict.append(v2)
        v3_dict.append(v3)
        v4_dict.append(v4)
        tm_dict.append(tm)
        P_dict.append(P)
        suitness_dict.append(suitness)
        SASA_2_dict.append(SASA_2)
        nt_num_dict.append(nt_num)
        motif_dict.append(motif)
        H_bond_N_dict.append(H_bond_N)
        H_bond_all_dict.append(H_bond_all)

    Csv["motif"] = motif_dict
    Csv["nt"] = nt_dict
    Csv["nt_num"] = nt_num_dict
    Csv["alpha"] = alpha_dict
    Csv["beta"] = beta_dict
    Csv["chi"] = chi_dict
    Csv["gamma"] = gamma_dict
    Csv["delta"] = delta_dict
    Csv["epsilon"] = epsilon_dict
    Csv["zeta"] = zeta_dict
    Csv["e_z"] = e_z_dict
    Csv["ssZp"] = ssZp_dict
    Csv["Dp"] = Dp_dict
    Csv["splay"] = splay_dict
    Csv["eta"] = eta_dict
    Csv["theta"] = theta_dict
    Csv["eta_1"] = eta_1_dict
    Csv["theta_1"] = theta_1_dict
    Csv["eta_2"] = eta_2_dict
    Csv["theta_2"] = theta_2_dict
    Csv["v0"] = v0_dict
    Csv["v1"] = v1_dict
    Csv["v2"] = v2_dict
    Csv["v3"] = v3_dict
    Csv["v4"] = v4_dict
    Csv["tm"] = tm_dict
    Csv["P"] = P_dict
    Csv["suitness"] = suitness_dict
    Csv["SASA_2"] = SASA_2_dict
    Csv["H-bond strength_for_N"] = H_bond_N_dict
    Csv["H-bond strength_mean"] = H_bond_all_dict

    Csv.to_csv(f'{path}/test/resources/{pdb_file[:-4]}/average_{pdb_file[:-4]}.csv', index=False)

def concat_average_files(path):
    os.chdir(f'{path}/test/resources/')
    filename = glob.glob(f'{path}/test/resources/*/average_*.csv')
    frame = []
    N = []
    for fn in filename:
        df = pd.read_csv(fn)
        frame.append(df)

    df_final = pd.concat(frame)
    df_final.to_csv(f'{path}/test/results/Averages_concat_file.csv',index=False)

def string_pos(df):
    string = df['motif']

    for st in string:
        st1 = st.split('_')[0]
        st2 = st.split('_')[1]

        position_start1 = 3 ## Because when modeling two GC pairs are added on each ends of the motif.
        position_end1 = position_start1 + len(st1) - 1

        position_start2 = position_end1 + 5
        position_end2 = position_start2 + len(st2) - 1

        return position_start1, position_end1, position_start2, position_end2

def DMS(row, df):
    motif = row['motif']
    num_test = row['nt_num']
    pos = 0
    if row['start1'] <= num_test <= row['end1']:
        pos = num_test - row['start1']
    elif row['start2'] <= num_test <= row['end2']:
        pos = num_test - row['start2'] + (row['end1'] - row['start1'] + 1)
    else:
        pass
    for i in range(len(df.index)):
        if df['sequence'][i] == motif:
            print(df['sequence'][i], motif)
            result = df['norm_avg'][i + pos]
            return result

def safe_log(x):
    if x is None or x <= 0:
        return None  # Return None or handle differently depending on your requirement
    return log(x)

def add_DMS(path):
    df = pd.read_csv(f'{path}/test/results/Averages_concat_file.csv')
    res = string_pos(df)

    x1 = res[0]
    x2 = res[1]
    x3 = res[2]
    x4 = res[3]

    df['start1'] = x1
    df['end1'] = x2
    df['start2'] = x3
    df['end2'] = x4

    df_react = pd.read_csv(f'{path}/bin/Reactivity_values.csv')

    motif_test = df['motif']
    motif_react = df_react['sequence']

    df['DMS'] = df.apply(lambda row: DMS(row, df_react), axis=1)
    df['ln(DMS)'] = df['DMS'].apply(safe_log)

    df['ln(DMS)'].fillna(value=np.nan, inplace=True)
    
    df.to_csv(f'{path}/test/results/Average_concat_file_with_DMS.csv', index=False)

def get_keras_model(num_hidden_layers, num_neurons_per_layer, activation,input_layers):

    inputs = tf.keras.Input(shape=(input_layers.shape[1]))
    x = keras.layers.Dropout(0)(inputs)

    for i in range(num_hidden_layers):
        x = keras.layers.Dense(num_neurons_per_layer, activation=activation)(x)
        x = keras.layers.Dropout(0)(x)

    outputs = keras.layers.Dense(1, activation="linear")(x)

    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    return model

#### I changed all the paths to pdb/{model_pdb[:-4]}/{model_pdb}
def report(y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = stats.pearsonr(y_true, y_pred)[0] ** 2
    r2_value = '{:.3f}'.format(r2.item())
    return rmse, r2_value

def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b

def main():
    path = '/Users/sandunideenalattha/Desktop/RNA_twoway_junction_prediction'
    os.mkdir(f'{path}/test/results')
    os.chdir(f'{path}/test/resources/')
    fn = glob.glob('*.pdb') ## pdb file names
    print(fn)

    for f in fn:
        print(f)
        prep_fasta_secstruct(f,path)
        m = pdb_seq.sequence(f).lower()
        renumber_pdb(f,m)
        farfar(m,f)
        os.chdir(f'{path}/test/resources{f[:-4]}/pdb/')
        run_dssr(f,path)
        cal_hbond(f,path)
        convert_txt_to_csv(f,path)
        rmsd(f,path)
        read_dssr_file(f,path)
        add_RMSD(f,path)
        sasa_final = SASA(f,path)
        hbond(f,path)
        group_hbond(f,path)
        add_hbond(f,path)
        get_structures_less_than_RMSD(f,path)
        get_average(f,path)
    concat_average_files(path)
    add_DMS(path)

    test_df = pd.read_csv(f'{path}/test/results/Average_concat_file_with_DMS.csv')

    ind_var_test = test_df.drop(columns=['motif','nt','nt_num','start1','start2','end1','end2','DMS','ln(DMS)'])

    dep_var_test = test_df["ln(DMS)"]

    test_layers = []

    for i in ind_var_test.T:
        test_layers.append(ind_var_test.T[i])

    test_layers = np.array(test_layers)

    early_stopping = EarlyStopping(monitor='val_loss', mode='max', patience=10, restore_best_weights=True)
    keras_model = get_keras_model(1,217,'sigmoid',test_layers)
    keras_model.load_weights(f'{path}/bin/keras_model_regression_weights.h5')

    test_pred = keras_model.predict(test_layers)
    df_pred = pd.DataFrame()
    df_pred['test_pred'] = test_pred.ravel()
    df_pred['target_value'] = dep_var_test.values

    df_pred.to_csv(f'{path}/test/results/Predictions.csv', index=False)

    plt.figure(figsize=(8, 8))
    pred = df_pred['test_pred']
    target = df_pred['target_value']
    a, b = best_fit(target, pred)
    rmse, r2 = report(target, pred)
    target_np = target.to_numpy()

    plt.scatter(target_np, pred, color='green', label=f'r2 = {r2}')
    x = np.linspace(min(target_np), max(target_np), 100)
    yfit = [a + b * xi for xi in target_np]
    plt.plot(target_np, yfit, 'r-', lw=2, alpha=0.6, label=f'y = {yfit[0]} + {yfit[1]}x')
    plt.xticks(font='Arial', fontsize=18)
    plt.xlabel('Experimental ln(mutation fraction)', font='Arial', fontsize=20)
    plt.yticks(font='Arial', fontsize=18)
    plt.ylabel('Predicted ln(mutation fraction)', font='Arial', fontsize=20)
    plt.legend()
    plt.savefig(f'{path}/test/results/Regression_plot.png')

if __name__ == "__main__":
    main()