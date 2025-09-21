#!/usr/bin/env python

import os
import sys
import math
from Bio.PDB import PDBParser, Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np
import warnings

# suprimir warning do Biopython
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# Critérios para ligações de hidrogênio baseados no exercício anterior (COCαDA, Lemos et al., 2024)
MAX_DISTANCE_DA = 3.5  # Distância máxima entre doador e aceptor (D_ad ≤ 3.9 Å)
MAX_DISTANCE_HA = 2.5  # Distância máxima entre hidrogênio e aceptor (mantida)
MIN_ANGLE_DHA = 90.0   # Ângulo mínimo D-H-A conforme Mercury (Macrae et al., 2008)

# Mapeamento de átomos doadores e aceptores com base na tabela fornecida
# A tabela identifica quais átomos em cada resíduo são doadores e aceptores

# Lista de átomos doadores por resíduo
DONORS_BY_RESIDUE = {
    "ALA": [],
    "ARG": ["NE", "NH1", "NH2"],
    "ASN": ["ND2"],
    "ASP": [],
    "CYS": ["SG"],
    "GLN": ["NE2"],
    "GLU": [],
    "GLY": [],
    "HIS": ["ND1", "NE2"],
    "ILE": [],
    "LEU": [],
    "LYS": ["NZ"],
    "MET": [],
    "PHE": [],
    "PRO": [],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TRP": ["NE1"],
    "TYR": ["OH"],
    "VAL": [],
    # Qualquer resíduo (incluindo não padrão e backbone)
    "ANY": ["N"]
}

# Lista de átomos aceptores por resíduo
ACCEPTORS_BY_RESIDUE = {
    "ALA": [],
    "ARG": [],
    "ASN": ["OD1"],
    "ASP": ["OD1", "OD2"],
    "CYS": [],
    "GLN": ["OE1"],
    "GLU": ["OE1", "OE2"],
    "GLY": [],
    "HIS": ["ND1", "NE2"],
    "ILE": [],
    "LEU": [],
    "LYS": [],
    "MET": ["SD"],
    "PHE": [],
    "PRO": ["O"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TRP": [],
    "TYR": ["OH"],
    "VAL": [],
    # Qualquer resíduo (incluindo não padrão e backbone)
    "ANY": ["O"]
}

def is_donor(atom):
    """
    Verifica se um átomo é um doador de hidrogênio com base na sua identidade
    e no resíduo ao qual pertence.
    """
    residue = atom.parent
    residue_name = residue.resname
    atom_name = atom.name
    
    # Verifica se o átomo está na lista de doadores para seu resíduo específico
    if residue_name in DONORS_BY_RESIDUE and atom_name in DONORS_BY_RESIDUE[residue_name]:
        return True
    
    # Verifica se o átomo está na lista geral de doadores (backbone)
    if atom_name in DONORS_BY_RESIDUE["ANY"]:
        return True
    
    # Átomos especiais de nucleotídeos
    if atom_name in ["N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]:
        if residue_name in ["A", "G", "C", "T", "U", "DA", "DG", "DC", "DT"]:
            return True
    
    return False

def is_acceptor(atom):
    """
    Verifica se um átomo é um aceptor de hidrogênio com base na sua identidade
    e no resíduo ao qual pertence.
    """
    residue = atom.parent
    residue_name = residue.resname
    atom_name = atom.name
    
    # Verifica se o átomo está na lista de aceptores para seu resíduo específico
    if residue_name in ACCEPTORS_BY_RESIDUE and atom_name in ACCEPTORS_BY_RESIDUE[residue_name]:
        return True
    
    # Verifica se o átomo está na lista geral de aceptores (backbone)
    if atom_name in ACCEPTORS_BY_RESIDUE["ANY"]:
        return True
    
    # Átomos especiais de nucleotídeos
    if atom_name in ["O2", "O4", "O6", "OP1", "OP2", "N1", "N3", "N7"]:
        if residue_name in ["A", "G", "C", "T", "U", "DA", "DG", "DC", "DT"]:
            return True
    
    return False

def calculate_distance(atom1, atom2):
    """Calcula a distância euclidiana entre dois átomos."""
    coord1 = np.array(atom1.coord)
    coord2 = np.array(atom2.coord)
    return np.linalg.norm(coord1 - coord2)

def calculate_angle(donor, hydrogen, acceptor):
    """
    Calcula o ângulo θ (em graus) do ângulo D-H-A (doador-hidrogênio-aceptor).
    O ângulo é calculado usando os vetores D-H e A-H que se encontram no hidrogênio.
    """
    # Vetores de coordenadas
    d_coord = np.array(donor.coord)
    h_coord = np.array(hydrogen.coord)
    a_coord = np.array(acceptor.coord)
    
    # Vetores D-H e A-H (ambos direcionados para H)
    dh_vector = h_coord - d_coord
    ah_vector = h_coord - a_coord
    
    # Normalizar os vetores
    dh_norm = np.linalg.norm(dh_vector)
    ah_norm = np.linalg.norm(ah_vector)
    
    if dh_norm == 0 or ah_norm == 0:
        return 0.0
    
    # Normalizar os vetores para obter os vetores unitários
    dh_unit = dh_vector / dh_norm
    ah_unit = ah_vector / ah_norm
    
    # Calcular o produto escalar
    dot_product = np.dot(dh_unit, ah_unit)
    
    # Garantir que o valor esteja dentro do intervalo válido [-1, 1]
    dot_product = max(-1.0, min(1.0, dot_product))
    
    # Calcular o ângulo em radianos e converter para graus
    angle_rad = math.acos(dot_product)
    angle_deg = math.degrees(angle_rad)
    
    return angle_deg

def find_hydrogens_for_donor(donor, structure):
    """Encontra átomos de hidrogênio ligados a um átomo doador."""
    hydrogens = []
    
    # Considerar todos os átomos de hidrogênio
    for atom in structure.get_atoms():
        if atom.element == "H":
            distance = calculate_distance(donor, atom)
            # Distância típica de ligação covalente D-H (um pouco mais permissiva)
            if distance < 1.6:  
                hydrogens.append(atom)
    
    # Se não encontrou hidrogênios, podemos estar lidando com um arquivo PDB
    # que não tem explicitamente os hidrogênios. Neste caso, adicionamos um
    # hidrogênio virtual com coordenadas calculadas para possibilitar os cálculos.
    """
    if not hydrogens and donor.name == "N":  # Backbone N
        try:
            # Tenta encontrar o CA e C da mesma residue para calcular a posição do H
            residue = donor.parent
            if "CA" in residue and "C" in residue:
                # Encontrar o C da residue anterior
                prev_residue = None
                chain = residue.parent
                residues = list(chain)
                for i, res in enumerate(residues):
                    if res == residue and i > 0:
                        prev_residue = residues[i-1]
                        break
                
                if prev_residue and "C" in prev_residue:
                    # Calcular a posição do H baseada na geometria conhecida
                    n_coord = np.array(donor.coord)
                    ca_coord = np.array(residue["CA"].coord)
                    c_prev_coord = np.array(prev_residue["C"].coord)
                    
                    # Calcular o vetor normalizado da direção N-H
                    v1 = n_coord - ca_coord
                    v1 = v1 / np.linalg.norm(v1)
                    v2 = n_coord - c_prev_coord
                    v2 = v2 / np.linalg.norm(v2)
                    
                    # Direção N-H é aproximadamente a soma dos vetores normalizados
                    h_direction = v1 + v2
                    h_direction = h_direction / np.linalg.norm(h_direction)
                    
                    # Distância típica N-H: ~1.0 Å
                    h_coord = n_coord + h_direction * 1.0
                    
                    # Criar um átomo virtual H
                    from Bio.PDB.Atom import Atom
                    virtual_h = Atom("H", h_coord, 0.0, 1.0, " ", " H ", 0, "H")
                    hydrogens.append(virtual_h)
        except Exception as e:
            # Em caso de erro, simplesmente não adiciona hidrogênio virtual
            pass

    """
    return hydrogens

def analyze_pdb(pdb_id, pdb_dir="./"):
    """
    Analisa um arquivo PDB e conta as ligações de hidrogênio usando critérios 
    de distância e ângulo.
    """
    pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    
    if not os.path.exists(pdb_path):
        print(f"Erro: Arquivo PDB {pdb_path} não encontrado.")
        return 0, 0
    
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, pdb_path)
    except Exception as e:
        print(f"Erro ao analisar o arquivo PDB {pdb_id}: {e}")
        return 0, 0
    
    hbonds_distance_only = 0
    hbonds_with_angle = 0
    
    # Considerar todas as possíveis ligações entre as cadeias e modelos
    for model in structure:
        # Identificar todos os potenciais doadores e aceptores na estrutura
        all_donors = []
        all_acceptors = []
        
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if is_donor(atom):
                        all_donors.append(atom)
                    if is_acceptor(atom):
                        all_acceptors.append(atom)
        
        # Para cada doador, verificar possíveis ligações com cada aceptor
        for donor in all_donors:
            # Encontrar átomos de hidrogênio ligados ao doador
            hydrogens = find_hydrogens_for_donor(donor, structure)
            
            for acceptor in all_acceptors:
                # Ignorar ligações dentro do mesmo resíduo
                if donor.parent is acceptor.parent:
                    continue
                
                # Calcular distância doador-aceptor
                da_distance = calculate_distance(donor, acceptor)
                
                # Verificar critério de distância COCαDA (D_ad ≤ 3.9 Å)
                if da_distance <= MAX_DISTANCE_DA:
                    # Contabilizar uma ligação baseada apenas em distância
                    hbonds_distance_only += 1
                    
                    # Se não houver hidrogênios explícitos, não podemos verificar critério angular
                    if not hydrogens:
                        continue
                        
                    # Para cada átomo de hidrogênio associado com este doador
                    for hydrogen in hydrogens:
                        # Calcular o ângulo D-H-A
                        angle = calculate_angle(donor, hydrogen, acceptor)
                        
                        # Verificar o critério de ângulo (θ > 90º conforme Mercury)
                        if angle > MIN_ANGLE_DHA:
                            hbonds_with_angle += 1
                            # Uma vez que encontramos uma ligação satisfatória com este par D-A,
                            # não precisamos verificar outros hidrogênios
                            break
                                
    return hbonds_distance_only, hbonds_with_angle

def process_pdb_list(output_file, pdb_dir="./"):
    """
    Processa uma lista de IDs PDB e grava os resultados 
    em um arquivo de saída.
    """
    
    pdb_ids = ["1k0p", "5jxv", "7quu"]
    results = []
    for pdb_id in pdb_ids:
        print(f"Processando PDB: {pdb_id}")
        hbonds_distance, hbonds_angle = analyze_pdb(pdb_id, pdb_dir)
        results.append(f"{hbonds_distance},{hbonds_angle}")
    
    with open(output_file, 'w') as f:
        for result in results:
            f.write(f"{result}\n")
    
    print(f"Resultados salvos em {output_file}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python hydrogen_bond_calculator.py output.txt [pdb_dir]")
        sys.exit(1)
    
    output_file = sys.argv[1]
    pdb_dir = sys.argv[2] if len(sys.argv) > 2 else "./dataset/"
    
    # Verificar se o diretório de PDBs existe
    if not os.path.exists(pdb_dir):
        print(f"Erro: Diretório de PDBs '{pdb_dir}' não encontrado.")
        sys.exit(1)
        
    try:
        process_pdb_list(output_file, pdb_dir)
    except Exception as e:
        print(f"Erro ao processar a lista de PDBs: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
