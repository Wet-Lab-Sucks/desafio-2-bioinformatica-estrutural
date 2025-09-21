import os
import pandas as pd
import glob

# three letters 
#AMINOACIDOS_CANONICOS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# single letter representation 
AMINOACIDOS_CANONICOS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def calcular_razao_hb_arquivo(arquivo_csv):
    """
    Calcula a razão entre o número de ligações de hidrogênio intra-cadeia
    e inter-cadeia para um único arquivo CSV, considerando apenas aminoácidos canônicos.

    Args:
        arquivo_csv (str): O caminho para o arquivo CSV.

    Returns:
        str: A razão formatada (float com uma casa decimal ou "0.0").
             Retorna None em caso de erro ao processar o arquivo.
    """
    total_hb_intra = 0
    total_hb_inter = 0

    try:
        df = pd.read_csv(arquivo_csv)
        for index, row in df.iterrows():
            tipo = row['Type']
            resname1 = row['ResName1'].upper()
            resname2 = row['ResName2'].upper()

            if tipo == 'HB' and resname1 in AMINOACIDOS_CANONICOS and resname2 in AMINOACIDOS_CANONICOS:
                chain1 = row['Chain1']
                chain2 = row['Chain2']
                res1 = row['Res1']
                res2 = row['Res2']

                # Desconsidera contatos entre átomos do mesmo resíduo
                if chain1 == chain2 and res1 != res2:
                    total_hb_intra += 1
                elif chain1 != chain2:
                    total_hb_inter += 1

        if total_hb_inter == 0:
            return "0.0"
        else:
            razao = total_hb_intra / total_hb_inter
            return f"{razao:.1f}"

    except FileNotFoundError:
        print(f"Arquivo não encontrado: {arquivo_csv}")
        return None
    except Exception as e:
        print(f"Erro ao processar o arquivo {arquivo_csv}: {e}")
        return None

def processar_lista_pdbs(arquivo_entrada="lista_pdbs.txt", arquivo_saida="razoes_aa_canonicos_hb.txt", diretorio_output="output"):
    """
    Lê uma lista de PDBs de um arquivo, processa os arquivos CSV correspondentes
    e salva as razões de ligações de hidrogênio (apenas entre aminoácidos canônicos)
    em um arquivo de saída.

    Args:
        arquivo_entrada (str): O caminho para o arquivo de texto com a lista de PDBs.
        arquivo_saida (str): O caminho para o arquivo de texto onde as razões serão salvas.
        diretorio_output (str): O diretório base contendo as pastas dos PDBs.
    """
    try:
        with open(arquivo_entrada, 'r') as f_entrada, open(arquivo_saida, 'w') as f_saida:
            for linha in f_entrada:
                pdb_id = linha.strip()
                # Encontra o arquivo CSV correspondente ao PDB ID
                caminho_csv = glob.glob(f"{diretorio_output}/{pdb_id}/*_contacts.csv")
                if caminho_csv:
                    # Assumindo que há apenas um arquivo CSV correspondente por PDB ID
                    razao = calcular_razao_hb_arquivo(caminho_csv[0])
                    if razao is not None:
                        f_saida.write(f"{razao}\n")
                    else:
                        f_saida.write("Erro ao calcular a razão\n") # Ou outra mensagem de erro
                else:
                    f_saida.write("Arquivo CSV não encontrado\n") # Ou outra mensagem indicando arquivo não encontrado
        print(f"Razões das ligações de hidrogênio (apenas canônicos) salvas em: {arquivo_saida}")

    except FileNotFoundError:
        print(f"Arquivo de entrada '{arquivo_entrada}' não encontrado.")
    except Exception as e:
        print(f"Ocorreu um erro: {e}")

if __name__ == "__main__":
    processar_lista_pdbs()
