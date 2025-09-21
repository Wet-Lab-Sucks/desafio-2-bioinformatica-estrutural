# Desafio 2.3: Análise de Precisão Local de Modelagem de Proteínas (RMSD)

## Descrição do Desafio

A modelagem tridimensional de proteínas é um campo fundamental para a bioinformática estrutural.  

Neste desafio, fomos convidados a avaliar a precisão de um modelo predito da proteína PotD, comparando-o com a sua estrutura experimental ([PDB ID: 2V84](https://www.rcsb.org/structure/2V84)). 

Ao invés de usar uma métrica global como o RMSD ([Root Mean Square Deviation](https://en.wikipedia.org/wiki/Root_mean_square_deviation)) de toda a proteína, a tarefa exigia o cálculo de RMSD local para o resíduo de glicina na posição 55 (`Gly55`). 
>O cálculo deveria ser feito com base nas coordenadas de quatro átomos específicos (`nitrogênio (N)`, `carbono alfa (CA)`, `carbono (C)` e `oxigênio (O)`) do resíduo, fornecidos em dois arquivos `.pdb`: [PotD de referência](dataset/2v84.pdb) e [PotD modelada pela LBB](dataset/PotD_Modelada_LBB.pdb).

## Nossa estratégia

A nossa estratégia para a solução deste desafio foi a seguinte:

### 1. Extração de Coordenadas
* Desenvolvemos um script em Python ([calculate_RMSD.py](src/calculate_RMSD.py)) para ler os dois arquivos `.pdb` ([2v84.pdb](src/2v84.pdb) e [PotD_Modelada_LBB.pdb](src/PotD_Modelada_LBB.pdb)).

### 2. Identificação do Resíduo
* No script, implementamos a lógica para buscar especificamente o resíduo `Gly55` em ambos os arquivos.

### 3. Extração de Coordenadas Atômicas
* Extraímos as coordenadas `X`, `Y` e `Z` para os quatro átomos principais (`N`, `CA`, `C` e `O`) do `Gly55` em cada uma das estruturas.

### 4. Cálculo do RMSD
* Utilizamos a fórmula do RMSD fornecida no enunciado para calcular a raiz quadrada da média das distâncias quadráticas entre os átomos correspondentes das duas estruturas.

### 5. Formatação do Resultado
* O valor final foi formatado com duas casas decimais e salvo no arquivo [RMSD_residuo_Gly55.txt](results/RMSD_residuo_Gly55.txt).

 ## Resultado
 
 O resultado `31.09` para o RMSD do `Gly55` foi a nossa resposta final.  
 Tivemos a pontuação máxima neste desafio!
