# Desafio 2.1 - Análise de Estruturas do PDB 

## Descrição do Desafio

Neste desafio de bioinformática estrutural, a tarefa era automatizar a análise de mais de 900 arquivos `.pdb`, a partir de uma lista de códigos identificadores do `Protein Data Bank (PDB)`. 

O objetivo era calcular a razão entre o número total de átomos e o número de resíduos para cada uma das estruturas. 
>Um PDB armazena as coordenadas 3D de átomos de uma proteína, e o cálculo exigia extrair essas informações. A razão final deveria ser formatada com duas casas decimais.

## Nossa estratégia

Nolssa abordagem para este desafio envolveu a manipulação de arquivos `.pdb`. 

Ao inspecionar os arquivos `.pdb`, identificamos que a contagem total de átomos para uma cadeia poderia ser facilmente obtida a partir da última linha que começa com a string `ATOM`. 

O número do átomo nessa linha (`ATOM 785...`) indicava o total.

Para o número de resíduos, era necessário somar os valores únicos da coluna de resíduos (a sexta coluna). A combinação dessas duas contagens nos fornecia os dados necessários para o cálculo, veja um exemplo simples:

```bash 
# identificar útima substring "ATOM"
tail 1a5k_B.pdb
ATOM    785  OXT LEU B 101      41.082  62.329 114.584  1.00 43.23           O
TER     796      HOH B 111
END

# resultado esperado é a razão entre átomos pelo número de resíduos com duas casas decimais 
 echo "scale=2; 785/101" | bc
7.77

# feito! só automatizar com um script!
```

Partindo desta interpretação, nosso script em Python foi desenvolvido para:

1. Ler a lista de códigos PDB de um arquivo de entrada (todos os 900 arquivos fornecidos pela LBB).
2. Para cada código, abrir e processar o arquivo `.pdb` correspondente.
3. Utilizar a identificação da última linha `"ATOM"` para extrair a contagem total de átomos.
4. Identificar e somar as ocorrências únicas de resíduos.
5. Calcular a razão `átomos / resíduos`.
6. Formatar o resultado com duas casas decimais e salvá-lo em um arquivo de saída `.txt`, conforme o formato solicitado.
