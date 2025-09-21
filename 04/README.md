# Desafio 2.4 - Análise de Ligações de Hidrogênio com Critérios de Ângulo

## Descrição do Desafio

Continuando o tema de interações moleculares, este desafio envolveu a análise de interações moleculares, com foco na detecção de ligações de hidrogênio (`HB`).

A tarefa era ir além do critério de distância, que é comum em muitos softwares, e implementar uma análise mais acurada que considerasse também o ângulo formado entre o átomo doador, o átomo de hidrogênio e o átomo aceptor.

O objetivo era comparar a contagem de HBs em um conjunto de arquivos `.pdb` usando dois métodos:
* Apenas o critério de distância (similar ao desafio anterior).
* A combinação de critério de distância e angular (θ > 90º), conforme o modelo do software Mercury.
>[Macrae et al., 2008](https://onlinelibrary.wiley.com/doi/abs/10.1107/S0021889807067908) define que o contato ocorre apenas se os critérios de distância foram permissivos e se θ > 90º.

O desafio ressaltou a importância de parâmetros geométricos para uma avaliação precisa, especialmente em estruturas onde a posição do hidrogênio pode ser inferida. 

O objetivo final era retornar ambas as contagens, demonstrando a diferença entre as duas abordagens.

## Nossa estratégia

Para resolver este desafio, construímos um [script em python](src/hydrogen_bond_calculator.py) que implementou a lógica de detecção de `HBs`. A nossa estratégia foi dividida em três fases principais:

### Análise e Interpretação
O primeiro passo foi aprofundar nossa compreensão dos critérios geométricos para a formação de ligações de hidrogênio. Dedicamos tempo à leitura de referências como Macrae et al. (2008), que forneceu as diretrizes para a nossa abordagem angular.
> Macrae, C.F., et al. (2008) Mercury CSD 2.0-New Features for the Visualization and Investigation of Crystal Structures. Journal of Applied Crystallography, 41, 466-470. https://doi.org/10.1107/S0021889807067908 

### Desenvolvimento do Script
Implementamos um script que, para cada arquivo `.pdb`, identificava potenciais doadores (`N`, `O`, `F`) e aceptores. 

O script calculava a [distância euclidiana](https://pt.wikipedia.org/wiki/Dist%C3%A2ncia_euclidiana) entre eles para o primeiro critério. 

Para o segundo critério, o desafio era mais complexo: foi necessário calcular o ângulo entre as ligações `doador-hidrogênio` e `aceptor-hidrogênio`.

### Filtragem e Geração de Relatório
Por fim, o script aplicava os dois conjuntos de critérios para filtrar e contar as ligações de hidrogênio. 

O resultado de cada cálculo era formatado e salvo em um arquivo de saída, contendo as duas contagens separadas por uma vírgula, exatamente como solicitado:

```bash
cat results/result_2_4.txt
61,23
115,46
269,132
```

## Resultados

Enfrentamos incertezas e chegamos a resultados diferentes com base nos parâmetros que utilizamos. 

Mesmo assim, tivemos um bom resultado e também um grande aprendizado neste desafio.
