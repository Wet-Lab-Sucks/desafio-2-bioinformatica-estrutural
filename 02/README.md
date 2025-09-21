# Desafio 2.2 - Análise de Ligações de Hidrogênio Intra e Inter-Cadeia

## Descrição do Desafio

Neste desafio, exploramos abordagens de bioinformática estrutural, focando na importância das ligações de hidrogênio (`HB`) para a estabilidade e função de proteínas. 

A tarefa era analisar uma lista de complexos `proteína-peptídeo` para calcular a razão entre as `ligações de hidrogênio intra-cadeia` (dentro da mesma cadeia proteica) e as `ligações inter-cadeia` (entre cadeias diferentes), utilizando critérios de distância como o do software [COCaDA](https://github.com/rplemos/COCaDA).

>O desafio nos permitiu demonstrar nossa proficiência na manipulação de arquivos `.pdb` e na aplicação de ferramentas de análise molecular, destacando como interações aparentemente "fracas" são, na verdade, essenciais para a estrutura tridimensional e a função biológica.

## Nossa estratégia

A nossa abordagem para este problema seguiu a seguinte metodologia:

### Análise Estrutural e Filtragem
Para cada complexo, [nosso script](src/parse_cocada_distances_canonical_aa.py) utilizou o [COCaDA](https://github.com/rplemos/COCaDA) para processar e identificar todas as ligações de hidrogênio. As interações foram então filtradas em dois grupos: `intra-cadeia` (entre resíduos da mesma cadeia) e `inter-cadeia` (entre resíduos de cadeias diferentes).

### Cálculo da Razão
Com as contagens separadas, nosso script calculou a razão entre as `ligações intra-cadeia` e `inter-cadeia`.

### Geração do Relatório Final 
O resultado final foi um arquivo de texto com as razões calculadas para cada complexo, pronto para a submissão:

```bash
head results/LBB_2_2_saida.txt
7.2
32.5
29.0
4.8
28.0
5.9
29.0
30.5
30.0
3.6
```
