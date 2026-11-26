# A vantagem de jogar em casa entre seleções nacionais

Código e dados do trabalho empírico de **Econometria I** (FGV EESP, 2º semestre de 2026),
de **Arthur Souza de Almeida Taveira** e **Bruno Lopes Peixoto**, sob orientação dos
professores Gabriela do Couto Fonseca e Leon Labre.

**Pergunta de pesquisa:** jogar em casa aumenta a probabilidade de vitória e o saldo de
gols em partidas internacionais de futebol entre seleções nacionais?

O trabalho compara 45.679 jogos entre seleções filiadas à FIFA (1881–2026) disputados no
país do mandante com jogos em campo neutro. O saldo de gols é estimado por mínimos
quadrados ordinários (MQO) e a vitória por probit, com controle de força por um rating Elo
calculado só com jogos anteriores e erros-padrão agrupados por mandante e visitante. As
extensões examinam a altitude, a evolução do efeito no tempo, os gols de pênalti e as
disputas de pênaltis.

## Conteúdo

```
.
├── Vantagem Mandante - Econometria I - Arthur e Bruno.R   script único com toda a análise
├── dados/
│   ├── brutos/      dados originais, como baixados em 18/11/2026
│   └── finais/      bases geradas pelo script
└── resultados/      tabelas, figuras e diagnósticos gerados pelo script
```

**`dados/brutos/`** (não editados)

| Arquivo | Conteúdo |
|---|---|
| `results.csv` | uma linha por jogo entre seleções (data, times, placar, torneio, cidade, país, campo neutro) |
| `goalscorers.csv` | uma linha por gol (autor, gol contra, pênalti) |
| `shootouts.csv` | uma linha por disputa de pênaltis (vencedor, quem cobrou primeiro) |
| `former_names.csv` | nomes antigos das seleções |
| `cities500.zip` | cidades do GeoNames com altitude (o script descompacta) |
| `countryInfo.txt` | países do GeoNames (código ISO, nome, capital) |
| `log_download.csv` | data e tamanho de cada arquivo baixado |

**`dados/finais/`**

| Arquivo | Conteúdo |
|---|---|
| `results_final.csv` | **base final da estimação** (amostra principal, 45.679 jogos), com todas as variáveis do modelo e os marcadores das extensões |
| `results_tratado.csv` | todos os jogos após o tratamento e a criação das variáveis, antes dos filtros de amostra |

Em CSV, os fatores viram texto: ao reler, refaça `torneio` (categoria de referência
"Amistoso"), `decada_jogo` e `periodo` como fatores.

**`resultados/`**: `descritivas/` (Tabelas 1 e 2, Figuras 1 e 2), `estimacao/`
(Tabelas 3 e 4), `diagnostico/` (testes), `robustez/` (Tabela 5), `extensoes/`
(Tabelas 6 e 7), além de `fluxo_amostra.csv`, `cidades_sem_altitude.csv` e
`sessao_R.txt` (versões do R e dos pacotes).

## Como reproduzir

1. Baixe ou clone o repositório.
2. Abra o R (ou o RStudio) com a **pasta do repositório como diretório de trabalho**.
   A linha `setwd()` no início do script está comentada; se preferir, use-a com o
   caminho da pasta no seu computador.
3. Instale os pacotes, uma vez:
   ```r
   install.packages(c("sandwich", "lmtest", "marginaleffects"))
   ```
4. Rode o script inteiro numa sessão limpa.

O script **não baixa os dados de novo** se eles já estiverem em `dados/brutos/`
(`atualizar_dados <- FALSE`). Isso fixa os resultados na versão de 18/11/2026: a base
original é atualizada com frequência, e um novo download muda os números. Nada no
script é aleatório.

Versões usadas: R 4.5.1, sandwich 3.1-3, lmtest 0.9-40, marginaleffects 1.0.0
(detalhes em `resultados/sessao_R.txt`).

## Fontes dos dados

- JÜRISOO, M. *International football results from 1872 to 2026*. Disponível em:
  <https://github.com/martj42/international_results>. Acesso em: 18 nov. 2026.
- GEONAMES. *GeoNames geographical database* (`cities500`, `countryInfo`). Licença
  CC BY 4.0. Disponível em: <https://www.geonames.org>. Acesso em: 18 nov. 2026.

## Uso de inteligência artificial

O uso de IA no trabalho, inclusive no código, está descrito na seção 7 do paper, que
traz o link público da conversa.
