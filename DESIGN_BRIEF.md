# Marley -- Design Brief

> Documento de handoff para o designer da interface web.
> Ultima atualizacao: 2026-04-09

---

## 1. O que e o Marley?

Marley e uma **plataforma computacional de descoberta de vacinas e farmacos** contra a **leishmaniose visceral canina** (calazar) -- uma doenca tropical negligenciada que afeta caes e humanos.

O pipeline processa o proteoma inteiro do parasita (*Leishmania infantum*) atraves de 5 modulos sequenciais, desde a analise de proteinas ate a simulacao de resposta imune. O resultado final e um candidato a vacina de mRNA e moleculas candidatas a farmacos.

### Publico-alvo da interface

| Perfil | O que busca |
|--------|------------|
| **Pesquisadores/academia** | Dados detalhados, tabelas, metricas de validacao |
| **Banca de mestrado** | Visualizacao clara do pipeline, resultados reprodutiveis |
| **Investidores/pharma** | KPIs de impacto, comparacoes com mercado, potencial comercial |
| **Publico geral** | Entender o proposito do projeto de forma acessivel |

---

## 2. Stack tecnico

| Tecnologia | Versao | Papel |
|-----------|--------|-------|
| **Next.js** | 14 (App Router) | Framework React com SSR |
| **TypeScript** | Strict mode | Tipagem |
| **Tailwind CSS** | 3.4.4 | Estilizacao utility-first |
| **Supabase** | Configurado | Backend/DB (ainda nao usado ativamente) |

O designer **nao precisa programar** -- pode entregar wireframes, mockups ou prototipos em qualquer ferramenta (Figma, Sketch, Adobe XD, etc). A implementacao fica com o dev.

---

## 3. Sitemap atual

```
/                   Dashboard (home)
  |
  +-- /vaccine      Vaccine Construct (modulos v1-v2)
  |
  +-- /drug         Drug Targets (modulo v3)
  |
  +-- /docking      Molecular Docking (modulo v4)
  |
  +-- /simulation   Immune Simulation (modulo v5)
```

### Layout atual

```
+------------------+------------------------------------------+
|                  |                                          |
|    Sidebar       |          Conteudo principal              |
|    (w: 224px)    |          (scroll vertical)               |
|                  |                                          |
|  [M] Marley      |   Header + KPIs + Secoes                |
|                  |                                          |
|  H  Home         |                                          |
|  V  Vaccine      |                                          |
|  D  Drug Targets |                                          |
|  K  Docking      |                                          |
|  S  Simulation   |                                          |
|                  |                                          |
|  ──────────      |                                          |
|  Reverse Vacc.   |                                          |
|  Pipeline        |                                          |
+------------------+------------------------------------------+
```

---

## 4. Inventario de conteudo por pagina

### 4.1 Dashboard (Home `/`)

**Objetivo:** Visao geral do pipeline inteiro.

**Conteudo:**
- **6 KPI cards** em grid (2 colunas mobile, 6 desktop):
  - Proteins Analyzed: 8,527
  - Epitopes Selected: 11
  - Drug Targets: 52
  - Compounds Docked: 77
  - Custom Molecules: 20
  - Best Affinity: -7.74 kcal/mol
- **5 cards de pipeline** (links para modulos):
  - Cada card tem: tag do modulo (v1-v5), titulo, descricao curta
  - Cores diferentes por area tematica

### 4.2 Vaccine Construct (`/vaccine`)

**Objetivo:** Detalhes do design da vacina multi-epitopo.

**Conteudo:**
- **4 KPI cards**: Epitopes (11), Protein Length (379 aa), Mol. Weight (kDa), Instability Index
- **Tabela de epitopos** (11 linhas): Sequence, Source Gene, Allele, IC50, Rank
- **Barras comparativas** de antigenicidade: Marley vs Leish-Tec (referencia)
- **Imagem 3D** da estrutura da vacina (ESMFold prediction)

### 4.3 Drug Targets (`/drug`)

**Objetivo:** Alvos terapeuticos no parasita.

**Conteudo:**
- **4 KPI cards**: Total Targets (20), Priority Targets, Pathways, Best Score
- **Tabela top 20** alvos: Gene, Pathway, Identity, Druggability, Priority badge
- **Barras de distribuicao** por pathway metabolica
- **Card destaque** da molecula MRL-003 (custom, -7.74 kcal/mol)

### 4.4 Molecular Docking (`/docking`)

**Objetivo:** Resultados de triagem virtual de compostos.

**Conteudo:**
- **4 KPI cards**: Total Dockings, Unique Targets, Compounds, Best Affinity
- **Tabela top 10** docking hits: Target, Compound, Affinity, RMSD
- **Tabela de seletividade**: Parasite TryR vs Human GR (3 compostos)
- **2 imagens 3D**: GMPS + Methotrexate pose, TryR + MRL-003 pose

### 4.5 Immune Simulation (`/simulation`)

**Objetivo:** Modelagem da resposta imune a vacina.

**Conteudo:**
- **4 KPI cards**: Th1 Dominance (82.4%), Memory Stability, Est. Protection (>693 days), Dose Schedule (3 doses)
- **Grid de picos imunes** (5 tipos celulares): Th, Tc, B, Ab, Memory
- **Grafico de barras** de formacao de celulas de memoria (por dose)
- **2 graficos/imagens**: Cinetica imune (365 dias), Balanco Th1/Th2

---

## 5. Design system atual

### 5.1 Paleta de cores

| Token | Hex | Uso |
|-------|-----|-----|
| **Background** | `#0a0a0a` | Fundo principal |
| **Foreground** | `#ededed` | Texto principal |
| **Surface** | `zinc-900/50` | Cards, paineis |
| **Border** | `zinc-800` | Bordas, divisores |
| **Text muted** | `zinc-400` / `zinc-500` | Subtitulos, legendas |

### Cores por modulo

| Modulo | Cor primaria | Uso |
|--------|-------------|-----|
| Vaccine (v1-v2) | **Cyan** (`cyan-400` a `cyan-900`) | KPI borders, tags, dados |
| Drug Targets (v3) | **Orange** (`orange-400` a `orange-900`) | KPI borders, badges |
| Docking (v4) | **Green** (`green-400` a `green-900`) | KPI borders, dados |
| Simulation (v5) | **Purple** (`purple-400` a `purple-900`) | KPI borders, graficos |

### 5.2 Tipografia

- **Font family:** System UI (`system-ui, -apple-system, sans-serif`)
- **Dados numericos:** Monospace (font-mono do Tailwind)
- **Hierarquia:**
  - H1: `text-2xl` ou `text-3xl`, bold, white
  - H2 (section): `text-sm`, uppercase, tracking-wider, zinc-500
  - Body: `text-sm`, zinc-400
  - Data values: `text-2xl` bold white (KPIs), `text-xs` mono (tabelas)

### 5.3 Componentes existentes

#### KPI Card
```
+---+----------------------------+
| | |  PROTEINS ANALYZED          |  <- xs, uppercase, zinc-500
| | |  8,527                      |  <- 2xl, bold, white
| | |  L. infantum proteome       |  <- xs, zinc-500
+---+----------------------------+
 ^
 Borda esquerda colorida (4px, cor do modulo)
```

#### Sidebar Nav
```
+------+
| [M]  |  Marley         <- Logo (cyan-600, 32x32px)
|      |
| [H]  |  Home           <- Icon badge (zinc-800, 24x24)
| [V]  |  Vaccine        <- Hover: bg-zinc-800/60
| [D]  |  Drug Targets
| [K]  |  Docking
| [S]  |  Simulation
|      |
| ──── |
| Reverse Vaccinology     <- Footer info
| Pipeline
+------+
```

#### Tabela padrao
- Header: `bg-zinc-900/60`, uppercase, tracking-wider
- Rows: divider `zinc-800`, hover `zinc-900/40`
- Dados numericos: alinhados a direita, font-mono

#### Pipeline Card (Home)
- Border + background com cor do modulo (ex: `bg-cyan-900/30 border-cyan-800`)
- Tag do modulo (v1, v2...) em badge monospacado
- Titulo + descricao
- Hover: brightness-125

---

## 6. Assets existentes

### Imagens (em `/web/public/images/`)

| Arquivo | Conteudo | Pagina |
|---------|----------|--------|
| `vaccine_3d.png` | Estrutura 3D da vacina (ESMFold) | /vaccine |
| `docking_3d.png` | Pose de docking GMPS + Methotrexate | /docking |
| `tryr_mrl003.png` | Pose TryR + MRL-003 | /docking |
| `immune_kinetics.png` | Grafico cinetica imune 365 dias | /simulation |
| `th1_th2.png` | Grafico balanco Th1/Th2 | /simulation |

---

## 7. O que precisa de design

### 7.1 Melhorias no que ja existe

- [ ] **Responsividade mobile** -- sidebar colapsa? Drawer? Bottom nav?
- [ ] **Icones reais** -- atualmente sao letras (H, V, D, K, S). Precisam de icones
- [ ] **Hierarquia visual** -- melhorar escaneabilidade das tabelas grandes
- [ ] **Empty/loading states** -- nao existem ainda
- [ ] **Graficos interativos** -- atualmente sao imagens estaticas (PNG)
- [ ] **Feedback visual** -- nav ativa, breadcrumbs, estado atual do pipeline

### 7.2 Novas telas/funcionalidades planejadas

- [ ] **Landing page publica** -- apresentacao do projeto para quem nao e tecnico
- [ ] **Pagina ASO (Oligonucleotideo Antisense)** -- novo modulo de validacao matematica
- [ ] **Comparacao multi-organismo** -- pipeline roda para diferentes parasitas (Leishmania, T. cruzi)
- [ ] **Export/relatorio** -- gerar relatorio PDF dos resultados
- [ ] **Visualizador 3D interativo** -- rodar estruturas moleculares no browser (WebGL/Mol*)

### 7.3 Decisoes de design abertas

1. **Dark mode only ou light mode tambem?** -- atualmente so dark
2. **Idioma:** Portugues, ingles, ou bilingual?
3. **Branding:** Logo formal do Marley? Paleta definitiva?
4. **Tom visual:** Cientifico/serio? Moderno/startup? Acessivel/educativo?

---

## 8. Tipos de dados que o design precisa acomodar

| Tipo | Exemplo | Frequencia |
|------|---------|------------|
| **Numeros grandes** | 8,527 proteinas | Toda pagina (KPIs) |
| **Numeros decimais** | -7.74 kcal/mol | Tabelas de docking |
| **Porcentagens** | 82.4% Th1 dominance | KPIs, barras |
| **Sequencias biologicas** | `YAFGLYTSY` (9 chars) | Tabela de epitopos |
| **Nomes de genes** | GMPS, TryR, DHFR-TS | Tabelas |
| **Barras de progresso** | Score comparativo | Vaccine, drug pages |
| **Imagens 3D** | Estruturas moleculares (PNG) | Vaccine, docking |
| **Graficos temporais** | 365 dias de simulacao | Simulation |
| **Badges/tags** | Priority, v1-v5, Yes/No | Diversas |

---

## 9. Como rodar o projeto localmente

```bash
cd web
npm install
npm run dev
# Abre em http://localhost:3000
```

Requer Node.js 18+ e os arquivos de dados em `../results/`.

---

## 10. Glossario rapido

| Termo | O que e |
|-------|---------|
| **Epitopo** | Pedaco de proteina que o sistema imune reconhece |
| **Docking** | Simulacao computacional de como uma molecula encaixa num alvo |
| **KPI** | Key Performance Indicator -- metrica principal |
| **IC50** | Concentracao que inibe 50% da atividade (menor = melhor) |
| **kcal/mol** | Unidade de energia de ligacao (mais negativo = ligacao mais forte) |
| **mRNA** | RNA mensageiro -- a tecnologia da vacina (mesma do COVID) |
| **Th1/Th2** | Tipos de resposta imune (Th1 = celular, melhor contra parasitas) |
| **ASO** | Oligonucleotideo antisense -- terapia genica para silenciar genes |
| **TryR** | Tripanotiona redutase -- enzima essencial do parasita |
| **VaxiJen** | Score de antigenicidade (>0.4 = provavel antigeno) |
| **ESMFold** | IA que prediz estrutura 3D de proteinas (Meta AI) |
