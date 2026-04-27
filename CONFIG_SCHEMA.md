# DeltaG Input Schema Proposal

Questa nota propone una struttura piu' chiara e scalabile per i file di input di `DeltaG`.

L'idea e' separare:

- il tipo di workflow;
- i dati strutturali;
- i parametri termodinamici;
- le opzioni di plotting e output.

Questa proposta e' ora anche il formato richiesto per i file di configurazione.

Stato attuale del parser:

- `mode: periodic` supportato
- `mode: patch` supportato
- `mode: molecule` supportato
- `mode: qha-post` supportato
- `mode: mixed` non ancora implementato

## Obiettivi

Lo schema dovrebbe:

- essere leggibile a mano;
- evitare command line troppo lunghe;
- distinguere bene tra input comuni e input specifici del workflow;
- gestire in modo naturale sistemi `bulk`, `surface`, `patch`, `molecule` e `mixed`;
- rendere facile aggiungere nuovi workflow in futuro.

## Scelta consigliata

Per file scritti a mano, la scelta migliore e' `YAML`.

Motivi:

- e' piu' leggibile di JSON;
- gestisce bene strutture annidate;
- e' comodo per liste, blocchi `plots`, componenti multiple, ecc.

JSON puo' comunque restare supportato come alternativa tecnica.

## Struttura top-level proposta

```yaml
mode: periodic | patch | molecule | mixed | qha-post

model:
  path: /path/to/model.model
  device: cuda

thermo:
  temperature: 298.15
  pressure: 0.0

plots:
  enabled: true
  phonon_dos: true
  band_structure: false
  thermo_curves: true
  format: png
  dpi: 200

output:
  dir: results
  write_summary: true
```

Per il workflow `periodic`, puoi anche aggiungere un blocco opzionale `qha`
per lanciare in fondo il post-processing di un `phonopy-qha.out` gia'
esistente.

## Significato dei blocchi comuni

### `mode`

Propongo quattro modalita':

- `periodic`
- `patch`
- `molecule`
- `mixed`

In piu', per il post-processing QHA il parser supporta anche:

- `qha-post`

### `model`

Contiene il calcolatore elettronico/ML:

```yaml
model:
  path: /path/to/model.model
  device: cuda
```

### `thermo`

Contiene temperatura, pressione e in futuro eventuali griglie termiche:

```yaml
thermo:
  temperature: 298.15
  pressure: 0.0
```

Nota:

- per sistemi periodici la pressione puo' essere espressa in `GPa`
- per molecole la pressione puo' essere espressa in `Pa`

Questa differenza andra' poi chiarita bene nell'implementazione o resa esplicita con un campo unita'.

### `plots`

Secondo me e' meglio tenerlo strutturato nel config, non come molte flag CLI.

```yaml
plots:
  enabled: true
  phonon_dos: true
  band_structure: true
  thermo_curves: true
  format: png
  dpi: 200
```

### `output`

Centralizza tutte le scelte di scrittura dei risultati:

```yaml
output:
  dir: results
  write_summary: true
```

In futuro qui si possono aggiungere:

- `save_intermediate_structures`
- `save_trajectories`
- `overwrite`

## Modalita' 1: `periodic`

Questa modalita' copre sia `bulk` sia `surface`.

Per evitare una flag ambigua come `Bulk/surface`, propongo:

```yaml
mode: periodic

periodic:
  kind: bulk | surface
```

### Blocchi consigliati

```yaml
mode: periodic

periodic:
  kind: bulk
  structure:
    geometry_file: bulk.xyz

  formula_units: 4

  relax:
    fmax: 0.01

  phonons:
    supercell: [5, 5, 5]
    delta: 0.05
    dos_kpts: [40, 40, 40]
    bandpath: GMKG
    emax: 0.035

  temperature_grid:
    t_min: 0.0
    t_max: 1000.0
    t_step: 50.0

qha:
  enabled: false
  phonopy_qha: phonopy-qha.out
  output_dir: qha_tables
```

### Nota su `structure`

Per strutture periodiche conviene supportare due alternative:

```yaml
structure:
  geometry_file: structure.xyz
```

oppure:

```yaml
structure:
  lattice_geometry:
    a: 3.0
    b: 3.0
    c: 4.5
    alpha: 90
    beta: 90
    gamma: 120
    Inum: 194
  lattice_basis:
    - symbol: Mo
      position: [0.0, 0.0, 0.0]
    - symbol: S
      position: [0.333333, 0.666667, 0.25]
```

## Modalita' 2: `patch`

Questa modalita' e' per frammenti finiti in vuoto trattati con
`HarmonicThermo`, senza contributi di gas ideale.

### Blocchi consigliati

```yaml
mode: patch

patch:
  geometry_file: patch.xyz
  vacuum: 15.0

  relax:
    fmax: 0.01

  constraints:
    fixed_indices: [0, 1, 2, 3]

  vibrations:
    indices: [4, 5, 6, 7]
    delta: 0.01
    clean: false

  temperature_grid:
    t_min: 0.0
    t_max: 1000.0
    t_step: 50.0
```

### Nota concettuale

Qui non va usato `IdealGasThermo`: il workflow aggiunge solo contributi
vibrazionali armonici. Se il frammento rappresenta una patch supportata,
conviene fissare alcuni atomi di bordo oppure limitare gli atomi vibrati.

## Modalita' 3: `molecule`

Questa modalita' e' per molecole isolate e termochimica tipo `IdealGasThermo`.

### Blocchi consigliati

```yaml
mode: molecule

molecule:
  geometry_file: molecule.xyz
  mol_geometry: nonlinear
  symmetry_number: 1
  spin: 0.0
  vacuum: 15.0

  relax:
    fmax: 0.01

  vibrations:
    clean: false
```

### Nota concettuale

Qui non userei `bulk/surface/molecule` come singola categoria, perche' la fisica e il workflow molecolare sono davvero diversi.

## Modalita' 4: `mixed`

Questa e' la parte piu' importante della tua idea.

Se `mixed` significa casi tipo:

- `surface + molecule`
- `surface + adsorbate`
- `surface + molecule + adsorbate`

allora non conviene far derivare tutto da una singola flag implicita. E' meglio dichiarare esplicitamente i componenti presenti.

### Struttura consigliata

```yaml
mode: mixed

mixed:
  kind: adsorption

  components:
    clean_surface:
      type: surface
      geometry_file: surf.xyz

    gas_molecule:
      type: molecule
      geometry_file: mol.xyz
      mol_geometry: nonlinear
      symmetry_number: 1
      spin: 0.0
      vacuum: 15.0

    adsorbed_state:
      type: surface_adsorbate
      geometry_file: surf_ads.xyz

  references:
    use_clean_surface: true
    use_gas_molecule: true

  analysis:
    adsorption_energy: true
    adsorption_free_energy: true
```

### Perche' cosi' funziona meglio

Con questa struttura:

- il file dice esattamente quali stati hai;
- puoi aggiungere nuovi stati senza cambiare il concetto di base;
- il codice puo' controllare meglio coerenza e campi mancanti;
- diventa naturale estendere il workflow a reazioni piu' complesse.

## Esempio completo 1: periodic bulk

```yaml
mode: periodic

model:
  path: /path/to/model.model
  device: cuda

thermo:
  temperature: 298.15
  pressure: 0.0

plots:
  enabled: true
  phonon_dos: true
  band_structure: true
  thermo_curves: true
  format: png
  dpi: 200

output:
  dir: bulk_results
  write_summary: true

qha:
  enabled: false
  phonopy_qha: phonopy-qha.out
  output_dir: qha_tables

periodic:
  kind: bulk
  structure:
    geometry_file: bulk.xyz

  formula_units: 4

  relax:
    fmax: 0.01

  phonons:
    supercell: [5, 5, 5]
    delta: 0.05
    dos_kpts: [40, 40, 40]
    bandpath: GMKG
    emax: 0.035

  temperature_grid:
    t_min: 0.0
    t_max: 1000.0
    t_step: 50.0
```

## Esempio completo 2: molecule

```yaml
mode: molecule

model:
  path: /path/to/model.model
  device: cuda

thermo:
  temperature: 298.15
  pressure: 101325.0

plots:
  enabled: false

output:
  dir: molecule_results
  write_summary: true

molecule:
  geometry_file: molecule.xyz
  mol_geometry: nonlinear
  symmetry_number: 1
  spin: 0.0
  vacuum: 15.0

  relax:
    fmax: 0.01

  vibrations:
    clean: false
```

## Esempio completo 3: mixed adsorption

```yaml
mode: mixed

model:
  path: /path/to/model.model
  device: cuda

thermo:
  temperature: 298.15
  pressure: 101325.0

plots:
  enabled: true
  thermo_curves: true
  format: png
  dpi: 200

output:
  dir: adsorption_results
  write_summary: true

mixed:
  kind: adsorption

  components:
    clean_surface:
      type: surface
      geometry_file: surf.xyz

    gas_molecule:
      type: molecule
      geometry_file: mol.xyz
      mol_geometry: nonlinear
      symmetry_number: 1
      spin: 0.0
      vacuum: 15.0

    adsorbed_state:
      type: surface_adsorbate
      geometry_file: surf_plus_ads.xyz

  relax:
    fmax: 0.01

  analysis:
    adsorption_energy: true
    adsorption_free_energy: true
```

## Consiglio architetturale

Io implementerei questa evoluzione in due passi.

### Passo 1

Mantenere i config file semplici gia' esistenti e aggiungere solo un nuovo formato annidato documentato.

### Passo 2

Introdurre una funzione di normalizzazione che trasformi il config annidato in una struttura interna standard, per esempio:

- `workflow`
- `model`
- `structures`
- `calculation`
- `plots`
- `output`

Cosi' i moduli scientifici lavorano sempre su uno schema interno coerente.

## Scelta che ti consiglierei oggi

Se devo scegliere la direzione giusta per `DeltaG`, farei cosi':

- `YAML` come formato umano principale
- `JSON` come alternativa supportata
- `mode: periodic | molecule | mixed`
- `periodic.kind: bulk | surface`
- `mixed.components` esplicito
- `plots` come blocco nel config, non come tante flag CLI

## Prossimo passo pratico

Il prossimo passo migliore non e' ancora implementare tutto, ma fissare:

1. i nomi definitivi dei blocchi top-level;
2. le unita' di `pressure` nei diversi mode;
3. quali componenti di `mixed` vuoi considerare obbligatori.

Dopo questo, si puo' tradurre lo schema in:

- esempi YAML realistici;
- validazione dei config;
- normalizzazione interna verso i workflow attuali.
