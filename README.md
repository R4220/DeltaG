# DeltaG

Piccolo toolkit Python per workflow di termochimica atomistica basati su ASE.

Al momento il progetto copre tre casi d'uso principali:

- bulk armonico a cella fissa con phonons e termodinamica;
- termochimica di molecole isolate in approssimazione di gas ideale;
- post-processing di output `phonopy-qha`.

## Stato del progetto

Il codice e' gia' utilizzabile come prototipo locale, ma e' ancora in fase di consolidamento come pacchetto autosufficiente. La roadmap di sviluppo e' in:

- [Gibbs_energy/ROADMAP.md](./Gibbs_energy/ROADMAP.md)

## Struttura attuale

```text
DeltaG/
├── pyproject.toml
├── README.md
└── Gibbs_energy/
    ├── __init__.py
    ├── __main__.py
    ├── main.py
    ├── calculators.py
    ├── structures.py
    ├── phonons_fixed_cell.py
    ├── thermo_fixed_cell.py
    ├── molecule_thermo.py
    ├── qha_post.py
    ├── io_utils.py
    ├── constants.py
    └── ROADMAP.md
```

Nota: il nome del progetto installabile e del comando CLI e' `deltaG`, mentre il nome del pacchetto Python interno e' `Gibbs_energy`.

## Requisiti

- Python `>=3.10`
- `numpy`
- `matplotlib`
- `ase`

Per i workflow che usano il calcolatore MACE serve anche un ambiente con `mace-torch` disponibile, perche' [Gibbs_energy/calculators.py](./Gibbs_energy/calculators.py) usa `MACECalculator`.

Il workflow `qha-post` invece non richiede MACE.

## Installazione

Dalla directory radice del progetto:

```bash
cd /leonardo_work/IscrC_TDSSI/lorenzo/DeltaG
pip install -e .
```

Se il tuo ambiente non ha gia' MACE:

```bash
pip install mace-torch
```

Su cluster questa dipendenza potrebbe essere gestita separatamente via `conda` o ambiente modulare locale.

## Avvio rapido

Dopo l'installazione puoi usare:

```bash
deltaG --help
```

oppure, in alternativa:

```bash
python -m Gibbs_energy --help
```

## File di configurazione JSON/YAML

Per evitare command line molto lunghe, `deltaG` puo' leggere i parametri da un file di configurazione:

```bash
deltaG --config examples/fixed_cell.yaml
```

Sono supportati:

- file `.yaml` e `.yml`
- file `.json`

I file di configurazione devono usare la struttura annidata con:

- `mode: periodic`, `mode: molecule`, `mode: qha-post` oppure `mode: mixed`

Per ora il parser supporta:

- `mode: periodic`
- `mode: molecule`
- `mode: qha-post`

`mode: mixed` e' gia' definito nello schema, ma non ancora implementato nel parser.

Esempio YAML:

```yaml
mode: periodic

model:
  path: /path/to/mace.model
  device: cuda

thermo:
  temperature: 298.15
  pressure: 0.0

output:
  dir: fixed_cell_results

periodic:
  kind: bulk
  structure:
    geometry_file: structure.xyz
  phonons:
    supercell: [5, 5, 5]
    dos_kpts: [40, 40, 40]
```

Esempio JSON:

```json
{
  "mode": "qha-post",
  "output": {
    "dir": "qha_tables"
  },
  "qha_post": {
    "phonopy_qha": "phonopy-qha.out",
    "qha_summary": "qha_summary.out"
  }
}
```

Nota pratica:

- i path relativi nel file di config vengono interpretati rispetto alla cartella del config file
- il parser normalizza internamente il file verso i workflow CLI esistenti
- gli override da CLI restano volutamente semplici e usano chiavi flat interne, per esempio `temperature=500`

Se vuoi cambiare solo uno o due parametri senza modificare il file, puoi usare override leggeri da riga di comando:

```bash
deltaG --config examples/fixed_cell.yaml \
  --override temperature=500 \
  --override output_dir=run_500K
```

Gli override accettano anche valori JSON:

```bash
deltaG --config examples/fixed_cell.yaml \
  --override supercell=[4,4,4] \
  --override dos_kpts=[24,24,24]
```

Esempi pronti sono disponibili in:

- [examples/fixed_cell.yaml](./examples/fixed_cell.yaml)
- [examples/molecule.yaml](./examples/molecule.yaml)
- [examples/qha_post.yaml](./examples/qha_post.yaml)

## Workflow disponibili

### 1. `fixed-cell`

Calcola termodinamica armonica per un cristallo periodico a cella fissa:

- lettura della struttura o costruzione da reticolo + base;
- rilassamento geometrico;
- calcolo fononi con spostamenti finiti;
- DOS fononica;
- quantita' termodinamiche in funzione della temperatura.

Help rapido:

```bash
deltaG fixed-cell --help
```

Esempio minimale da file:

```bash
deltaG fixed-cell \
  --model-path /path/to/model.model \
  --geometry-file structure.xyz \
  --output-dir fixed_cell_results
```

Output tipici:

- `initial_crystal.xyz`
- `relaxed_crystal.xyz`
- `relax.log`
- `relax.traj`
- `phonon_BS_and_DOS.png`
- `fixed_cell_thermo_temperature.dat`
- `fixed_cell_summary.out`

### 2. `molecule`

Calcola termochimica di una molecola isolata con `IdealGasThermo`:

- costruzione della scatola di vuoto;
- rilassamento geometrico;
- vibrazioni molecolari;
- energia interna, entalpia, entropia, Helmholtz e Gibbs.

Help rapido:

```bash
deltaG molecule --help
```

Esempio minimale:

```bash
deltaG molecule \
  --model-path /path/to/model.model \
  --geometry-file molecule.xyz \
  --mol-geometry nonlinear \
  --symmetry-number 1 \
  --output-dir molecule_results
```

Output tipici:

- `initial_molecule.xyz`
- `relaxed_molecule.xyz`
- `relax_molecule.log`
- `relax_molecule.traj`
- `vibrations_summary.txt`
- `molecule_thermo_summary.out`

### 3. `qha-post`

Post-processa un file `phonopy-qha.out` e costruisce tabelle termodinamiche utili.

Help rapido:

```bash
deltaG qha-post --help
```

Esempio minimale:

```bash
deltaG qha-post \
  --phonopy-qha phonopy-qha.out \
  --output-dir qha_tables
```

Se disponibile, puoi passare anche un file riassuntivo aggiuntivo:

```bash
deltaG qha-post \
  --phonopy-qha phonopy-qha.out \
  --qha-summary qha_summary.out \
  --output-dir qha_tables
```

Output tipici:

- `enthalpy-temperature.dat`
- `entropy-temperature.dat`
- `gibbs-temperature.dat`
- `qha_thermo_temperature.dat`

## Filosofia del progetto

L'idea e' mantenere separati:

- il livello CLI in [Gibbs_energy/main.py](./Gibbs_energy/main.py);
- i workflow scientifici principali;
- i moduli di supporto per strutture, I/O e costanti.

Questo dovrebbe rendere piu' semplice:

- aggiungere nuovi workflow;
- migliorare validazione e test;
- trasformare il progetto in un pacchetto piu' robusto.

## Limitazioni attuali

- `fixed-cell` non include espansione termica: per quello serve la QHA.
- `molecule` e' appropriato per molecole isolate in fase gas, non per specie adsorbite.
- L'installazione di `mace-torch` puo' dipendere dall'ambiente di calcolo.
- Non ci sono ancora test automatici e casi esempio completi inclusi nel repository.

## Prossimi passi consigliati

- aggiungere un `README` piu' completo con esempi reali di input;
- inserire una cartella `examples/`;
- aggiungere test leggeri;
- rendere piu' robusta la validazione degli input;
- chiarire meglio la distinzione tra quantita' per cella e per formula unitaria.
