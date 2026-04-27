# bulk_thermo roadmap

## Obiettivo
Trasformare questa cartella in un piccolo programma autosufficiente per termochimica con tre workflow:

- bulk armonico a cella fissa;
- termochimica di molecole isolate in approssimazione di gas ideale;
- post-processing QHA.

## Stato attuale
La struttura del codice e gia' abbastanza buona:

- `main.py` espone una CLI con tre sottocomandi;
- `thermo_fixed_cell.py` contiene il workflow bulk;
- `molecule_thermo.py` contiene il workflow molecolare;
- `qha_post.py` contiene il post-processing QHA;
- `structures.py`, `calculators.py`, `io_utils.py`, `constants.py` fanno da supporto.

In pratica, il progetto e' gia' un prototipo funzionante, ma non ancora un programma veramente autosufficiente e facile da distribuire.

## Priorita' consigliata

### 1. Igiene da pacchetto
- Rinominare `__init.py` in `__init__.py`.
- Passare dagli import locali assoluti a import relativi, per esempio `from .qha_post import ...`.
- Aggiungere un `pyproject.toml` minimale con nome pacchetto e dipendenze.
- Aggiungere un `README.md` con installazione, prerequisiti ed esempi.

### 2. Configurazione chiara degli input
- Decidere un solo stile di input principale:
  - CLI pura con molti flag, oppure
  - file di configurazione `yaml/json` con pochi flag da riga di comando.
- Definire uno schema chiaro per i parametri comuni:
  - `model_path`
  - `device`
  - `temperature` o griglia di temperature
  - `pressure`
  - directory di output
- Separare bene input per bulk, molecola e QHA.

### 3. Robustezza del workflow
- Validare in anticipo i file di input e stampare errori leggibili.
- Verificare che supercelle, mesh DOS e bandpath abbiano formati corretti.
- Aggiungere controlli su valori fisicamente sospetti:
  - frequenze immaginarie;
  - formula units nulle o inconsistenti;
  - geometrie molecolari incompatibili con `IdealGasThermo`.
- Salvare un file riassuntivo dei parametri usati in ogni run.

### 4. Esperienza d'uso
- Uniformare i nomi dei file prodotti tra i workflow.
- Scrivere help CLI piu' espliciti con esempi reali.
- Aggiungere una cartella `examples/` con:
  - un esempio bulk da file;
  - un esempio bulk da reticolo e base;
  - un esempio molecolare;
  - un esempio QHA post.

### 5. Validazione scientifica
- Confrontare almeno un caso bulk e un caso molecolare con riferimenti noti.
- Verificare le unita' in tutti gli output.
- Decidere chiaramente cosa viene riportato per cella e cosa per formula unitaria.
- Documentare i limiti del metodo:
  - harmonic approximation;
  - fixed-cell vs QHA;
  - validita' di `IdealGasThermo` solo per molecole isolate.

### 6. Test automatici
- Aggiungere test unitari per:
  - `infer_formula_units`
  - `formula_strings`
  - parsing di `phonopy-qha.out`
  - derivate finite in `qha_post.py`
- Aggiungere almeno un smoke test per la CLI senza eseguire calcoli costosi.

## Ordine pratico di lavoro
1. Sistemare packaging e import.
2. Scrivere README e uno o due esempi minimi.
3. Rendere piu' robusti input e messaggi di errore.
4. Aggiungere test leggeri.
5. Solo dopo, ampliare le funzionalita' scientifiche.

## Prima lista di miglioramenti concreti
- Creare `pyproject.toml`.
- Rinominare `__init.py`.
- Convertire gli import a forma package-safe.
- Aggiungere `README.md`.
- Aggiungere una cartella `examples/`.
- Aggiungere una cartella `tests/`.

## Nota importante
La base e' gia' sufficientemente modulare. Il rischio principale, prima ancora dei dettagli scientifici, e' che oggi il codice sembri un insieme di script locali. Il passaggio chiave per renderlo autosufficiente e' farlo diventare un piccolo pacchetto con input, output e dipendenze ben definiti.
