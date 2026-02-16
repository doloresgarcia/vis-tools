# Dump MC Tree Tool

## Requirements

Key4hep software stack

```bash
# stable release
source /cvmfs/sw.hsf.org/key4hep/setup.sh

# or nighlies
# source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

## Usage

```bash
./dump-mc-tree.py filename.edm4hep.root evt_number
```

## Output

SVG MC record tree `./doctest-output/Event{evt_number}.gv.svg`

![MC record of the Z->ss event](image.png)