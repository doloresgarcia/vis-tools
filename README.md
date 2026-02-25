## Usage

### Environment

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

### MC table

```bash
python3 ./dump_mc_table.py /path/to/file.edm4hep.root evt_number
```

### MC tree

```bash
python3 ./dump_mc_tree.py /path/to/file.edm4hep.root evt_number
```


### Event Display

```bash
glced &
python3 ./event_display.py $K4GEO/path/to/detector/compact/file.xml /path/to/file.edm4hep.root
```
