
## metaparse.py

Usage:

      metaparse.py cmd station


*cmd*

 :    command to be executed:

- `c` or `compact`: compact print of metadata set
- `e` or `extended`: extended output of metadat set
- `p` or `plain`: plain dump of python structures
- `x` or `xml`: XML output of metadata set
- `n` or `normcheck`: formal check of normalization and poles&zeros


*station*

 :    Station specifier like `GR.BFO`. For Station code `*` as wild card is allowed, e.g. `GR.*`, but quotes are needed on the command line.


Examples:

      metaparse.py c ge.eil
      metaparse.py n 'gr.*'
