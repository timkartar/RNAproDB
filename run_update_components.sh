#!/bin/bash
#cd /srv/www/rnascape/rnaview/
## update chemical components
rm components.cif
wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
python process_components.py
