## RNAproDB backend setup
1. Ensure Python 3 is installed
2. Run: `git checkout apache-server`
3. `pip install -e pydpb/`
4. Install other Python dependencies as required
5. Ensure x3dna-dssr is available in PATH
6. Other external dependencies are available in ./external/ and can be added to the path
7. Run: `python rna_vis.py 1ivs` (you may need to adjust pdb_path and output paths in rna_vis.py and run_dssr.py as needed)

Example output is available at [https://rohslab.usc.edu/rnaprodb/1ivs](https://rohslab.usc.edu/rnaprodb/1ivs)
