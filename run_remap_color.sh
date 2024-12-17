#!/bin/bash

for filename in /home/db/rnaprodb/output/old_electrostatics/*.ply; do
	python remap_ply_color.py $filename /home/db/rnaprodb/output/electrostatics/$(basename "$filename")
	break
done

