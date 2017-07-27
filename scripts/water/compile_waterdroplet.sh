#!/bin/bash
d="$SCRATCH/droplet/gen_droplet"
g++ `root-config --glibs --cflags` $d/waterdroplet_tip4p_new.C -o $d/waterdroplet_tip4p_new.out
