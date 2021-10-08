#!/bin/bash
#SBATCH --job-name=ICE_DYNAMICS_TEST
#SBATCH --ntasks=32
#SBATCH --time=12:00:00
#SBATCH --partition=mix

mpirun /data90t/geosci/spetrov/DYNAMICS_TEST_NEW/build/DYNAMICS_TESTS ./config.json > ./output.txt
