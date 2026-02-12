#!/bin/sh
#
# run_mc_coulomb.sh
# Introduction to Scientific Computing
# Marta Xiulan Aribó Herrera - (c) GNU GPLv3 2025
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#

set -e

# -----------------------------------------------------------------------------
# Sweep parameters
# -----------------------------------------------------------------------------

N_LIST="64 125"

# SC n³: 64, 125, 216, 343...
# BCC 2n³: 54, 128, 250, 434...
# FCC 4n³: 108, 256, 500, 864...

# Density loop.
initial_density=0.1
final_density=2.0
density_step=0.1
DENSITY_LIST=$(LC_NUMERIC=C awk "BEGIN{for(d=$initial_density; d<=$final_density+1e-9; d+=$density_step) printf \"%.2f \", d}")

# Physical / algorithmic parameters
FREEZE_MC_STEPS_SCALE=0.5
ALPHA=0.5
CHARGE=1.0

T_INITIAL=1.0
T_FINAL=0.0001
T_STEP=0.98

# Output base directory
BASE_OUTDIR=output
mkdir -p "$BASE_OUTDIR"

#
# Loop over N and Density
#
for N in $N_LIST; do
    for DENSITY in $DENSITY_LIST; do
        RUN_NAME="n${N}_density${DENSITY}_t${T_INITIAL}-${T_FINAL}-${T_STEP}"
        OUTDIR="${BASE_OUTDIR}/${RUN_NAME}"

        mkdir -p "$OUTDIR"

        L=$(LC_NUMERIC=C awk "BEGIN { print ($N / $DENSITY)^(1/3) }")

        cat << _header_
-------------------------------------------------------------------------------
Running N      = $N
Density        = $DENSITY
Temperatures   = $T_INITIAL - $T_FINAL (Step: $T_STEP)
Output dir     = $OUTDIR
-------------------------------------------------------------------------------
_header_

        # Generate input file
        cat > input_parameters.in << EOF
$N                      ! Number of particles
$FREEZE_MC_STEPS_SCALE  ! Freeze MC steps
$ALPHA                  ! Parameter alpha
$DENSITY                ! Density
$CHARGE                 ! Particle charge
$T_INITIAL              ! Initial temperature
$T_FINAL                ! Final temperature
$T_STEP                 ! Temperature step
$L                      ! Box size
EOF

        # Run simulation
        ./mc_annealing_wigner_crystal | tee -a "$OUTDIR/run.log"

        #
        # Collect outputs
        #
        mv input_parameters.in "$OUTDIR/"

        for f in *.out *.xyz
        do
            [ -f "$f" ] && mv "$f" "$OUTDIR/"
        done
    done
done

printf "\\n%s\\n" 'All sweeps completed successfully.'
