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

N_LIST="16 32 64"

# Density loop.
initial_density=0.1
final_density=2.0
density_step=0.1
DENSITY_LIST="$(seq $initial_density $density_step $final_density)"

# Physical / algorithmic parameters
FREEZE_MC_STEPS_SCALE=0.5
ALPHA=1.0
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
for DENSITY in $DENSITY_LIST; do
    for N in $N_LIST; do
        RUN_NAME="n${N}_density${DENSITY}_t${T_INITIAL}-${T_FINAL}-${T_STEP}"
        OUTDIR="${BASE_OUTDIR}/${RUN_NAME}"

        mkdir -p "$OUTDIR"

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
EOF

        # Run simulation
        ./wigner_crystal | tee -a "$OUTDIR/run.log"

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
