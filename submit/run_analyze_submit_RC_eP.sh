#!/bin/sh

#Add full patch to your eic-shell here:
/gpfs/mnt/gpfs02/eic/janvanek/eic/eic-shell << EOF

root -b -q './submit_workdir/analyze_eD_DIS_eicrecon_RC_eP.C+("${1}","${2}")'

exit
EOF
