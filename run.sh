#! /bin/bash


workspace=$1
hostname=$2
nproc=$3
ntask_per_node=$5
nodes=$4
branch=$(git ls-remote --heads origin | grep $(git rev-parse HEAD) | cut -d / -f 3)

echo "Directory: $workspace"
echo "Machine: $hostname"
echo "Branch name: $branch"

# Get the branch of pdata

echo "RUN numerics test"

cd ..
source $HOME/openfpm_vars_$branch

echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
pwd
echo "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB "

mpirun -np $3 ../build/openfpm_numerics/src/numerics
if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_numerics test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
fi


