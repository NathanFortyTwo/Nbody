#!/bin/bash

reachable=""
unreachable=""

# Set the SSH timeout in seconds (adjust as needed)
ssh_timeout=2

for i in {01..13}; do
    machine="3a401-$i"
    echo "Testing:" $machine
    if ssh -q -o BatchMode=yes -o ConnectTimeout=$ssh_timeout $machine exit; then
        # Use parameter expansion to remove the leading space
        reachable="${reachable} $machine"
    else
        unreachable="${unreachable} $machine"
    fi
done

# Remove the leading space from the reachable machines list
reachable="${reachable# }"

echo "Reachable machines:$reachable"
echo "Unreachable machines:$unreachable"

# Generate a hosts file with reachable machines
if [ -n "$reachable" ]; then
    echo "Generating hosts file with reachable machines..."
    echo "$reachable" | tr ' ' '\n' >hosts
    echo "Hosts file 'hosts' created with reachable machines."
else
    echo "No reachable machines found."
fi
