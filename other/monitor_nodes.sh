#!/bin/bash

# Get unique nodes used by the user, excluding the header and empty lines
nodes=$(squeue -u kubacki.michal -h -o "%N" | tr -s " " "\n" | sort -u | grep -v "^$")

if [ -z "$nodes" ]; then
    echo "No active nodes found"
    exit 1
fi

# Create the command string for watch
cmd='nodes=$(squeue -u kubacki.michal -h -o "%N" | tr -s " " "\n" | sort -u | grep -v "^$"); '
cmd+='for node in $nodes; do '
cmd+='echo "Node: $node"; '
cmd+='ssh $node "top -bn1 | grep '"'"'%Cpu'"'"'" | sed "s/\x1B\[[0-9;]*[a-zA-Z]//g"; '
cmd+='echo; '
cmd+='done'

# Execute with watch
watch -n 2 "$cmd"