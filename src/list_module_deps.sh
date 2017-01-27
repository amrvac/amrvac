#!/usr/bin/env bash

# Get all occurrences of use mod_... in .t files
deps="$(grep -r -e "use mod_" --include \*.t)"

# Remove the INCLUDES (are compiled first)
deps=$(echo "$deps" | sed 's/use mod_global_parameters//')
deps=$(echo "$deps" | sed 's/use mod_usr_methods//')
deps=$(echo "$deps" | sed 's/use mod_forest//')
deps=$(echo "$deps" | sed 's/use mod_physicaldata//')
deps=$(echo "$deps" | sed 's/use mod_connectivity//')

# Remove comments
deps=$(echo "$deps" | sed 's/!.*$//')

# Remove 'only: ...'
deps=$(echo "$deps" | sed 's/, *only.*$//')

# Remove lines without dependencies
deps=$(echo "$deps" | sed 's/^.*: *$//')

# Remove directories
deps=$(echo "$deps" | sed 's/^.*[/]//')

# Fix spacing around ':'
deps=$(echo "$deps" | sed 's/ *: */:/')

# Replace extension
deps=$(echo "$deps" | sed 's/[.]t/.o/')

# Replace 'use ' by ' '
deps=$(echo "$deps" | sed 's/use / /')

# Sort lines and remove duplicates
deps=$(echo "$deps" | sort -u)

# Print results
echo "$deps"
