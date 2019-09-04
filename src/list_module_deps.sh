#!/usr/bin/env bash

# Make Linux sort order match mac OS X
export LC_ALL=C
# Get all occurrences of use mod_... in .t files
# - The '.' is for Mac compatibility
# - The 'sort' is to ensure the order is identical on different systems
deps="$(grep -r -e "^\s*use mod_" --include \*.t . | sort)"

# Remove the INCLUDES (are compiled first)
deps=$(echo "$deps" | sed 's/use mod_global_parameters//')
deps=$(echo "$deps" | sed 's/use mod_usr_methods//')
deps=$(echo "$deps" | sed 's/use mod_forest//')
deps=$(echo "$deps" | sed 's/use mod_physicaldata//')
deps=$(echo "$deps" | sed 's/use mod_connectivity//')
deps=$(echo "$deps" | sed 's/use mod_constants//')
deps=$(echo "$deps" | sed 's/use mod_variables//')

# Remove old_physics entries
deps=$(echo "$deps" | sed 's/^.*old_physics[/].*$//')

# Remove amrvac.o entry (it is compiled later)
deps=$(echo "$deps" | sed 's/^amrvac[.]t.*$//')

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

# Replace 'use mod_xxx' by ' mod_xxx.mod'
deps=$(echo "$deps" | sed 's/use \(.*\)$/ \1.mod/')

# Sort lines and remove duplicates
deps=$(echo "$deps" | sort -u)

# Print results
echo "$deps"
