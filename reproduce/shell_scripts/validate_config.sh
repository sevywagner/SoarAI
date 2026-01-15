#!/usr/bin/bash
set -euo pipefail

validate_all=0

usage() { echo "usage: $0 <path to yaml config> [-a]"; }

while getopts ":a" opt; do
  case "$opt" in
    a) validate_all=1 ;;
    \?) echo "Unknown flag: -$OPTARG"; usage ; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

validate_dirs() {
	while IFS= read -r line; do
		if [ "${line: -1}" != "/" ]; then
			echo "${line} ** add a '/' at the end of all dir names **"
		fi
	done <<< "$1"
}

paths=`yq -r '
  .paths.core
  | to_entries[]
  | select(.key | test("_dir$"))
  | "\(.key)=\(.value)"
' "$1"`

validate_dirs "$paths"

if [ $validate_all -eq 1 ]; then
	paths=`yq -r '
  		.paths.api
  		| to_entries[]
  		| select(.key | test("_dir$"))
  		| "\(.key)=\(.value)"
	' "$1"`

	validate_dirs "$paths"

	paths=`yq -r '
  		.paths.analysis
  		| to_entries[]
  		| select(.key | test("_dir$"))
  		| "\(.key)=\(.value)"
	' "$1"`

	validate_dirs "$paths"
fi
