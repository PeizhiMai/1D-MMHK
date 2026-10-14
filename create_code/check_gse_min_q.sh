#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 L Q"
  exit 1
fi

L="$1"
target_q="$2"
half_fill=$((L / 2))
pattern="L${L}Nup${half_fill}Ndn${half_fill}f*_gse.dat"
found_any=0
all_ok=1

while IFS= read -r f; do
  found_any=1

  awk -v target_q="$target_q" '
    /^Q[[:space:]]/ {
      q = $2
    }
    /^E0[[:space:]]/ {
      e = $2 + 0
      if (!found || e < min) {
        min = e
        minq = q
        found = 1
      }
    }
    END {
      if (minq != target_q) {
        printf "%s  NOT_OK  min_E0=% .18e  Q=%s\n", FILENAME, min, minq
        exit 2
      }
    }
  ' "$f"
  awk_status=$?
  if [ "$awk_status" -eq 2 ]; then
    all_ok=0
  elif [ "$awk_status" -ne 0 ]; then
    exit "$awk_status"
  fi
done < <(find . -maxdepth 1 -type f -name "$pattern" -print | sed 's|^\./||' | sort)

if [ "$found_any" -eq 0 ]; then
  echo "No files matched $pattern"
  exit 1
fi

if [ "$all_ok" -eq 1 ]; then
  echo "All matching files for L=$L have minimum E0 at Q=$target_q"
fi
