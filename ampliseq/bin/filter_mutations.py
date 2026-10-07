#!/usr/bin/env python3
"""
Usage: python3 filter_mutations.py <input_tsv> <output_mutations>
"""

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_tsv> <linking_file>", file=sys.stderr)
        sys.exit(1)

    input_tsv, output_file = sys.argv[1], sys.argv[2]

    df      = pd.read_csv(input_tsv, sep="\t", dtype=str)

    

    out_file = os.path.join(os.getcwd(), os.path.basename(output_file))
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Written: {out_path}")

if __name__ == "__main__":
    main()
