#    ANARCI - Antibody Numbering and Antigen Receptor ClassIfication
#    Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the BSD 3-Clause License.

import argparse
import os
import sys

DESCRIPTION = '''

ANARCI                                                 \\\    //
Antibody Numbering and Antigen Receptor ClassIfication  \\\  //
                                                          ||
(c) Oxford Protein Informatics Group (OPIG). 2015-17      ||

Usage:

ANARCI -i <inputsequence or fasta file>

e.g.
    ANARCI -i Example_sequence_files/12e8.fasta
    This will number the files in 12e8.fasta with imgt numbering scheme and print to stdout.

    ANARCI -i Example_sequence_files/sequences.fasta -o Numbered_sequences.anarci -ht hit_tables.txt -s chothia -r ig
    This will number the files in sequences.fasta with chothia numbering scheme only if they are an antibody chain (ignore TCRs).
    It will put the numbered sequences in Numbered_sequences.anarci and the alignment statistics in hit_tables.txt

    ANARCI -i Example_sequence_files/lysozyme.fasta
    No antigen receptors should be found. The program will just list the names of the sequences.

    ANARCI -i EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA
    Or just give a single sequence to be numbered.
'''

EPILOGUE = '''
Author: James Dunbar (dunbar@stats.ox.ac.uk)
        Charlotte Deane (deane@stats.ox.ac.uk)

Contact: opig@stats.ox.ac.uk

Copyright (C) 2017 Oxford Protein Informatics Group (OPIG)
Freely distributed under the BSD 3-Clause Licence.
'''


def which(name, flags=os.X_OK):
    """Search PATH for executable files with the given name."""
    result = []
    exts = [_f for _f in os.environ.get('PATHEXT', '').split(os.pathsep) if _f]
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return list(set(result))


def main(argv=None):
    try:
        from anarci import scheme_names, wrap, all_species, run_anarci
    except ImportError as e:
        print("Fatal Error:", e, file=sys.stderr)
        sys.exit(1)

    parser = argparse.ArgumentParser(
        prog="ANARCI",
        description=DESCRIPTION,
        epilog=EPILOGUE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--sequence', '-i', type=str, help="A sequence or an input fasta file", dest="inputsequence")
    parser.add_argument('--outfile', '-o', type=str, default=False, help="The output file to use. Default is stdout", dest="outfile")
    parser.add_argument('--scheme', '-s', type=str, choices=scheme_names, default="imgt", help="Which numbering scheme should be used. i, k, c, m, w and a are shorthand for IMGT, Kabat, Chothia, Martin (Extended Chothia), Wolfguy and Aho respectively. Default IMGT", dest="scheme")
    parser.add_argument('--restrict', '-r', type=str, nargs="+", choices=["ig", "tr", "heavy", "light", "H", "K", "L", "A", "B"], default=False, help="Restrict ANARCI to only recognise certain types of receptor chains.", dest="restrict")
    parser.add_argument('--csv', action='store_true', default=False, help="Write the output in csv format. Outfile must be specified. A csv file is written for each chain type <outfile>_<chain_type>.csv. Kappa and lambda are considered together.", dest="csv")
    parser.add_argument('--outfile_hits', '-ht', type=str, default=False, help="Output file for domain hit tables for each sequence. Otherwise not output.", dest="hitfile")
    parser.add_argument('--hmmerpath', '-hp', type=str, default="", help="(Legacy) path to hmmer programs; ignored — this fork uses pyhmmer.", dest="hmmerpath")
    parser.add_argument('--ncpu', '-p', type=int, default=1, help="Number of parallel processes to use. Default is 1.", dest="ncpu")
    parser.add_argument('--assign_germline', action='store_true', default=False, help="Assign the v and j germlines to the sequence. The most sequence identical germline is assigned.", dest="assign_germline")
    parser.add_argument('--use_species', type=str, help="Use a specific species in the germline assignment. If not specified, only human and mouse germlines will be considered.", choices=all_species, dest="use_species")
    parser.add_argument('--bit_score_threshold', type=int, default=80, help="Change the bit score threshold used to confirm an alignment should be used.", dest="bit_score_threshold")

    args = parser.parse_args(argv)

    if argv is None and len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    outfile = False
    if args.outfile:
        path, fname = os.path.split(args.outfile)
        if (not path) or os.path.exists(path):
            outfile = args.outfile
        else:
            print("Error: Output file path does not exist", file=sys.stderr)
            sys.exit(1)

    if args.csv and not args.outfile:
        print("Error: When --csv option is used an ouput file name must be given.", file=sys.stderr)
        sys.exit(1)

    hitfile = False
    if args.hitfile:
        path, fname = os.path.split(args.hitfile)
        if (not path) or os.path.exists(path):
            hitfile = args.hitfile
        else:
            print("Error: Hit output file path does not exist", file=sys.stderr)
            sys.exit(1)

    if args.hmmerpath:
        print('Note: ANARCI is using pyhmmer, ignoring --hmmerpath argument', file=sys.stderr)

    types_to_chains = {"ig": ["H", "K", "L"], "tr": ["A", "B", "G", "D"], "heavy": ["H"], "light": ["K", "L"]}
    if args.restrict:
        allow = []
        for r in args.restrict:
            try:
                allow += types_to_chains[r]
            except KeyError:
                allow.append(r)
        allow = set(allow)
    else:
        allow = set(["H", "K", "L", "A", "B", "G", "D"])
    if args.scheme not in ("imgt", "i", "aho", "a") and allow - set(["H", "K", "L"]):
        print("Warning: Non IG chains cannot be numbered with the %s scheme. These will be ignored." % args.scheme, file=sys.stderr)
        allow = allow - set(["A", "B", "G", "D"])

    allowed_species = ['human', 'mouse']

    if args.use_species:
        assert args.use_species in all_species, 'Unknown species'
        allowed_species = [args.use_species]

    try:
        sequences, numbered, alignment_details, hit_tables = run_anarci(
            args.inputsequence,
            scheme=args.scheme,
            output=True,
            outfile=outfile,
            csv=args.csv,
            allow=allow,
            ncpu=args.ncpu,
            assign_germline=args.assign_germline,
            allowed_species=allowed_species,
            bit_score_threshold=args.bit_score_threshold,
        )

        if hitfile:
            with open(hitfile, "w") as hit_out:
                print("# Hit file for ANARCI", file=hit_out)
                for i in range(len(sequences)):
                    name, sequence = sequences[i]
                    print("NAME    ", name, file=hit_out)
                    print("SEQUENCE " + '\n'.join(['\nSEQUENCE '.join(wrap(block, width=71)) for block in sequence.splitlines()]), file=hit_out)

                    pad = max(list(map(len, hit_tables[i][0])))
                    for line in hit_tables[i]:
                        print(" ".join([str(_).rjust(pad) for _ in line]), file=hit_out)
                    print("//", file=hit_out)
    except Exception as e:
        print("Error: ", e, file=sys.stderr)
        sys.exit(1)

    sys.exit(0)
