"""Back-propagate PSI-MS CV accessions into the column database.

Fills the `psi_ms_id` column of resources/column_database/column_database.tsv from the
PSI-MS id ledger (.github/repo-rt/psi-ms-column-ids.tsv in HUPO-PSI/psi-ms-CV), which is
the only authority for which (company, column) model has which MS:5xxxxxx accession.
PSI-MS reads this column back on its weekly sync: a model whose name was corrected keeps
its psi_ms_id here, and that is what lets PSI-MS treat the change as a rename (id kept,
old name kept as a synonym) instead of retiring the old term and minting a new one.

Rules (they mirror apply_mirror() in PSI-MS' generate_psi_ms_columns.py):
  - key in the ledger                -> write the ledger id on every row of the model,
                                        overwriting whatever was there
  - key not in the ledger, rows carry the id of a model that has LEFT the catalog
                                     -> pending rename: keep the id (and fill it on any
                                        blank row of that model)
  - key not in the ledger, any other id (a live model's id, e.g. rows copied from an
    existing model; a vendor id; an id the ledger never minted)
                                     -> blank it, otherwise the PSI-MS sync aborts
  - key not in the ledger, no id     -> leave blank (new model, minted on the next sync)
  - blank company or column          -> row left untouched (PSI-MS drops these rows)
The join key is normalised exactly as PSI-MS does (clean()); the company/column cells
themselves are never rewritten.

Usage:
    python scripts/Python/add_psi_ms_id.py [--ledger URL_OR_PATH] [--catalog TSV] [--check]
"""

import argparse
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

# TODO: switch to the HUPO-PSI master URL once HUPO-PSI/psi-ms-CV#544 is merged:
# https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/.github/repo-rt/psi-ms-column-ids.tsv
LEDGER_DEFAULT = ('https://raw.githubusercontent.com/jonathan-hunter/psi-ms-CV/'
                  'repo-rt-column-generator-2/.github/repo-rt/psi-ms-column-ids.tsv')
CATALOG_DEFAULT = 'resources/column_database/column_database.tsv'
LEDGER_HEADER = 'company\tcolumn\tpsi_ms_id'
ID_COLUMN = 'psi_ms_id'
MS_ID = re.compile(r'MS:\d{7}')

# copied verbatim from PSI-MS scripts/chrom_columns/generate_psi_ms_columns.py: the ledger
# keys are catalog cells after this normalisation, so the join must apply it too
INVISIBLE = str.maketrans(
    {
        " ": " ",  # no-break space
        "​": "",  # zero-width space
        "‌": "",  # zero-width non-joiner
        "‍": "",  # zero-width joiner
        "⁠": "",  # word joiner
        "﻿": "",  # zero-width no-break space / BOM
    }
)


def clean(text):
    """Normalize a catalog cell for use as an OBO name/def value and as an identity key:
    drop control characters, fold invisible/non-breaking characters, collapse whitespace."""
    text = "".join(ch for ch in text if ch >= " " and ch != "\x7f")
    return re.sub(r"\s+", " ", text.translate(INVISIBLE)).strip()


def read_ledger(source):
    """{(company, column): psi_ms_id} from a local path or URL, validated like PSI-MS does."""
    if re.match(r'https?://', source):
        with urllib.request.urlopen(source, timeout=60) as response:
            text = response.read().decode('utf-8')
    else:
        text = Path(source).read_text(encoding='utf-8')
    lines = text.split('\n')
    if lines[0] != LEDGER_HEADER:
        sys.exit(f'{source}: ledger header must be {LEDGER_HEADER!r}, got {lines[0]!r}')
    ledger, seen_ids = {}, set()
    for n, line in enumerate(lines[1:], start=2):
        if line == '':
            continue
        fields = line.split('\t')
        if len(fields) != 3 or not MS_ID.fullmatch(fields[2]):
            sys.exit(f'{source} line {n}: expected company<TAB>column<TAB>MS:nnnnnnn, got {line!r}')
        company, column, ms_id = fields
        if (company, column) in ledger or ms_id in seen_ids:
            sys.exit(f'{source} line {n}: duplicate key or id: {line!r}')
        ledger[(company, column)] = ms_id
        seen_ids.add(ms_id)
    return ledger


def read_catalog(path):
    """(header, rows) with every cell kept verbatim; no quote processing (as PSI-MS reads it)."""
    try:
        text = Path(path).read_text(encoding='utf-8')
    except UnicodeDecodeError as e:
        sys.exit(f'{path} is not valid UTF-8: {e}')
    lines = text.split('\n')
    if lines[-1] == '':
        lines.pop()
    header = lines[0].split('\t')
    rows = []
    for n, line in enumerate(lines[1:], start=2):
        fields = line.split('\t')
        if len(fields) > len(header):
            sys.exit(f'{path} line {n}: {len(fields)} fields but the header has {len(header)}')
        rows.append(fields + [''] * (len(header) - len(fields)))
    return header, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--ledger', default=LEDGER_DEFAULT, help='PSI-MS id ledger (URL or path)')
    parser.add_argument('--catalog', default=CATALOG_DEFAULT, help='column database TSV to update in place')
    parser.add_argument('--check', action='store_true',
                        help='write nothing; exit 1 if the catalog would change')
    args = parser.parse_args()

    ledger = read_ledger(args.ledger)
    holder = {ms_id: key for key, ms_id in ledger.items()}
    header, rows = read_catalog(args.catalog)
    if ID_COLUMN not in header:
        header.append(ID_COLUMN)
        for row in rows:
            row.append('')
    company_i, column_i, id_i = (header.index(c) for c in ('company', 'column', ID_COLUMN))

    models = defaultdict(list)  # normalised key -> row indices
    for i, row in enumerate(rows):
        key = (clean(row[company_i].strip()), clean(row[column_i].strip()))
        if key[0] and key[1]:
            models[key].append(i)

    new_ids = {}  # row index -> psi_ms_id to write
    filled, kept, blanked, conflicts = [], [], [], []
    for key, indices in models.items():
        if key in ledger:
            filled.append(key)
            new_ids.update((i, ledger[key]) for i in indices)
            continue
        ids = {rows[i][id_i].strip() for i in indices} - {''}
        if len(ids) > 1:
            conflicts.append((key, sorted(ids)))
            continue
        old_key = holder.get(next(iter(ids))) if ids else None
        if old_key is not None and old_key[1] != '' and old_key not in models:
            kept.append((old_key, key, ledger[old_key]))
            new_ids.update((i, ledger[old_key]) for i in indices)
        else:
            if ids:
                blanked.append((key, next(iter(ids))))
            new_ids.update((i, '') for i in indices)

    # two new keys claiming the same retired id: PSI-MS can only rename one of them
    claims = defaultdict(list)
    for old_key, key, ms_id in kept:
        claims[ms_id].append(key)
    conflicts += [(key, [ms_id]) for ms_id, keys in claims.items() if len(keys) > 1 for key in keys]
    if conflicts:
        print('psi_ms_id conflicts; fix these rows by hand, nothing was written:', file=sys.stderr)
        for key, ids in conflicts:
            print(f'  {key!r}: {", ".join(ids)}', file=sys.stderr)
        sys.exit(2)

    changed = 0
    for i, ms_id in new_ids.items():
        if rows[i][id_i] != ms_id:
            rows[i][id_i] = ms_id
            changed += 1
    output = '\n'.join('\t'.join(fields) for fields in [header, *rows]) + '\n'
    original = Path(args.catalog).read_text(encoding='utf-8')

    retired = [key for key in ledger if key[1] != '' and key not in models]
    print(f'models filled from the ledger: {len(filled)} ({sum(len(models[k]) for k in filled)} rows)')
    print(f'ledger models not in the catalog (retired or renamed): {len(retired)}')
    print(f'pending renames (id kept for a key not yet in the ledger): {len(kept)}')
    for old_key, key, ms_id in kept:
        print(f'  {ms_id}: {old_key!r} -> {key!r}')
    print(f'ids blanked for keys not in the ledger: {len(blanked)}')
    for key, ms_id in blanked:
        why = f'id of live model {holder[ms_id]!r}' if holder.get(ms_id) in models else 'not a model id in the ledger'
        print(f'  {key!r}: {ms_id} ({why})')
    print(f'new models left blank: {len(models) - len(filled) - len(kept) - len(blanked)}')
    print(f'cells changed: {changed}')

    if output == original:
        return
    if args.check:
        print(f'{args.catalog} is not up to date with the ledger', file=sys.stderr)
        sys.exit(1)
    tmp = Path(args.catalog).with_suffix('.tmp')
    tmp.write_text(output, encoding='utf-8', newline='')
    tmp.replace(args.catalog)


if __name__ == '__main__':
    main()
