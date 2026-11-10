# HMD data directory

Do **not** commit HMD archives to this repository.

For exact manuscript reproduction, download the HMD **Death Counts** archive for the snapshot dated **27 Aug 2026** from:

<https://www.mortality.org/Data/ZippedDataFiles>

Rename the archive to:

```text
HMD_Deaths.zip
```

and place it in this directory.

The archive must contain exactly one member whose path ends in `Deaths_1x1.txt`. The code reads that member directly and checks that its metadata reports `Last modified: 27 Aug 2026`.

No HMD exposure file is needed for the analyses in this repository.

See `../../docs/DATA.md` for licensing/citation notes.

For an exact file-content check, the manuscript `Deaths_1x1.txt` has MD5:

```text
08be9302bb6355fda6d9d89553507092
```
