# Human Mortality Database data

The real-data code uses the HMD period **Death Counts, 1x1** file. The response is the female share of all deaths in each year-age cell:

`Female / Total`, where `Total = Female + Male`.

The denominator-aware sensitivity uses the same death-count panel; no exposure file is required by this repository.

## Exact snapshot used for the manuscript

The manuscript analyses used the HMD Death Counts file with metadata:

```text
Deaths (period 1x1), Last modified: 27 Aug 2026; Methods Protocol: v6 (2017)
```

For exact reproduction, use the 27 Aug 2026 snapshot available from the HMD zipped-data page:

<https://www.mortality.org/Data/ZippedDataFiles>

Rename the Death Counts archive to `HMD_Deaths.zip` and place it at `real_data/data/HMD_Deaths.zip`, or set the `BGHMM_DATA_FILE` environment variable to its location.

The import code checks the metadata date and stops on a different snapshot by default, because HMD may revise historical values.

## Why the data are not committed

HMD asks users to obtain data directly from the HMD website so they receive the current licensing, source, processing, and citation information. HMD-constructed estimates are made available under CC BY 4.0, while original input data can remain subject to the original providers' distribution terms. Therefore this public software repository does not redistribute the HMD archive.

Current terms and citation guidance:

- <https://www.mortality.org/Data/UserAgreement>
- <https://www.mortality.org/Research/CitationGuidelines>

## Suggested acknowledgement

Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). <https://www.mortality.org/>. Data snapshot: 27 Aug 2026.

For an exact file-content check, the manuscript `Deaths_1x1.txt` has MD5:

```text
08be9302bb6355fda6d9d89553507092
```
