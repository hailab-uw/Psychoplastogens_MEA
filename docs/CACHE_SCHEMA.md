# Cache schema (v2.0)

This document defines the on-disk layout of the per-dataset cache files produced by `src/pipeline/preprocess_and_save.m` and consumed by everything downstream (figures, analysis modules, future notebooks). **Never touch raw TDT data outside the preprocessing pipeline.** Every figure and analysis script must load its inputs from these caches via `src/utils/load_cache.m`.

## Location

```
<repo>/cache/
├── cache_<safe_datasetName>.mat         # one file per dataset (this schema)
└── connectivity/                        # optional; written by src/analysis/run_connectivity.m
    └── connectivity_<study>_<safe_datasetName>.mat
```

`safe_datasetName` is the folder name with `-` and `#` replaced by `_` (see `src/utils/cache_filename.m`).

## Why spike times are cached

v1 cached only summary metrics (`rates`, `counts`, `burstRates`). Connectivity, rasters, ISI distributions, and burst-onset alignment all need the raw spike times, so v2 caches them directly. This means any analysis can be re-run in seconds without ever re-reading the raw TDT blocks.

The trade-off is disk: at ~10 min × 60 channels × few spikes/s/channel, a full dataset cache is on the order of 2–5 MB. With 18 datasets total that is well under 100 MB, which is negligible next to the raw data.

## Top-level fields

A v2.0 cache file is saved with `save(path, '-struct', cache, '-v7.3')`, so all fields below appear as **top-level variables** inside the `.mat` file.

### Metadata

| Field | Type | Description |
|---|---|---|
| `version` | `char` | Always `'2.0'`. Loaders validate this. |
| `datasetName` | `char` | Original dataset folder name (with `-` and `#` intact). |
| `createdAt` | `char` | ISO-8601 timestamp of when the cache was built. |
| `fs` | `double` | Sampling rate in Hz, taken from the TDT block. |
| `durationSec` | `double` | Recording duration in seconds (`numel(samples) / fs`). |
| `channelsUsed` | `1 x nCh` int | Channel indices actually processed (typically `1:60`; see `cfg.channels.default`). |
| `cfgUsed` | `struct` | Full snapshot of `project_config()` at build time. Lets any downstream script recover exact filter / detection parameters. |

### Spike-level (new in v2.0)

| Field | Type | Description |
|---|---|---|
| `spikeTimes` | `1 x nCh` cell | Each cell is a column vector of spike timestamps in seconds, sorted ascending. Empty `zeros(0,1)` for channels with no detected spikes or that failed to load. |

### Spike summary

| Field | Type | Description |
|---|---|---|
| `spikeCounts` | `nCh x 1` double | Total number of spikes per channel. |
| `spikeRates` | `nCh x 1` double | Firing rate in **spikes/min** per channel. The field name remains `spikeRates` for schema compatibility. `NaN` if `durationSec <= 0`. |

### Burst-level (new in v2.0)

Bursts are detected with the ISI-threshold method (see `src/utils/detect_bursts.m`). A burst is a run of spikes with within-burst ISI <= `cfg.burst.isi_max_ms` (100 ms), terminated when ISI > `cfg.burst.isi_end_ms` (200 ms), requiring at least `cfg.burst.min_spikes` (5) spikes.

| Field | Type | Description |
|---|---|---|
| `burstStartIdx` | `1 x nCh` cell | Each cell: row vector of start indices into `spikeTimes{c}` (1 entry per burst). Empty when no bursts. |
| `burstEndIdx` | `1 x nCh` cell | Matching row vector of end indices (inclusive) into `spikeTimes{c}`. |
| `burstStartTimes` | `1 x nCh` cell | Burst start timestamps in seconds, derived as `spikeTimes{c}(burstStartIdx{c})`. Stored to avoid re-derivation. |
| `burstEndTimes` | `1 x nCh` cell | Burst end timestamps in seconds. |

### Burst summary

| Field | Type | Description |
|---|---|---|
| `burstCounts` | `nCh x 1` double | Number of bursts per channel. |
| `burstRates` | `nCh x 1` double | Burst rate in **bursts/min** per channel. `NaN` if `durationSec <= 0`. |

## Deriving other quantities

Everything else commonly needed in analysis is a one-liner off this schema:

| Quantity | Formula |
|---|---|
| ISI distribution for channel `c` | `diff(cache.spikeTimes{c})` |
| Spikes per burst on channel `c` | `cache.burstEndIdx{c} - cache.burstStartIdx{c} + 1` |
| Burst durations (s) on channel `c` | `cache.burstEndTimes{c} - cache.burstStartTimes{c}` |
| In-burst spike fraction | `sum(spikesPerBurst) / cache.spikeCounts(c)` |
| Binned spike matrix (nCh × nBins) | `histcounts(cache.spikeTimes{c}, edges)` over channels |
| Raster row for channel `c` | `cache.spikeTimes{c}` — plot as ticks |

## Loading caches

**Always** go through `load_cache.m` — do not call `load(path)` directly. It:

1. Validates that the file exists (with an actionable error).
2. Checks the schema `version` string and rejects anything unknown.
3. Verifies that every required field is present.

```matlab
cfg   = project_config();
cache = load_cache('IdoDOI-230914-142502_#1', cfg);

spikes_ch12 = cache.spikeTimes{12};                 % col vector, seconds
rate_ch12   = cache.spikeRates(12);                 % spikes/min
burst_ends  = cache.burstEndTimes{12};              % row vector, seconds
```

For paired baseline/treatment access use `load_pair_cache` (two structs) or the convenience pool loaders:

```matlab
[pairs, ~] = get_pairs_and_labels(cfg, 'doi');

% Both cache structs for pair 3:
[bCache, tCache] = load_pair_cache(pairs(3), cfg);

% Pool a summary metric across all pairs, channel-aligned:
[bRates, tRates] = load_pair_metric(pairs, 1:60, 'spikeRates', cfg);

% Spike-time cells for pair 3 (e.g. for connectivity/raster work):
[bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pairs(3), 1:60, cfg);
```

## Allowed metric names (for `load_pair_metric`)

- `'spikeRates'` — spikes/min per channel
- `'spikeCounts'` — total spikes per channel
- `'burstRates'` — bursts/min per channel
- `'burstCounts'` — total bursts per channel

Anything beyond these scalar-per-channel summaries should go through `load_cache` / `load_pair_cache` and pick the field directly.

## Rebuilding after a schema change

Because caches store a `version` string, any change to the schema must bump the version and will cause `load_cache` to raise `load_cache:UnknownSchema` for legacy files. To rebuild:

```matlab
preprocess_and_save('overwrite', true);           % all datasets
preprocess_and_save('overwrite', true, 'study', 'doi');
```

## File format notes

- Saved with `-v7.3` (HDF5) so cell arrays of variable-length spike-time vectors round-trip cleanly and files can be opened from Python (via `h5py`) if needed.
- Field ordering inside the file is not significant — loaders reference fields by name.
- `cfgUsed` is a deep snapshot. If you ever need to know exactly which filter parameters produced a cached spike train, inspect `cache.cfgUsed.filter`, `cache.cfgUsed.spike`, and `cache.cfgUsed.burst`.
