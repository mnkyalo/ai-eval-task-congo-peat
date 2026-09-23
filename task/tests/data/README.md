# Core MBK-2 -- pollen and chronology dataset

`site_metadata.txt`
: Site description, coordinates, coring date, sampling and counting protocol.

`radiocarbon_dates.csv`
: `lab_code`, `depth_cm`, `material`, `pretreatment` (laboratory preparation
  carried out on that sample), `c14_age_bp` (conventional radiocarbon
  age, yr BP), `c14_error` (1 sigma, yr), `c_to_n` (atomic C:N of the dated
  fraction), `d13c_permille`. Uncalibrated.

`pollen_counts.csv`
: `count_id`, `depth_cm`, one column per pollen taxon (raw grains counted),
  `lycopodium_count` (exotic marker grains counted alongside), and `notes`
  from the counting record.

`spike_and_volume.csv`
: `depth_cm`, `sample_volume_cm3` (volume of peat prepared),
  `lycopodium_tablets` (tablets added to that sub-sample), `tablet_batch`,
  `spores_per_tablet` and `spores_per_tablet_sd` for the batch.

`curves/`
: The laboratory's calibration curve library, in the standard `.14c` format:
  `cal BP, 14C age, sigma, Delta14C, sigma`, with comment lines prefixed `#`.
