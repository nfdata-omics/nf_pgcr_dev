# nf-core/variantprioritization: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## dev

### Added

- [#94](https://github.com/nf-core/variantprioritization/pull/94) - added optional column 'sex' to samplesheet (@HomoPolyethylen)

### Changed

- [#93](https://github.com/nf-core/variantprioritization/pull/93) - Template update for nf-core/tools v4.0.2 (@HomoPolyethylen)
- [#96](https://github.com/nf-core/variantprioritization/pull/96) - moved sample-specific PCGR parameters to samplesheet (@HomoPolyethylen)
- [#100](https://github.com/nf-core/variantprioritization/pull/100) - updated modules and subworkflows (@HomoPolyethylen)

### Fixed

### Removed

### Dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| bcftools   | 1.22        | 1.23.1      |
| htslib     | 1.21        | 1.24        |
| MultiQC    | 1.32        | 1.35        |

### Parameters

| Params         | status                         |
| -------------- | ------------------------------ |
| `tumor_site`   | removed / moved to samplesheet |
| `tumor_purity` | removed / moved to samplesheet |
| `tumor_ploidy` | removed / moved to samplesheet |

## v1.0.0 - 09.04.2026 - Jane Addams

Initial release of nf-core/variantprioritization, created with the [nf-core](https://nf-co.re/) template.

Originally written by @barrydigby, @yussab and @matbonfanti. Adapted to nf-core standards by @famosab.

Team from October 2025 nf-core Hackathon: @georgiakes, @petanska, @maxibor, @eolaniru, @colorstorm, @mkatsanto, @dhtt and @famosab.

Release names are chosen from the [list of peace activits on wikipedia](https://en.wikipedia.org/wiki/List_of_peace_activists).
