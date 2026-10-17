# Data provenance and reuse

This repository contains derived, harmonised SNP-level subsets and aggregate outputs used to reproduce the reported analyses. It contains no individual-level data and no authentication credential.

## Dataset identifiers

The complete source crosswalk, including verified public access pages, is in `DATASET_SOURCE_MANIFEST.csv`. The main MR analyses use the following OpenGWAS identifiers:

- `ebi-a-GCST90013923`: mood swings, UK Biobank field 1920, SPA analysis;
- `ebi-a-GCST90013959`: help-seeking, UK Biobank field 2090, SPA analysis;
- `ieu-b-102`: major-depression meta-analysis, UK Biobank plus PGC, excluding 23andMe in the public release;
- `ebi-a-GCST90013938`: UK Biobank asthma, SPA analysis;
- `finn-b-J10_ASTHMA`: FinnGen asthma; and
- `ukb-b-6991`: the help-seeking release used in the superseded original analysis (mixed scale); retained only to reproduce the scale diagnostic and the arithmetic-rescaling illustration. The repository does not state or imply that the release has been withdrawn from OpenGWAS.

The calibration analyses additionally use `ukb-a-246` and `ukb-b-14180` as linear-probability releases of the same UK Biobank fields represented by `ebi-a-GCST90013959` and `ebi-a-GCST90013923`.

The phenotype-architecture analysis uses seven public source files recorded in `scheme_c/code/ldsc_input_manifest.tsv`: the three depression-related exposures, ADHD, autism spectrum disorder, neuroticism and UK Biobank asthma. Raw genome-wide files are not redistributed. The exact files used in the completed run are identified by SHA-256 hashes in `scheme_c/provenance/SOFTWARE_AND_INPUT_PROVENANCE.md`.

## Attribution

OpenGWAS service use should cite both the OpenGWAS infrastructure report and the MR-Base platform paper, as required by the OpenGWAS end-user licence. Each underlying GWAS must also be cited according to its original provider's requirements. The access pages listed in `DATASET_SOURCE_MANIFEST.csv` and the seven direct Scheme C source files in `scheme_c/code/ldsc_input_manifest.tsv` were checked on 12 September 2026. These links establish availability and provenance; they do not override the original providers' access or redistribution conditions.

## Reuse boundary

OpenGWAS does not grant a single blanket licence over all hosted third-party datasets. Reuse and redistribution remain subject to the terms of the original dataset owners. The small derived subsets in `data/calibration/`, `data/mr/`, `case_study/direct_diagnostic/` and `scheme_c/source_version_check/` are supplied for transparent academic verification of this analysis; their presence does not relicense the underlying GWAS data.

Licences by file category are listed in `DATA_LICENSES.md`. Users who redistribute these derived subsets further are responsible for complying with the terms of the original UK Biobank, PGC, FinnGen, GWAS Catalog and OpenGWAS providers. The dataset identifiers and acquisition routes in `DATASET_SOURCE_MANIFEST.csv` allow the subsets to be re-extracted from the original sources.

The public-source ADHD file used for genome-wide LDSC matched SNP and allele identities but was not numerically identical to the older OpenGWAS `ieu-a-1183` instrument snapshot retained from the original analyses: 1 of 39 overlapping effects agreed within 1e-4. The retained instrument-level effects used for this comparison (SNP, alleles, beta and P value for six source releases, including `ieu-a-1183` and `ebi-a-GCST90029028`) are in `scheme_c/source_version_check/retained_instrument_effects.csv`; the comparison output is `scheme_c/source_version_check/Scheme_C_source_version_reverse_check.csv`. The ADHD results describe the recorded public file, not an exact re-download of the historical snapshot.

OpenGWAS policies: <https://opengwas.io/policies/eula>
