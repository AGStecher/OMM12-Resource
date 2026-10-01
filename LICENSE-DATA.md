# Data license

The code in this repository (app.R, the scripts in build_*/, content_generators/,
omm12_docker/ and the other .R/.py/.sh files) is licensed under the MIT License
(see LICENSE).

The data files in data/ and www/ that were generated for this resource
(Prokka, COGclassifier, eggNOG-mapper, DeepLocPro, CRISPRCasFinder, antiSMASH,
AMRFinderPlus and 16S results, summary tables, figures and taxonomy diagrams)
are licensed under the Creative Commons Attribution 4.0 International license
(CC BY 4.0, https://creativecommons.org/licenses/by/4.0/). Please cite the
OMM12 Resource when you reuse them.

## Third-party data (not covered by CC BY 4.0)

These files come from other sources and keep their original terms:

- Genome sequences, assembly accessions and taxonomy: NCBI GenBank/RefSeq and
  NCBI Taxonomy (https://www.ncbi.nlm.nih.gov/).
- Strain phenotype fields (oxygen requirement, cell shape, DSM numbers): BacDive
  (https://bacdive.dsmz.de/).
- Genome-scale metabolic models (data/*_model.xml and the derived *_model.rds):
  gapseq models from Zimmermann & Burrichter 2025, Zenodo
  https://doi.org/10.5281/zenodo.17358311.
- Publication metadata (data/omm12_publications*.tsv): derived from a
  Dimensions (https://www.dimensions.ai/) export, screened manually.
- KEGG pathway and module names used in data/*_kegg.tsv: KEGG
  (https://www.kegg.jp/), subject to the KEGG terms of use.
