# NK cells in Parkinson's disease CSF: a single-cell look for α-synuclein-clearance candidate genes

*Short project report · Prokash · 2026*

## Background
Parkinson's disease (PD) is marked by the loss of dopaminergic neurons in the substantia nigra and the build-up of
misfolded α-synuclein. Mouse work shows that natural killer (NK) cells internalise and degrade α-synuclein, and that depleting
NK cells worsens the pathology (Earls et al., 2020). In human CSF, NK cells have been reported to expand in PD, together
with a RAC1⁺ subset (Guan et al., 2022). If NK cells help clear α-synuclein, the genes they use for this should be detectable in
the NK cells found in patient CSF. They could then be measured cheaply by qRT-PCR in local cohorts.

**Aim:** to characterise NK cells in PD CSF at single-cell resolution and short-list about five NK genes for experimental
follow-up.

## Data and methods
- **Data:** CSF 10x scRNA-seq from GEO GSE141578 (Gate et al., 2021). PD and healthy-control donors were kept, and Lewy body
  dementia donors were excluded.
- **Pre-processing (scanpy):** per-sample doublet removal with scVI-SOLO; cells with < 200 genes, ≥ 20% mitochondrial reads or
  ≥ 55% ribosomal reads were removed; genes present in < 10 cells were dropped.
- **Integration and clustering:** scVI (batch = donor) on 2,000 highly variable genes, then a UMAP and Leiden graph built on
  the scVI latent space. Clusters were annotated with canonical markers (T, NK, B, monocyte, macrophage, DC and pDC panels).
- **NK genes:** Wilcoxon test of NK cells against all other cells, then against CD8⁺ T cells. A gene was kept if adjusted
  p < 0.05, log2FC > 1 and it was expressed in > 25% of NK cells.
- **NK states:** NK cells were re-clustered in the scVI latent space and each state was compared against the other.
- **Pathways:** Enrichr GO Biological Process 2023 on the NK-enriched genes, plus a Hallmark GSEA.
- **Donor-level checks:** NK frequency and NK module scores were compared per donor (Mann–Whitney), along with an NK
  pseudobulk DESeq2 analysis (6 PD vs 8 HC donors with ≥ 5 NK cells).

## Results
1. **NK cells are a distinct CSF population.** 504 of 19,265 immune cells (2.6%) were NK cells. They express 75 NK-enriched
   genes, and 58 of these are also higher than in CD8⁺ T cells. The top genes were *CTSW*, *TYROBP*, *XCL1*, *NKG7*, *KLRD1*
   and *GNLY* (Fig. 1).
2. **Two NK states.** One is a SELL⁺IL7R⁺TCF7⁺ state (74%, CD56^bright-like). The other is a CCL5⁺GZMA/H/M⁺PRF1⁺ cytotoxic
   effector state (26%) (Fig. 2).
3. **Function.** GO enrichment is led by NK-mediated cytotoxicity, NK activation, interferon response and antiviral
   defence (Fig. 3a).
4. **Candidate genes for follow-up** (`results/candidate_genes.csv`):

| Gene | log2FC | % NK | % other | Rationale (literature) |
|---|---|---|---|---|
| *TYROBP* | 5.4 | 95 | 11 | DAP12; the TREM2–TYROBP axis drives phagocytic clearance |
| *FCER1G* | 4.5 | 85 | 10 | Fc-receptor γ chain, a possible route for uptake of antibody-bound α-synuclein |
| *CD63* | 2.4 | 67 | 28 | late endosome/lysosome, the pathway NK cells use to degrade α-synuclein |
| *CTSW* | 3.7 | 96 | 45 | NK/CTL-restricted cathepsin, the top NK gene |
| *RHOC* | 1.6 | 46 | 22 | Rho GTPase, in the same family as RAC1 (migration and engulfment) |

## Limitations
- Few donors and few NK cells (6 PD and 9 HC donors; 7–55 NK cells each). NK frequency was nominally higher in PD
  (p = 0.05), but no NK gene or pathway differed between PD and HC after correction. The candidates are therefore
  NK-specific, not PD-specific.
- The links to α-synuclein are literature-based hypotheses. Nothing here measures uptake or degradation.
- All of this is from one public, non-Bangladeshi cohort.

## Planned validation
The planned validation is qRT-PCR of the five candidates (plus NK reference genes *NKG7* and *KLRD1*) in paired blood and CSF
from PD patients and controls at NINS, Dhaka. The full design is in `../proposal/Neuroprotective_Cells_PD_Proposal.docx`.
Sample access was arranged, but funding was not secured, so this step has not been done.

## References
1. Earls RH et al. *PNAS* 2020 — NK cells clear α-synuclein and their depletion exacerbates synuclein pathology.
2. Gate D et al. *Science* 2021;374:868–874 — CD4⁺ T cells in Lewy body dementia (source of GSE141578).
3. Guan L et al. *Front Immunol* 2022;13:992505 — RAC1⁺ NK cells in PD CSF.
