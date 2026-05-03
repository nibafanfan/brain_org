# Project proposal

## Background

The Human Neural Organoid Cell Atlas (HNOCA) established a foundational
benchmark for brain organoid protocols by integrating 36 single-cell RNA
sequencing datasets across 26 differentiation methods and systematically
comparing cellular composition to primary human brain references. However,
HNOCA mostly evaluated single-lineage protocols deriving from pluripotent stem
cells, which produce primarily neuroectodermal cell types while lacking
vascular endothelium, immune cells, and multi-regional integration present in
the developing brain. A growing number of multi-lineage protocols now
incorporate vascularization, multi-region fusion, and enhanced regional
patterning, but no systematic benchmark has evaluated whether these approaches
genuinely improve correspondence with primary brain cellular diversity.

## Aim

We propose an integrated benchmark of multi-lineage human brain organoid
protocols, extending the HNOCA framework to evaluate emerging differentiation
methods published in 2024–2025. Our benchmark will integrate scRNA-seq data
from multiple multi-lineage protocol categories. Examples include:

1. Multi-region organoids combining cerebral, mid-hindbrain, and endothelial
   systems [1].
2. Vascularized midbrain organoids engineered with diffusible scaffolds [2].
3. Neuroimmune assembloids incorporating iPSC-derived microglia into cortical
   organoids [3].

We will train deep learning models (scVI, SAE) on primary human brain
reference data and quantify organoid-to-primary correspondence using cell type
mapping accuracy, embedding alignment, and cell type coverage metrics.

## Questions

1. Do multi-lineage protocols show improved cell type coverage compared to
   single-lineage HNOCA protocols?
2. Do these added cell types show transcriptomic correspondence with their
   primary brain counterparts?
3. What cell type gaps remain even in multi-lineage systems, informing future
   protocol development priorities?

## Outcome

The proposed framework will establish quantitative criteria for evaluating
multi-lineage brain organoid protocols and provide a computational pipeline
for systematic protocol comparison. By extending HNOCA to emerging
multi-lineage methods, this work will support evidence-based protocol
selection for MPS applications requiring specific cell type coverage,
particularly for studies involving neuroinflammation, blood-brain barrier
modeling, or regional brain specificity.

## References

1. Kshirsagar, A., Mnatsakanyan, H., Kulkarni, S. et al. (2025). Advanced
   Science 12, e03768. doi:10.1002/advs.202503768
2. Cai, H., Tian, C., Chen, L. et al. (2025). Cell Stem Cell 32, 824-837.e5.
   doi:10.1016/j.stem.2025.02.010
3. Kalpana, K., Rao, C., Semrau, S. et al. (2024). Methods in Molecular
   Biology 2951, 139-158. doi:10.1007/7651_2024_554
