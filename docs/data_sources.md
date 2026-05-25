# Data sources

All accessions confirmed via "Data Availability" statements unless noted.
Last updated: 2026-05-15 (Phase B schema migration recovered all 8 originally-failed deposits; corrected gse253230 cell count).

---

## Compiled atlas — 137 deposits (2026-05-15)

Final state of the GEO/SRA compilation pipeline. Each row is one deposit converted to a single AnnData h5ad in `data/processed/` (v1) and `data/processed_v2/` (HNOCA-compatible). Master manifests: `data/manifest.tsv`, `data/manifest_v2.tsv`.

**Totals:** 137 deposits | 10,905,385 cells | 7,641,866 control cells (70.1%) | 1,032 samples (722 control)

The detailed prose entries for the historical anchor papers (A1–M plus ad-hoc additions) appear below this table; they remain the source of truth for paper-to-accession pairings and any per-deposit quirks not captured by the manifest. The table here is a one-line view of the entire processed collection.

| Accession | Slug | Organoid type | Cells | Ctrl cells | Ctrl/All | Filter | Title (truncated) |
|-----------|------|--------------|------:|----------:|---------:|--------|-------------------|
| GSE303735 | `gse303735_gwi` | 3mo_brain | 54,478 | 27,758 | 4/8 | authors | Unraveling Single-Cell Regulatory Networks in iPSC-Derived Neurons fro |
| GSE197887 | `gse197887_typical_autism` | 8wk_forebrain | 77,819 | 40,796 | 5/10 | unknown | Single-cell RNA sequencing of human forebrain organoids from typical a |
| GSE241453 | `gse241453_apoe3` | apoe_organoid | 110,625 | 16,070 | 1/8 | authors | APOE3Christchurch modulates tau phosphorylation and β-catenin/Wnt sign |
| GSE266155 | `gse266155_alg13` | brain | 47,231 | 23,866 | 3/6 | authors | Deciphering the Neurological Puzzle in ALG13-CDG Through Cutting-Edge |
| GSE295097 | `gse295097_radiation` | brain | 33,562 | 5,675 | 1/3 | authors | Temporal Mapping of Radiation-Induced Changes in Human Brain Organoids |
| GSE283473 | `gse283473_ws` | brain_d56 | 96,969 | 96,969 | 6/6 | authors | Gene expression profile at single cell level of brain organoid at day5 |
| GSE301348 | `gse301348_h5n1` | brain_organoid | 70,089 | 19,962 | 2/10 | authors | Contemporary highly pathogenic avian influenza (H5N1) viruses retain n |
| GSE273907 | `gse273907_oxygen` | brain_oxygen | 278,137 | 76,857 | 2/7 | authors | Gene expression responses to oxygen challenge in a panel of human brai |
| GSE145306 | `gse145306` | brainstem | 2,345 | 2,345 | 1/1 | authors | Brainstem organoids from human pluripotent stem cells contain neural c |
| GSE150153 | `gse150153_ipsc` | cerebellar | 2,512 | 2,512 | 2/2 | authors | Human iPSC-derived organoids recapitulate development of the cerebellu |
| GSE263652 | `gse263652_cerebellar` | cerebellar | 214,075 | 214,075 | 16/16 | authors | Human cerebellar organoids to study medulloblastoma pathogenesis and t |
| GSE310490 | `gse310490` | cerebellar | 87,315 | 14,149 | 4/9 | authors | PTEN variant and genetic backgrounds combine to modify cerebellar neur |
| GSE108571 | `gse108571` | cerebral | 15,239 | 15,239 | 4/4 | authors | Single-cell RNA sequencing of vascularized brain organoids |
| GSE113089 | `gse113089` | cerebral | 3,491 | 3,491 | 2/2 | authors | Spontaneous functional network activity in organoids resembles program |
| GSE117512 | `gse117512` | cerebral | 57,707 | 57,707 | 1/1 | authors | Inhibition of BET proteins rescues neural defects in Rett syndrome [sc |
| GSE118697 | `gse118697` | cerebral | 66,234 | 66,234 | 10/10 | min500c/200g | Modeling human telencephalic development and autism-associated SHANK3 |
| GSE119861 | `gse119861` | cerebral | 20,418 | 20,418 | 14/14 | authors | Single cell gene expression analysis of stem cell-derived cerebral org |
| GSE122342 | `gse122342` | cerebral | 11,277 | 11,277 | 1/1 | authors | Fusion of Human Brain Organoids Assembles Bidirectional Connectivity B |
| GSE124031 | `gse124031` | cerebral | 1,873 | 724 | 1/2 | authors | Altered migratory trajectories in cerebral organoids derived from indi |
| GSE124174 | `gse124174` | cerebral | 13,280 | 13,280 | 1/1 | authors | Single cell RNA-sequencing of 75-day old human H1 and H9 embryonic ste |
| GSE131434 | `gse131434` | cerebral | 1,566 | 1,042 | 1/2 | authors | Virus-induced apoptosis in human cerebral organoids occurs primarily i |
| GSE132105 | `gse132105` | cerebral | 23,353 | 23,353 | 2/2 | authors | Transcriptome regulation by ALK in cerebral organoids revealed by sing |
| GSE133894 | `gse133894` | cerebral | 9,000 | 9,000 | 3/3 | authors | The Parkinson's disease associated mutation LRRK2-G2019S affects astro |
| GSE137877 | `gse137877` | cerebral | 10,985 | 10,985 | 6/6 | authors | An optimized platform to generate human cerebral organoids |
| GSE137941 | `gse137941` | cerebral | 7,287 | 7,287 | 1/1 | min500c/200g | Sliced Human Cortical Organoids for Modeling Distinct Cortical Neurona |
| GSE142143 | `gse142143` | cerebral | 33,939 | 33,939 | 1/1 | authors | NGLY1-deficient cerebral organiods display abnormal neuronal different |
| GSE146878 | `gse146878` | cerebral | 26,007 | 11,241 | 1/2 | authors | FMRP-Mediated Gene Regulation in Human Brain Development |
| GSE147047 | `gse147047` | cerebral | 24,503 | 13,192 | 1/2 | authors | Loss of function of the mitochondrial peptidase PITRM1 induces Alzheim |
| GSE157019 | `gse157019` | cerebral | 10,651 | 10,651 | 2/2 | authors | Electrophysiological Maturation of Cerebral Organoids Correlates with |
| GSE157525 | `gse157525` | cerebral | 37,665 | 17,872 | 3/6 | authors | SMARCB1 loss interacts with neuronal differentiation state to block ma |
| GSE163018 | `gse163018` | cerebral | 49,534 | 26,975 | 3/5 | authors | Single cell epigenomics reveals mechanisms of human cortical developme |
| GSE164089 | `gse164089` | cerebral | 27,082 | 16,100 | 2/4 | authors | Modeling sporadic Alzheimer’s disease in human brain organoids under s |
| GSE165577 | `gse165577` | cerebral | 58,586 | 30,691 | 3/6 | authors | Identification of neural oscillations and epileptiform changes in huma |
| GSE165975 | `gse165975` | cerebral | 11,085 | 11,085 | 1/1 | authors | In vivo development and single-cell transcriptome profiling of human b |
| GSE168323 | `gse168323` | cerebral | 123,294 | 123,294 | 27/27 | authors | Single cell transcriptomics captures features of developing and mature |
| GSE171263 | `gse171263` | cerebral | 27,910 | 20,418 | 1/2 | min500c/200g | Single cell profiling of the cerebral organoids derived from Schinzel- |
| GSE175719 | `gse175719` | cerebral | 67,203 | 44,298 | 4/6 | authors | Human cortical organoids with engineered microglia-like cells [scRNA-S |
| GSE180945 | `gse180945` | cerebral | 2,659 | 2,659 | 1/1 | authors | Neuro-immune organoid model reveals a role of microglia in cell stress |
| GSE182224 | `gse182224` | cerebral | 28,001 | 28,001 | 11/11 | authors | Structural variation at the ZNF558 locus controls a gene regulatory ne |
| GSE184409 | `gse184409` | cerebral | 12,877 | 12,877 | 6/6 | authors | Induction of inverted morphology in brain organoids by vertical-mixing |
| GSE184878 | `gse184878` | cerebral | 27,786 | 27,786 | 7/7 | authors | Cell-specific neuropathology and multiple morphogenic mechanisms in 3D |
| GSE185052 | `gse185052` | cerebral | 1,152 | 576 | 2/4 | authors | Multi-dimensional modeling disrupted synapse formation underlying psyc |
| GSE186814 | `gse186814` | cerebral | 28,116 | 19,304 | 2/4 | authors | Single cell RNA-seq post valproic acid exposure in human forebrain org |
| GSE187877 | `gse187877` | cerebral | 6,899 | 6,899 | 1/1 | authors | Androgens increase excitatory neurogenic potential in human brain orga |
| GSE189535 | `gse189535` | cerebral | 63,052 | 63,052 | 16/16 | authors | Cell-Type-Specific Impact of Glucocorticoid Receptor Activation on the |
| GSE190815 | `gse190815` | cerebral | 154,144 | 107,765 | 3/6 | authors | Maturation and circuit integration of transplanted human cortical orga |
| GSE192405 | `gse192405` | cerebral | 77,804 | 77,804 | 13/13 | authors | Single-cell transcriptional and functional analysis of human dopamine |
| GSE195510 | `gse195510` | cerebral | 140,709 | 140,709 | 1/1 | authors | FOXP1 Orchestrates Neurogenesis in Human Cortical Basal Progenitors |
| GSE195666 | `gse195666` | cerebral | 7,235 | 3,180 | 1/2 | authors | Human iPSC-derived cerebral organoids reveal progenitor pathology in E |
| GSE195692 | `gse195692` | cerebral | 72 | 48 | 12/18 | authors | Human PSCs determine the competency of cerebral organoid differentiati |
| GSE196423 | `gse196423` | cerebral | 18,129 | 18,129 | 3/3 | authors | Silk scaffolding drives self-assembly of functional and mature human b |
| GSE198927 | `gse198927` | cerebral | 14,276 | 14,276 | 2/2 | min500c/200g | Mitochondrial HSF1 triggers mitochondrial dysfunction and neurodegener |
| GSE207749 | `gse207749` | cerebral | 5,638 | 5,638 | 3/3 | min500c/200g | Functional neuronal circuitry and oscillatory dynamics in human brain |
| GSE208418 | `gse208418` | cerebral | 70,386 | 33,956 | 37/76 | min500c/200g | Human Tau Mutations in Cerebral Organoids Induce a Progressive Dyshome |
| GSE208438 | `gse208438` | cerebral | 22,204 | 22,204 | 3/3 | authors | Depressive patient-derived GABA interneurons reveal abnormal neural ac |
| GSE208710 | `gse208710` | cerebral | 43,378 | 43,378 | 3/3 | authors | Gene expression profile at single cell level of human cerebral organoi |
| GSE210720 | `gse210720` | cerebral | 42,321 | 8,519 | 1/4 | authors | Generation of ventralized thalamic organoids with inhibitory thalamic |
| GSE214422 | `gse214422` | cerebral | 97,606 | 66,372 | 8/12 | authors | Single cell gene expression profiles for isogenic PTEN panel iPSC-deri |
| GSE214538 | `gse214538` | cerebral | 18,823 | 8,547 | 1/2 | authors | FEZ1 participates in human embryonic brain development by modulating n |
| GSE219245 | `gse219245` | cerebral | 4,742 | 4,742 | 9/9 | authors | In vitro Modeling of the Human Dopaminergic System using spatially arr |
| GSE220690 | `gse220690` | cerebral | 17,989 | 8,510 | 1/2 | authors | Microglia gravitate toward amyloid plaques surrounded by externalized |
| GSE228315 | `gse228315` | cerebral | 247,009 | 70,877 | 8/26 | authors | Single cell transcriptomic profiling of human brain organoids reveals |
| GSE237855 | `gse237855` | cerebral | 15,110 | 15,110 | 11/11 | authors | Generation of advanced cerebellar organoids for neurogenesis and neuro |
| GSE241631 | `gse241631` | cerebral | 11,440 | 11,440 | 3/3 | authors | Modeling early phenotypes of Parkinson’s disease by age-induced midbra |
| GSE241743 | `gse241743` | cerebral | 30,070 | 18,342 | 2/4 | authors | The single cell transcriptomic profiling in iPSC derived Midbrain Orga |
| GSE253940 | `gse253940_ngb_timepoints` | cerebral | 32,016 | 19,055 | 2/4 | authors | Gene expression profile of control and NGB knockdown human cerebral or |
| GSE260532 | `gse260532_kcnj2` | cerebral | 7,098 | 7,098 | 1/1 | min500c/200g | KCNJ2 inhibition mitigates mechanical injury in human brain organoids |
| GSE273180 | `gse273180_foxg1` | cerebral | 41,565 | 25,929 | 2/4 | authors | scRNA-seq of human cerebral organoids derived from healthy and FOXG1 s |
| GSE286054 | `gse286054_meth` | cerebral | 33,647 | 17,802 | 8/14 | authors | Glial cell diversity and methamphetamine-induced neuroinflammation in |
| GSE297594 | `gse297594_mecp2` | cerebral | 164,235 | 164,235 | 4/4 | authors | MeCP2 regulates telencephalic development in human cerebral organoids |
| GSE306010 | `gse306010` | cerebral | 17,194 | 17,194 | 2/2 | authors | Increased Reproducibility of Brain Organoids through Controlled Fluid  |
| GSE320222 | `gse320222_cbp` | cerebral | 94,074 | 27,908 | 2/6 | min500c/200g | scRNA-seq of cerebral organoids derived from CBP mutant iPSCs |
| GSE324211 | `gse324211_xist` | cerebral | 127,086 | 72,864 | 2/4 | authors | Inducible XIST-mediated trisomy 21 correction uncovers a USP16-p16 sen |
| GSE86153 | `gse86153` | cerebral | 82,291 | 82,291 | 2/2 | authors | Identification of extensive cellular diversity and maturation of activ |
| GSE293664 | `gse293664_ngb65` | cerebral_65d | 16,675 | 13,636 | 1/2 | authors | Gene expression profile of control and NGB knockdown human cerebral or |
| GSE231546 | `gse231546_arid` | cerebral_arid | 54,480 | 54,480 | 1/1 | authors | ARID1B controls transcriptional programs of axon projection in an orga |
| GSE185472 | `gse185472_astrocytes` | cerebral_astrocyte | 216,708 | 216,708 | 8/8 | authors | Morphological diversification and functional maturation of human astro |
| GSE280812 | `gse280812_d30` | cerebral_d30 | 105,783 | 105,783 | 4/4 | authors | Gene expression profilling of day 30 cerebral organoids at the single |
| GSE299777 | `gse299777_endo` | cerebral_endo | 74,603 | 36,523 | 2/4 | min500c/200g | Endothelial cell-conditioned medium enhanced astrocyte differentiation |
| GSE163952 | `gse163952` | cerebral_hsv | 34,632 | 34,632 | 1/1 | authors | Neurodegeneration in herpes simplex virus 1 infected human brain organ |
| GSE224346 | `gse224346` | cerebral_pericyte | 177,929 | 75,874 | 7/15 | authors | Generation of human cerebral organoids with a structured outer subvent |
| GSE253230 | `gse253230_ube3a` | cerebral_ube3a | 151,948 | 151,948 | 9/9 | min500c/200g | Loss of UBE3A impacts both neuronal and non-neuronal cells in human ce |
| GSE309815 | `gse309815_measles` | cerebral_undirected | 24,189 | 6,233 | 2/7 | authors | Single-cell RNA-seq of undirected human cerebral organoids (HuCOs) inf |
| GSE216700 | `gse216700_rf_radiation` | cortical | 103,680 | 31,572 | 1/4 | authors | Radiofrequency radiation impairs the human cortical organoid developme |
| GSE239542 | `gse239542_d18_d35` | cortical | 20,614 | 20,614 | 4/4 | authors | Inferring pattern-driving intercellular flows from single-cell transcr |
| GSE277159 | `gse277159_lineage_trace` | cortical | 196,944 | 53,760 | 3/7 | min500c/200g | Neuronal lineage tracing from progenitors in human cortical organoids |
| GSE285074 | `gse285074_hypoxia` | cortical | 23,338 | 15,036 | 1/2 | authors | Single-cell transcriptomic analysis reveals distinct cellular response |
| GSE293717 | `gse293717_cep290` | cortical | 172,785 | 80,772 | 2/4 | min500c/200g | CEP290-deficiency disrupts ciliary axonemal architecture in human iPSC |
| GSE304516 | `gse304516_fabp7` | cortical | 38,119 | 18,958 | 1/2 | authors | FABP7 controls radial glial scaffold stability  during human cortical |
| GSE309759 | `gse309759_glucose` | cortical | 375,573 | 97,126 | 10/36 | authors | Glucose metabolism regulates human cortical cell diversification and m |
| GSE315443 | `gse315443_mbnl2` | cortical | 10,564 | 10,564 | 3/3 | authors | MBNL2 dysfunction in outer radial glial cells is associated with disru |
| GSE161550 | `gse161550` | cortical_assembloid | 8,669 | 8,669 | 4/4 | authors | ScRNAseq - Neurotransmitter signaling regulates distinct phases of multi |
| GSE231319 | `gse231319_polaroid` | cortical_assembloid | 71,221 | 38,400 | 1/2 | min500c/200g | Specification of a rostro-caudal axis in cortical assembloids through |
| GSE281622 | `gse281622_hippo` | cortical_hippocampal_assembloid | 40,162 | 20,263 | 3/6 | authors | Cortical versus hippocampal network dysfunction in a human brain assem |
| GSE98201 | `gse98201_mge_cortical` | cortical_mge | 59,236 | 59,236 | 8/8 | authors | Medial Ganglionic Eminence and Cortical Organoids Model Human Brain De |
| GSE183627 | `gse183627_kmt2d` | cortical_neural | 517,309 | 129,845 | 9/28 | min500c/200g | Investigating the regulatory role of KMT2D in neurodeveolopment using |
| GSE180122 | `gse180122_alsftd` | cortical_slice | 121,480 | 60,029 | 4/6 | authors | Single cell RNA-sequencing of long-term human control and C9ORF72 ALS/ |
| GSE171344 | `gse171344` | cerebral | 1,404,709 | 1,404,709 | 36/36 | min500c/200g | Glutamatergic dysfunction precedes neuron loss in cerebral organoids w |
| GSE305121 | `gse305121_somato` | cortical_somato | 30,925 | 30,925 | 5/5 | authors | Three-dimensional Co-culturing Reveals Human Stem Cell-Derived Somatos |
| GSE300486 | `gse300486_gaucher` | cortical_spinal | 253,228 | 125,199 | 8/16 | min500c/200g | Patient iPSC-derived brain organoids of neuronopathic Gaucher disease |
| GSE232448 | `gse232448_s6k1` | dorsal_forebrain | 72,536 | 36,973 | 4/8 | authors | Genomewide analysis of S6K1-depleted human dorsal forebrain organoids |
| GSE273941 | `gse273941_acm` | dorsal_forebrain_d90 | 20,000 | 20,000 | 2/2 | authors | Single cell RNA-sequencing of 90-day old astrocyte conditioned medium |
| GSE248480 | `gse248480_febo` | fetal_brain_expanding | 87,345 | 87,345 | 6/6 | min500c/200g | Human fetal brain self-organizes into long-term expanding organoids (s |
| GSE253889 | `gse253889_multiomic` | forebrain | 9,913 | 9,913 | 1/1 | authors | Multi-omic analysis of guided and unguided forebrain organoids reveal |
| GSE296775 | `gse296775_strada` | forebrain | 262,503 | 262,503 | 28/28 | authors | Delayed forebrain excitatory and inhibitory neurogenesis in STRADA-rel |
| GSE244281 | `gse244281_ganglionic` | ganglionic_eminence | 3,522 | 3,522 | 1/1 | authors | Patterning ganglionic eminences in developing human brain organoids us |
| GSE243015 | `gse243015_gbo` | hippocampal_pfc | 42,435 | 5,675 | 1/5 | authors | Single-cell Transcriptome Landscape and Cell Fate Decoding in Human Br |
| GSE215173 | `gse215173_hypothalamus_nsc` | hypothalamic | 11,813 | 11,813 | 1/1 | authors | Hypothalamus-specific NSCs derived from human brain organoids ameliora |
| GSE237274 | `gse237274_hypothalamus` | hypothalamus | 84,126 | 84,126 | 2/2 | authors | Novel human pluripotent stem cell-derived hypothalamus organoids demon |
| GSE251679 | `gse251679_medullary` | medullary_spinal_trigeminal | 20,069 | 20,069 | 2/2 | authors | Generation of human region-specific brain organoids with medullary spi |
| GSE312664 | `gse312664_microglia_assembloid` | microglia_brain_assembloid | 155,142 | 155,142 | 3/3 | min500c/200g | Human microglia in brain assembloids display region-specific diversity |
| GSE282644 | `gse282644_hiv` | microglia_cerebral | 12 | 12 | 1/1 | authors | Inflammatory responses revealed through HIV infection of microglia-con |
| GSE281452 | `gse281452_iMG` | microglia_cortical | 129,014 | 129,014 | 4/4 | authors | A microglia-containing cerebral organoid model to study early life imm |
| GSE216673 | `gse216673_trex1` | microglia_cortical_assembloid | 55,156 | 35,505 | 3/5 | authors | TREX1 is required for microglial cholesterol homeostasis and subsequen |
| GSE237133 | `gse237133_midbrain` | midbrain | 17,672 | 12,524 | 2/3 | authors | Single cell profile of human iPSC-derived midbrain organoids from heal |
| GSE271116 | `gse271116_pd_midbrain` | midbrain | 35,348 | 29,053 | 11/15 | authors | Single-cell transcriptomics revealed molecular vulnerability in a huma |
| GSE275820 | `gse275820_pd_midbrain` | midbrain | 231,501 | 225,963 | 3/4 | min500c/200g | Gene expression profile of healthy control organoids, a Parkinson dise |
| GSE281535 | `gse281535_brainstem` | midbrain | 34,702 | 34,702 | 1/1 | unknown | BrainSTEM: A multi-resolution fetal brain atlas to assess the fidelity |
| GSE277968 | `gse277968_4protocols` | multi_protocol | 91,260 | 91,260 | 69/69 | authors | Reconstitution of Human Brain Cell Diversity in Organoids via Four Pro |
| GSE279264 | `gse279264_multiregion` | multi_region | 96,938 | 96,938 | 12/12 | authors | Generation of human induced pluripotent stem cell-derived 3D human cor |
| GSE252522 | `gse252522_glucocorticoid` | neural_d70 | 124,636 | 43,346 | 8/32 | authors | Chronic exposure to glucocorticoids amplifies inhibitory neuron cell f |
| GSE247456 | `gse247456_ribosome` | neural_lineage | 106,083 | 54,928 | 5/10 | unknown | A programmed decline in ribosome levels governs human early neurodevel |
| GSE260711 | `gse260711_opioid` | neural_organoid | 25,510 | 10,383 | 2/5 | authors | Chronic opioid treatment arrests neurodevelopment and alters synaptic |
| GSE302350 | `gse302350_kmt2d` | neural_progenitor | 24,260 | 12,142 | 4/8 | authors | KMT2D-deficiency destabilizes lineage progression in immature neural p |
| GSE310824 | `gse310824_neuromod` | neuromodulatory_assembloid | 186,327 | 186,327 | 8/8 | authors | Human neuromodulatory assembloids to study serotonin signaling and dis |
| GSE241071 | `gse241071_montelukast` | neuronal_d60 | 91,485 | 32,808 | 3/9 | authors | Anti-asthma drug Montelukast induces autistic behaviors via disrupting |
| GSE286235 | `gse286235_nbm` | nucleus_basalis | 33,661 | 27,326 | 3/4 | min500c/200g | Construction of Human nucleus basalis Organoids and Cholinergic Projec |
| GSE242275 | `gse242275_oligo` | oligo_organoid | 12,727 | 12,727 | 1/1 | authors | Human oligodendrocyte progenitor cells mediate synapse elimination thr |
| GSE290048 | `gse290048_pineal` | pineal | 21,623 | 21,623 | 2/2 | authors | Generation of human pineal gland organoids with melatonin production f |
| GSE302899 | `gse302899_rf_bet` | radial_glia | 109,068 | 32,841 | 1/4 | authors | Radiofrequency regulates the BET-mediated pathways in radial glia diff |
| GSE181518 | `gse181518_sosrs` | single_rosette | 31,529 | 31,529 | 5/5 | authors | Self-organizing Single-Rosette Brain Organoids (SOSRS) from Human Plur |
| GSE290980 | `gse290980_aso_spinal` | spinal_organoid | 117,798 | 48,630 | 4/8 | authors | Targeted Antisense Oligonucleotide Treatment Rescues Developmental Alt |
| GSE183903 | `gse183903` | striatal | 29,187 | 29,187 | 2/2 | authors | Human striatal organoids derived from pluripotent stem cells recapitul |
| GSE251684 | `gse251684_striato_nigral` | striatal | 30,145 | 30,145 | 2/2 | min500c/200g | Construction of human 3D striato nigral circuitoids to recapitulate me |
| GSE220085 | `gse220085` | telencephalic | 9,458 | 9,458 | 1/1 | authors | Human telencephalic organoid development in the presence and absence o |
| GSE245719 | `gse245719_thalcortx` | thalamocortical | 142,210 | 142,210 | 31/31 | authors | Thalamocortical organoids reveal axonogenesis phenotypes of 22q11.2 mi |
| GSE325956 | `gse325956` | thalamocortical | 143,739 | 143,739 | 3/3 | authors | Thalamic NRXN1-Mediated Input to Human Cortical Progenitors Drives Exc |
| GSE287254 | `gse287254_3dprinted_vasc` | vascularized_cerebral | 388,584 | 215,382 | 4/8 | authors | Impact of 3D-printed perfusable synthetic vasculature on stem-cell bas |

---

## Quarantined / deferred (not in atlas)

Deposits we triaged but did NOT include in the final atlas. Each has a `REASON.md` in its directory.

### `data/_pending_fastq_reprocess/` — need FASTQ reprocessing for raw counts

- **gse124299**: deposit ships only normalized/processed values (no raw integer UMI counts).
- **gse131094**: SCT-scaled aggregated matrix (negative float values), no raw integer counts in deposit. No raw integer counts available in the deposit.
- **gse150903**: deposit ships only normalized/processed values (no raw integer UMI counts).
- **gse181290**: deposit ships only normalized/processed values (no raw integer UMI counts).
- **gse233295**: deposit ships only normalized/processed values (no raw integer UMI counts).
- **gse280341_msm**: normalized-only deposit; no raw integer counts
- **gse285126**: deposit ships only normalized/processed values (no raw integer UMI counts).
- **gse304918_thyroid**: deposit's `GSE304918_anndata.h5ad` ships only log-normalized counts
- **gse75140**: Smart-seq2 master data frame with normalized expression values (e.g. 6.395, 4.435 — clearly not integer UMI counts). No raw integer counts a

### `data/_pending_r_interop/` — recovered (was: need rpy2/anndata2ri to read Seurat RDS)

**Empty as of 2026-05-13.** All four RDS-only deposits (`gse277968_4protocols`, `gse219245`, `gse163952`, `gse231546_arid`) were recovered via an `Rscript → mtx → AnnData` pipeline. See `docs/atlas_compilation_notes.md` §4.6.

### `data/_blocked/` — only the 6 genuinely-blocked deposits

After the 2026-05-13 cleanup, `_blocked/` holds **only** deposits with no usable processed output. Previously-blocked deposits that were recovered keep their raw download in `data/raw/<slug>/`.

**Permanent rejects (not scRNA-seq — kept here as record-of-decision):**
- `gse97882` — ATAC-seq peaks/bed only. *Paper IS in atlas* via sibling `gse98201_mge_cortical` (59,236 cells).
- `gse135634` — bedGraph from FASTQ alignment only; no count matrix. SuperSeries `GSE106872` had only ChIP-seq + bulk + 79 SMART-seq2 NPC siblings (not organoids).
- `gse167208` — bulk RNA-seq (CLC Workbench TPM/RPKM/exon-intron metrics across 45 xlsx files). No scRNA-seq sibling.
- `gse271118_spatial` — MERFISH/Xenium spatial transcriptomics. *Paper IS in atlas* via sibling `gse271116_pd_midbrain` (35,348 cells).

**Previously deferred — now recovered (2026-05-14):**
- ~~`gse171344`~~ → **RECOVERED**: streaming-CSV per-file loader + n_counts≥500/n_genes≥200 cell-call. Final shape 1,404,709 cells × 52,584 genes (largest deposit in atlas). HTO multiplexing means all 36 samples flagged `is_control=True` at the sample level (mixed pools); cell-level demuxing not done.
- ~~`gse180122_alsftd`~~ → **RECOVERED in full**: batch 1 (105,070 cells) streamed via chunked CSV → sparse vstack; concatenated with the existing batch 2 (16,410 cells). Final shape 121,480 cells × 24,582 genes, 4/6 control samples (EpiC/EpiC2/WT/WTS1 = control; C9S1/C9S2 = C9ORF72 ALS/FTD patient).


---

## Primary brain reference — Braun et al. 2023, *Science*

Linnarsson lab, "Comprehensive cell atlas of the first-trimester developing
human brain." ~1.7M cells, 26 donors, 5–14 post-conceptional weeks, 12 major
classes / 600+ cell states. HNOCA itself uses Braun as its primary
first-trimester reference, which is the right developmental window for the
organoid protocols we're benchmarking.

- **Processed h5ad (open access — use this):** CELLxGENE collection
  https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077
- **Repo with notebooks + matrix links:** https://github.com/linnarsson-lab/developing-human-brain
- **HCA project page:** https://explore.data.humancellatlas.org/projects/cbd2911f-252b-4428-abde-69e270aefdfc
- **Raw FASTQ (controlled access — only if we need it):** EGA `EGAS00001004107`

**Status:** not yet downloaded.

---

## HNOCA atlas (He, Dony, Fleck et al. 2024, *Nature*)

Paper DOI: 10.1038/s41586-024-08172-8.

- **Integrated AnnData (recommended):** Zenodo `10.5281/zenodo.11203684`
  — `hnoca_cleanedmeta.h5ad` (18.7 GB), `hnoca_extended.h5ad` (20.2 GB),
  `disease_atlas.h5ad` (2.2 GB). ~1.7M cells, 36 datasets, 26 protocols.
- **Full intermediate metadata:** Zenodo `10.5281/zenodo.14160929` (48.9 GB).
- **Code:** https://github.com/theislab/neural_organoid_atlas — uses
  **scPoli** (not scVI). Snapseed marker hierarchy + per-figure notebooks.
- **CELLxGENE mirror:**
  https://cellxgene.cziscience.com/collections/de379e5f-52d0-498c-9801-0f850823c847

**Status:** not yet downloaded.

> **Modeling note:** HNOCA's published integration is scPoli, not scVI.
> Decision still pending: scVI on Braun (proposal default) vs scPoli for
> HNOCA-consistency.

---

## Original 3 anchor papers (from the proposal)

### A1. Kshirsagar et al. 2025, *Adv. Sci.* — multi-region organoid

- **Accession:** GEO `GSE288165` (BioProject `PRJNA1216308`)
- **Paper DOI:** 10.1002/advs.202503768
- **Experiment was snRNA-seq.** The GEO landing page describes 14 samples
  on a 10x platform, library_source=`transcriptomic single cell`. The
  underlying biology is single-cell.
- **⚠ Deposit problem.** What's actually on GEO (verified by inspecting
  every file 2026-05-03):
  - `GSE288165_deseq2_star_salmon_length_scaled_counts.tsv.gz` (1.4 MB) —
    per-sample matrix: 49,133 genes × **5 samples** (only `JHUBK2/3/4/5/7`).
    Length-scaled counts from a STAR + salmon + DESeq2 pipeline (cell-unaware).
  - `GSE288165_series_matrix.txt.gz` (3.3 KB) — metadata block only; the
    data table between `!series_matrix_table_begin / end` is **empty**
    (header line, zero data rows).
  - `GSE288165_family.soft.gz`, `GSE288165_family.xml.tgz` — series metadata.
  - Every GSM has `Sample_supplementary_file_1 = NONE`.
- **Why this matters for the benchmark:** the deposited TSV is post-aggregation
  per-sample sums (no cell barcodes, no UMIs, no per-cell vector). It is
  technically scRNA-seq-derived but cell information is gone. The paper's
  headline samples — the **9 NIBSC8 / WT4 fusion samples**, including the
  MRBO multi-region fusion organoid that's the entire point of the paper —
  are listed in the series but have **no deposited data of any kind**, in
  any file. Only 5 JHUBK samples are present in the TSV.
- **Sample list (from series_matrix.txt):**

  | Sample | TSV present? |
  |---|---|
  | NIBSC8 Cerebral | no |
  | NIBSC8 Cerebral + Mid/hindbrain | no |
  | NIBSC8 MRBO (the fusion) | no |
  | NIBSC8 Endothelial | no |
  | NIBSC8 Endothelial + Mid/hindbrain | no |
  | NIBSC8 Endothelial + Cerebral | no |
  | NIBSC8 Mid/hindbrain | no |
  | WT4 Cerebral | no |
  | WT4 Endothelial | no |
  | JHUBK2, JHUBK3, JHUBK4, JHUBK5, JHUBK7 | yes (TSV columns) |

- **Local copy:** `data/raw/gse288165_kshirsagar_pseudobulk/` (the 4 series-
  level files, ~1.4 MB total). Useful only for sample metadata / gene
  namespace; cannot drive cell-level analysis.
- **Paths to actual per-cell data:** (a) email Annie Kathuria
  (`kathu003@gmail.com`, JHU Biomedical Eng); (b) reprocess FASTQ from SRA
  `PRJNA1216308` — 14 SRR runs, ~122 GB total, multi-day cellranger;
  (c) Zenodo/figshare search for a companion AnnData (not yet done).

### A2. Cai et al. 2025, *Cell Stem Cell* — vascularized midbrain

- **Accession:** GEO `GSE276401` (also seen as Cell PII S1934-5909(25)00049-9)
- **Paper DOI:** 10.1016/j.stem.2025.02.010
- **Local:** `data/raw/cai2025/` (256 MB). 3 samples — `perfusable`,
  `ctrl`, `fentanyl`. Standard 10x trio per sample.

### A3. Kalpana et al. 2024, *MMB* — neuroimmune assembloid

- **Paper DOI:** 10.1007/7651_2024_554
- **No novel scRNA-seq deposited** (methods chapter only).
- **Action:** replace with a primary neuroimmune paper. Candidates: Buonfiglioli/De
  Witte 2024 (PMC11142229), Paşca lab assembloid papers, etc.

---

## Additional papers / datasets (added 2026-05-03)

### B. Bertacchi et al. 2024, *eLife* — FGF8 telencephalic

- **Accession:** GEO `GSE276558` (BioProject `PRJNA1157844`)
- **Primary paper:** Bertacchi M, Maharaux G, Loubat A, Jung M, Studer M.
  "FGF8-mediated gene regulation affects regional identity in human cerebral
  organoids." *eLife* 13:e98096 (2024-11-01). DOI `10.7554/eLife.98096`.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11581432/
- **Methods companion:** Bio-protocol (PMC11753195).
- **Platform:** 10x Genomics Chromium, Illumina HiSeq 4000.
- **Samples (4):** GSM8501109 WNTi (batch A), GSM8501110 WNTi (batch B),
  GSM8501111 WNTi+FGF8 (batch A), GSM8501112 WNTi+FGF8 (batch B). Treatment =
  XAV-939 WNT inhibition ± FGF8 from day 10–11, harvest day 50–60.
- **Local:** `data/raw/bertacchi2024_fgf8/` (268 MB). Lab aliases MLSR22–25
  appear to map onto the four GSMs — verify before merging.
- **Use:** single-lineage cerebral + regional patterning baseline.

### C. González-Sastre et al. 2024, *J. Tissue Eng.* — 2D-derived cortical organoids

- **Accessions:** GEO `GSE242329` (BioProject `PRJNA1013104`) and GEO
  `GSE266667` (same group, second timepoint).
- **Primary paper:** González-Sastre R, López-Alonso V, Liste I. "Efficient
  generation of human cerebral organoids directly from adherent cultures of
  pluripotent stem cells." *J. Tissue Eng.* 15 (2024). DOI
  `10.1177/20417314231226027`. PMID `38343770`.
  https://journals.sagepub.com/doi/10.1177/20417314231226027
- **Methods companion:** STAR Protocols, PMID `38933177` (PII S266616672500084X).
- **Platform:** 10x Genomics, Illumina NovaSeq 6000.
- **Samples:** GSE242329 → 1 sample (GSM7758184, "AND2 hCOs scRNAseq", 4-week);
  GSE266667 → GSM8253609 (45-day).
- **Locals:**
  - `data/raw/bertacchi2024_2d_cortical/` (140 MB) — GSM7758184, 4-week hCO.
    *Directory name is misleading — this is González-Sastre, not Bertacchi.*
  - `data/raw/bertacchi2024_45d/` (57 MB) — GSM8253609, 45-day hCO. Same naming
    issue.
- **Use:** AND2-line cerebral organoids; paper claims oligodendrocyte
  precursors, astrocytes, microglia, *and endothelial cells* appear at 4 weeks.
  Verify that claim against primary refs before treating as multi-lineage.

### D. Lampersperger et al. 2025, bioRxiv — mechanical impact

- **Accession:** ENA `PRJEB103796` (FASTQ-only)
- **Paper:** doi:10.1101/2025.01.08.631895
- **Skipped:** not multi-lineage (mechanical perturbation of cerebral organoid);
  FASTQ-only would need cellranger reprocessing.

### E. Wang et al. 2025, bioRxiv — endothelial + neural organoid

- **Paper:** doi:10.1101/2025.05.20.653559
- **Highly relevant** (multi-lineage by design — embryoid-body-derived
  endothelial + neural).
- **Status:** no scRNA-seq accession found yet. Page is Cloudflare-walled to
  curl. Manual investigation needed (or contact authors).

### F. Ullah/Shcheglovitov 2026, *Nature Protocols* — single neural rosette

- **Accession:** GEO `GSE210960`
- **Paper:** doi:10.1038/s41596-025-01197-x
- **Skipped:** single-lineage telencephalic SHANK3 disease modeling — out of
  scope for the multi-lineage benchmark.

### G. PMC11753195-paired — Pham et al. 2019, vascularized brain organoid

- **Accession:** GEO `GSE134049`
- **Paper title:** "Development of human brain organoids with functional
  vascular-like system" (PubMed 31591580)
- **Local:** `data/raw/gse134049_vascular_2019/` (177 MB). Series-level
  combined matrix (`GSE134049_barcodes/genes/matrix.gz`); 2 samples merged.
- **Use:** older but directly multi-lineage (vasc + neural) — useful as a
  pre-2024 vascularization reference. *Verify the paper-to-GEO pairing*: the
  PMC11753195 ID in the user's list is later than this 2019 dataset.

### H. Sun et al. 2022, *eLife* — vascularized brain organoid + microglia

- **Title:** "Generation of vascularized brain organoids to study
  neurovascular interactions"
- **Authors:** Sun XY, Ju XC, Li Y, Zeng PM, Wu J, Zhou YY, Shen LB, Dong J,
  Chen YJ, Luo ZG.
- **Paper:** https://elifesciences.org/articles/76707 (PMC9246368)
- **Accession:** SRA `SRP338043` — runs `SRR15992285`, `SRR15992286`
  (BioProject `PRJNA764860`, experiments `SRX12280793` + `SRX12280794`).
- **Format:** raw 10x FASTQ deposited in SRA. **Genuine scRNA-seq** —
  recoverable via 2× cellranger runs. *Not* pseudobulk; this is
  qualitatively different from the Kshirsagar problem (Kshirsagar's GEO
  deposit collapsed the data to bulk before upload; here the per-cell data
  is intact in SRA, just needs alignment).
- **Highly relevant:** vascularization + microglia + brain organoid — multi-
  lineage by construction. Strong candidate for the benchmark.
- **Status:** not yet downloaded. Needs FASTQ pull (probe ENA size first)
  → cellranger run against GRCh38.

### I. Wang et al. 2021, *Nat. Med.* — neural-perivascular assembloid + SARS-CoV-2

- **Accession:** SRA `PRJNA668200` only — **no GEO**, no Zenodo, no figshare.
- **Paper:** "A Human 3D neural-perivascular assembloid promotes astrocytic
  development and enables modeling of SARS-CoV-2 neuropathology."
  https://www.nature.com/articles/s41591-021-01443-1
- **Platform:** 10x Genomics + CellRanger v4.0 (per paper Methods).
- **Scope:** 2 biosamples — SAMN16401133 (CO, cortical organoids) and
  SAMN16401134 (PCCO, pericyte-containing cortical organoids), each with 4
  runs (ctrl + experiment, 2 batches). Total 8 SRR runs, ~39 GB FASTQ.
- **⚠ FASTQ structure anomaly:** spot check showed symmetric 101+101 bp reads
  in `_1`/`_2`, not standard 10x 28+90 bp. Some runs have a third unsuffixed
  `.fastq.gz` — likely the real cell-barcode read. Verify before reprocessing.
- **Skipped:** raw-only, 2 samples, 2021 SARS-CoV-2 focus — outside the
  multi-lineage 2024–2025 benchmark scope.

### J. Glia 2022 (PMC9314680) — LCSB direct download

- **Data:** https://webdav.lcsb.uni.lu/public/data/cx25-ht49/ (custom layout —
  Figure 1–5 + Supplementary subdirs)
- **Paper:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9314680/
- **Status:** not yet downloaded (subdir structure needs traversal). Likely
  midbrain organoid data from the Schwamborn/Schäfer lab (LCSB Luxembourg).

### K. Klein et al. 2018 (?), *Nature Methods* — drop-seq cerebral organoid

- **Accession:** GEO `GSE110006`
- **Paper:** https://www.nature.com/articles/s41592-018-0081-4
- **Local:** `data/raw/gse110006_2018/` (25 MB). Mostly bulk FPKM tracking
  files from 2018; 1 single-cell sample (`GSM3243667_H7noid_filtered_*`,
  10x format) is the only useful piece for scRNA-seq.

### L. Pollen et al. 2019 (?), *Nature Neuroscience* — Smart-seq cerebral

- **Accession:** GEO `GSE115011`
- **Paper:** https://www.nature.com/articles/s41593-018-0316-9
- **⚠ Format issue:** 296 GSMs of `*_sorted.bam.txt.gz` files (~4.5 MB each,
  1.34 GB total). Not standard 10x — likely Smart-seq2 plate-based with BAM
  exports. Custom loader needed.
- **Status:** not yet downloaded — confirm relevance before pulling 1.3 GB.

### M. Nature 2025 (s41586-025-08808-3) — sensory/spinal organoid

- **Accession:** GEO `GSE251892`
- **Paper:** https://www.nature.com/articles/s41586-025-08808-3
- **Local:** `data/raw/gse251892_2025_nature/` (2.5 GB). 24 samples:
  - 4× `hSeO_*` (sensory organoid?)
  - 5× `hdSpO_*` (dorsal spinal organoid?)
  - 4× `hASA_*` (assembly/assembloid?)
  - 4× `IndividualOrganoidPooled_*`
  - 7× `hDRG_*` (dorsal root ganglion)
- **⚠ Scope question:** this is **spinal cord + DRG**, not brain. Not in the
  brain-organoid benchmark scope as currently framed. Worth keeping if the
  benchmark widens to "neural organoid" generally; otherwise drop.
