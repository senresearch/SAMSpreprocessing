# SAMSstudypreprocessing
___

# Effect of Statin Treatment on Metabolites, Lipids and Prostanoids in Patients with Statin Associated Muscle Symptoms (SAMS)
---

## Abstract


*Background*: Between 5-10% of patients discontinue statin therapy due to statin-associated adverse reactions, primarily statin associated muscle symptoms (SAMS). The absence of a clear clinical phenotype or of biomarkers poses a challenge for diagnosis and management of SAMS. Similarly, our incomplete understanding of the pathogenesis of SAMS hinders the identification of treatments for SAMS. Metabolomics, the profiling of metabolites in biofluids, cells and tissues is an important tool for biomarker discovery and provides important insight into the origins of symptomatology. In order to better understand the pathophysiology of this common disorder and to identify biomarkers, we undertook comprehensive metabolomic and lipidomic profiling of plasma samples from patients with SAMS who were undergoing statin rechallenge as part of their clinical care.

*Methods and Findings*: We report our findings in 66 patients, 27 with SAMS (cases) and 39 statin-tolerant controls. SAMS patients were studied during statin rechallenge and statin tolerant controls were studied while on statin. Plasma samples were analyzed using untargeted LC-MS metabolomics and lipidomics to detect differences between cases and controls. Differences in lipid species in plasma were observed between cases and controls. These included higher levels of linoleic acid containing phospholipids and lower ether lipids and sphingolipids. Reduced levels of acylcarnitines and altered amino acid profile (tryptophan, tyrosine, proline, arginine, and taurine) were observed in cases relative to controls. Pathway analysis identified significant increase of urea cycle metabolites and arginine and proline metabolites among cases along with downregulation of pathways mediating oxidation of branched chain fatty acids, carnitine synthesis, and transfer of acetyl groups into mitochondria.

*Conclusions*: The plasma metabolome of patients with SAMS exhibited reduced content of long chain fatty acids and increased levels of linoleic acid (18:2) in phospholipids, altered energy production pathways (b-oxidation, citric acid cycle and urea cycles) as well as reduced levels of carnitine, an essential mediator of mitochondrial energy production. Our findings support the hypothesis that alterations in pro-inflammatory lipids (arachidonic acid pathway) and impaired mitochondrial energy metabolism underlie the muscle symptoms of patients with statin associated muscle symptoms (SAMS).


## Supplemental material

The directory contains the supporting code necessary to run the analysis contained in the publication:

Timothy J. Garrett, Michelle A. Puchowicz, Qingming Dong, Gregory Farage, Richard Childress, Edwards A. Park, Joy Guingab, Claire L. Simpson , Saunak Sen, Elizabeth C. Brogdon, Logan M. Buchanan, Rajendra Raghow, Marshall B. Elam (2023) (2021) *Effect of Statin Treatment on Metabolites, Lipids and Prostanoids in Patients with Statin Associated Muscle Symptoms (SAMS) (submitted)*


The `src` directory contains all the necessary code files, including functions required for processing. The `notebooks/preprocessing` directory includes two Jupyter notebooks, [`PreprocessingMeta.ipynb`](https://github.com/senresearch/SAMSpreprocessing/blob/main/notebooks/preprocessing/PreprocessingMeta.ipynb) and [`PreprocessingLipo.ipynb`](https://github.com/senresearch/SAMSpreprocessing/blob/main/notebooks/preprocessing/PreprocessingLipo.ipynb). These notebooks were used to preprocess the metabolomic and lipidomic datasets used in our statistical analysis, respectively. These notebooks include detailed explanations of the preprocessing steps taken and serve as a reference for reproducing our results.