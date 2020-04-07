## Ligand- and Receptor-Based Docking ##

In [LiBELa](libela.md), we developed an very simple, but powerful algorithm for the docking of a ligand in a protein receptor binding site. This approach is based in the [MolShaCS](https://linkinghub.elsevier.com/retrieve/pii/S0223-5234(12)00682-4) algorithm, previously developed in our group, and re-implemented in LiBELa.

Briefly, the first page of the docking procedure, the ligand initial placement, is carried out by maximize the overlap of shape and charge between a *search* ligand and a *reference* ligand. This procedure is very fast and, in the end, a Hodgkin's similarity index is computed. This index is a number between -1.0 and 1.0, where -1.0 indicates a perfect anti-similarity and 1.0 a perfect similarity between those ligands.

After the *search* ligand is placed in the binding site by comparing its 3D properties with a *reference* ligand (already placed), the search ligand pose (i.e., orientation and conformation) is optimized by minimizing the interaction energy with the receptor. Here, we use the [AMBER Force Field](http://ambermd.org/AmberModels.php) together with [GAFF2](http://ambermd.org/antechamber/gaff.html) atomic parameters and ligand-provided atomic charges to compute the binding energy.

Some additional references can be found below:
1.  [Ligand- and receptor-based docking with LiBELa. Journal of Computer-Aided Molecular Design, p. 713-723, 2015](http://dx.doi.org/10.1007/s10822-015-9856-1).
2. [Towards a critical evaluation of an empirical and volume-based solvation function for ligand docking. Plos One, v. 12, p. e0174336, 2017](http://dx.doi.org/10.1371/journal.pone.0174336).
3. [Comparative Analysis of Electrostatic Models for Ligand Docking. FRONTIERS IN MOLECULAR BIOSCIENCES, v. 6, p. 52, 2019.](http://dx.doi.org/10.3389/fmolb.2019.00052).
