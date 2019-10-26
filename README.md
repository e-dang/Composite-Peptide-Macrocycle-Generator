# UNDER DEVELOPMENT

Publication pending.

## Info
Python3.6.8 program for generating a macrocycle library using RDKit, and subsequently filtering out synthetically intractable compounds as predicted by RegioSQM (carbon EAS reactivity) and Jaguar (heteroatomic EAS reactivity). Also contains a conformational sampling tool for macrocycles that is an RDKit implementation of the algorithm used in ConfBuster<sup>[4]</sup>. A stand alone version of this conformational sampling tool can be found

The program is made to mimic the chemistry developed in our lab's prior experimental work<sup>[5][6][7][8][9][10]</sup> on synthesizing peptidomimetic macrocycles using a templated based approach, where the template in this case is the cinnamoyl cation. The macrocycle library is thus generated sequentially in the following manner: heterocycles -> sidechains -> monomers -> peptides -> template-peptide oligomers -> macrocycle. The library of heterocycles used in this program were hand selected from the VEHICLe database<sup>[11]</sup> based on percieved likelihood of undergoing EAS reactions. Also included are all natrual (D, L) amino acids, and a set of modified prolines. Side chain lengths are also made to vary from methylene, ethyl, and propyl. Alpha, beta2, and beta3 amino acid backbones are also included.

A stand-alone version of our RDKit implementation of ConfBuster's algorithm can also be found [here](https://github.com/e-dang/Macrocycle-Conformer-Generator).

## Dependencies (for full list see environment.yml):
- Python 3.6.8
- RDKit<sup>[1]</sup>
- RegioSQM<sup>[2]</sup>
- OpenBabel<sup>[3]</sup>

## Citations:
- [1] https://www.rdkit.org
- [2] Kromann, Jimmy, et al. “Fast and Accurate Prediction of the Regioselectivity of Electrophilic Aromatic Substitution Reactions.” 2017, doi:10.26434/chemrxiv.5435935.
- [3] http://openbabel.org/wiki/Main_Page
- [4] Barbeau, Xavier, et al. “ConfBuster: Open-Source Tools for Macrocycle Conformational Search and Analysis.” Journal of Open Research Software, vol. 6, 2018, doi:10.5334/jors.189.
- [5] Lawson, K. V; Rose, T. E.; Harran, P. G. Template-Constrained Macrocyclic Peptides Prepared from Native, Unprotected Precursors. Proc. Natl. Acad. Sci. 2013, 110 (40), E3753--E3760.
- [6] Lawson, K. V; Rose, T. E.; Harran, P. G. Template-Induced Macrocycle Diversity through Large Ring-Forming Alkylations of Tryptophan. Tetrahedron 2013, 69 (36), 7683–7691.
- [7] Rose, T. E.; Lawson, K. V; Harran, P. G. Large Ring-Forming Alkylations Provide Facile Access to Composite Macrocycles. Chem. Sci. 2015, 6 (4), 2219–2223.
- [8] Rose, T. E.; Curtin, B. H.; Lawson, K. V; Simon, A.; Houk, K. N.; Harran, P. G. On the Prevalence of Bridged Macrocyclic Pyrroloindolines Formed in Regiodivergent Alkylations of Tryptophan. Chem. Sci. 2016, 7 (7), 4158–4166.
- [9] Zhao, H.; Negash, L.; Wei, Q.; LaCour, T. G.; Estill, S. J.; Capota, E.; Pieper, A. A.; Harran, P. G. Acid Promoted Cinnamyl Ion Mobility within Peptide Derived Macrocycles. J. Am. Chem. Soc. 2008, 130 (42), 13864–13866.
- [10] Curtin, B. H.; Manoni, F.; Park, J.; Sisto, L. J.; Lam, Y.-H.; Gravel, M.; Roulston, A.; Harran, P. G. Assembly of Complex Macrocycles by Incrementally Amalgamating Unprotected Peptides with a Designed Four-Armed Insert. J. Org. Chem. 2018, 83 (6), 3090–3108.
- [11] Pitt, William R., et al. “Heteroaromatic Rings of the Future.” Journal of Medicinal Chemistry, vol. 52, no. 9, 2009, pp. 2952–2963., doi:10.1021/jm801513z.
