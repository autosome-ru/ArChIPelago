"""Cell-type overlap of the human and mouse ChIP-Seq experiments of a TF, derived by a written rule.

The first submission (results/round1_package/Manuscript_analysis/analysis/cross_species/cross_species_table.py, dict OVERLAP) labelled each TF by hand: "manual judgement of shared /
equivalent tissue between the human and mouse experiments", values yes / partial / no (a trailing '?' marked doubt; the transfer-driver analysis strips it).
That definition is applied here literally and reproducibly:
  1. every abbreviated cell label (the ABBR output of cross_species_table.py) is assigned one tissue / lineage (LIN below;
     written to 3_analysis/integration/celltype_lineage_map.csv); cell lines take their tissue of origin; labels without a known tissue
     are 'unknown' and never count as shared;
  2. shared lineages S = lineages present among the human and among the mouse experiments of the TF (sets actually used);
  3. label = 'yes' if S is not empty, 'no' if S is empty. The round-1 definition does not define 'partial', so it is not used.
Check on the first-submission cell lists (cross_species_table_round1.csv): the rule reproduces all 32 round-1 'yes' and 'no' labels;
the 4 differences are the round-1 'partial' labels (CEBPB, E2F4 '?', IRF1, JUND '?'), which become 'yes'. A majority rule
(shared lineage >= 50 % of the experiments of both species) reproduced only 25 of 36 and was rejected.
The note gives S and the share of experiments in S per species, so each label can be checked against Sup. Table 5's cell columns.
"""
import re

LIN = {
 "prostate": "LNCaP LNCaP-1F5 LNCaP-AR C4-2B DU145 VCaP RWPE-1 |prostate ca.| prostate |prostate(ERG-tg)| |prostate(Pten-/-)|",
 "epididymis": "epididymis",
 "B lymphoid": "GM12878 LCL BJAB BL41 CA46 Ramos OCI-LY10 OCI-LY3 P493-6 U-2932 |CD20+ B| Blue1 MM1.S H929 |plasma cell| 38B9(pre-B) |B cells| CH12 "
               "pre-B pre-pro-B pro-B |splenic B| plasmablasts",
 "T lymphoid": "Jurkat CCRF-CEM RPMI-8402(T-ALL) T-ALL Tconv Th1 Th2 Treg |CD4 T| |DN thymocytes| |DN1 thymocytes| |DN2b thymocytes| |DP thymocytes| |T cells| Th17 thymocytes",
 "erythroid / megakaryocytic": "K562 |CD36+ erythroid| erythroblasts proerythroblasts |CB megakaryocytes| CMK MEL G1E G1E-ER4 G1ME |GATA1s megakaryoblasts| "
               "L8057(megakaryoblast) |FL erythroblasts| |Ter119+ erythroid| megakaryocytes ES-EP",
 "myeloid": "|CD14+ monocytes| monocytes macrophages THP-1 HL-60 NB4 |APL blasts| Kasumi-1 BMDC BMDM RAW264.7",
 "haematopoietic progenitor": "CD34+ |CB CD133+| FDCP-mix HPC-7 |LSK HPC| |Lin- HPC| |BM iFLT3| |bone marrow| |diff. hematopoietic|",
 "liver": "HepG2 liver",
 "intestine / colon": "Caco-2 LoVo GP5d HCT-116 LS180 |intestinal villus| |intestinal villus(Cdx2-/-)| intestine |sm. intestine|",
 "mammary": "MCF-7 T-47D ZR-75-1 MDA-MB-134 |breast tumor| HMEC MCF10A-Src HMF AB32 H3396 mammary |mammary 3134| |mammary 7438|",
 "endometrium / uterus": "ECC-1 Ishikawa uterus",
 "adipocyte": "ASC-adipocyte MSC-adipocyte |SGBS adipocytes| adipocyte hASC 3T3-L1 |3T3-L1 adipocytes| |BAT adipocytes| |WAT adipocytes|",
 "fibroblast": "BJ-hTERT IMR90 AG09309 HFF HVMF NHLF WI-38 |fetal lung fibr.| fibroblast HCFaa MEF |MEF(OSK reprogr.)| NIH-3T3",
 "pluripotent": "H1-hESC H9-hESC ESC |diff. ESC| |F9 EC|",
 "neural": "NH-A SK-N-SH PFSK-1 ReN-VM(NPC) U87 brain cerebellum cortex |E14.5 brain| ES-NPC |neural EB| |STHdh striatal|",
 "muscle / heart": "HSMMtube C2C12 HL-1 heart |E14.5 heart|",
 "lung": "A549 lung HCC95 H2171 H128",
 "endothelium": "HUVEC", "kidney": "HEK293 kidney", "bone": "osteoblast U2OS MC3T3-E1 IDG-SW3(osteocyte)", "skin": "NHEK keratinocytes",
 "pancreas": "|panc. islets| pancreas PANC-1", "cervix": "HeLa HeLa-S3", "oesophagus": "KYSE70", "thyroid": "TT", "adrenal": "H295R",
 "pituitary": "AtT-20", "embryonic limb": "|E14.5 limb| limb", "trophoblast": "|TS cells| placenta", "testis": "testis", "spleen": "spleen",
 "unknown": "FB8470 RMWT_227 n/a n/a(FH) n/a(MH) tissue(C57BL6) |tissue(p53 mut.)|",
}
LINEAGE = {}
for lin, s in LIN.items():
    for lab in [x for x in re.findall(r"\|([^|]+)\||(\S+)", s)]:
        lab = lab[0] or lab[1]
        assert lab not in LINEAGE, lab
        LINEAGE[lab] = lin


def parse(cells):
    """'LNCaP(x8); VCaP(x8); C4-2B' -> {'LNCaP': 8, 'VCaP': 8, 'C4-2B': 1}"""
    out = {}
    for t in cells.split("; "):
        m = re.fullmatch(r"(.*)\(x(\d+)\)", t)
        lab, n = (m.group(1), int(m.group(2))) if m else (t, 1)
        if lab not in LINEAGE:
            raise SystemExit("CHECK FAILED cell label without lineage: %r" % lab)
        out[lab] = out.get(lab, 0) + n
    return out


def overlap(cells_h, cells_m):
    h, m = parse(cells_h), parse(cells_m)
    lh = {LINEAGE[k] for k in h} - {"unknown"}
    lm = {LINEAGE[k] for k in m} - {"unknown"}
    S = lh & lm
    fh = sum(n for k, n in h.items() if LINEAGE[k] in S) / sum(h.values())
    fm = sum(n for k, n in m.items() if LINEAGE[k] in S) / sum(m.values())
    lab = "yes" if S else "no"
    note = ("shared: %s; human %d/%d, mouse %d/%d experiments in shared lineages" % (", ".join(sorted(S)), round(fh * sum(h.values())), sum(h.values()),
            round(fm * sum(m.values())), sum(m.values()))) if S else "no shared lineage (human: %s; mouse: %s)" % (", ".join(sorted(lh)), ", ".join(sorted(lm)))
    return lab, note


def write_map(path):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["cell label", "tissue / lineage"])
        for k in sorted(LINEAGE, key=lambda x: (LINEAGE[x], x)):
            w.writerow([k, LINEAGE[k]])
