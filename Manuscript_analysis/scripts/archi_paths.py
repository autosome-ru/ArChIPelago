"""Locations shared by the table and figure generators of Manuscript_analysis/.

RES     Manuscript_analysis/            results table (results_table.csv), Sup_Tables/, Figures/
BASIS   Manuscript_analysis/analysis/   outputs of the per-TF analyses (model_evaluation/, mouse_transfer/, saturation/,
                                        operational/, slim_dichipmunk/, cross_species/, datasets/, motif_subtypes/)
INPUTS  Manuscript_analysis/analysis/inputs/   tables written by the pipeline notebooks and the GTRD metadata
ZENODO  folder with the Zenodo archive unpacked (PWMs_mono_HUMAN/, PWMs_di_HUMAN/, ...); environment variable
        ARCHI_ZENODO_DIR, default <repository>/14927304 (only needed by the scripts that read PWM files)
"""
import os

RES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASIS = os.path.join(RES, "analysis")
INPUTS = os.path.join(BASIS, "inputs")
ZENODO = os.environ.get("ARCHI_ZENODO_DIR", os.path.join(os.path.dirname(RES), "14927304"))
