from __future__ import annotations

import collections

import yaml
from pathlib import Path
from typing import Optional, Any, Dict, Literal
from pydantic import BaseModel, ValidationError, AnyUrl, computed_field


class Globals(BaseModel):
    threads: int = 4
    sort_memory_limit: str = "2G"


class Paths(BaseModel):
    """Holds the common high-level directory structure for the experiment"""
    root: Optional[Path] = None
    data_dir: Optional[Path] = None
    log_dir: Optional[Path] = None
    results_dir: Optional[Path] = None
    reference_genome_dir: Optional[Path] = None
    externals_dir: Optional[Path] = None


class RunSteps(BaseModel):
    basecalling: bool = False
    align: bool = False
    methylation_summary: bool = False
    deconvolution: bool = False


class PipelineControl(BaseModel):
    run_steps: RunSteps


class SetupDownloads(BaseModel):
    fast5_download_url: AnyUrl
    reference_genome_url: AnyUrl
    uxm_atlas_url: AnyUrl
    manifest_url: AnyUrl
    full_atlas_url: AnyUrl


class SetupPaths(BaseModel):
    reference_genome_name: str = "genome.fa"
    fast5_input_dir_name: str = "fast5_input"

    fast5_input_dir: Optional[Path] = None
    reference_genome_dir: Optional[Path] = None
    full_reference_genome_path: Optional[Path] = None
    def build_paths(self, common_paths: Paths):
        """Builds full paths for setup step """
        self.fast5_input_dir = common_paths.data_dir / self.fast5_input_dir_name
        self.reference_genome_dir = common_paths.reference_genome_dir
        self.full_reference_genome_path = self.reference_genome_dir / self.reference_genome_name

class SetupStep(BaseModel):
    params: dict
    downloads: SetupDownloads
    paths: SetupPaths


class BasecallingMethodComplex(BaseModel):
    """Settings for the model complex selection method"""
    model_speed: str
    basecalling_modifications: Optional[str] = None
    kit_name: Optional[str] = None


class BasecallingMethodExplicit(BaseModel):
    base_model_name: str
    mod_model_name: Optional[str] = None


class BasecallingParams(BaseModel):
    # Use EITHER the model complex OR the explicit base/mod model names

    # In the end, this will be where the dorado command building happens. That's for the future though.
    method: Literal["complex", "explicit"]
    complex_settings: BasecallingMethodComplex
    explicit_settings: BasecallingMethodExplicit

    # Other settings
    batch_size: int

    @computed_field
    @property
    def resolved_base_model(self) -> str:
        """Returns the final base model name"""
        if self.method == 'complex':
            if self.complex_settings.basecalling_modifications:
                return f"{self.complex_settings.kit_name}_{self.complex_settings.model_speed}.cfg"
            return None
        else:
            return self.explicit_settings.base_model_name

    @computed_field
    @property
    def resolved_mod_model(self) -> str:
        """Returns the final base model name"""
        if self.method == 'complex':
            if self.complex_settings.basecalling_modifications:
                base_model = f"{self.complex_settings.kit_name}_{self.complex_settings.model_speed}.cfg"
                return f"{base_model}_{self.complex_settings.basecalling_modifications}"
            return None
        else:
            return self.explicit_settings.mod_model_name


class BasecallingPaths(BaseModel):
    # File names
    basecalled_output_dir_name: str = "basecalled_output"
    demultiplexed_dir_name: str = "demultiplexed_output"
    pod5_dir_name: str = "pod5"
    dorado_model_dir_name: str = "models"
    unaligned_bam_name: str
    pod5_name: str

    basecalled_output_dir: Optional[Path] = None
    demultiplexed_output_dir: Optional[Path] = None
    pod5_dir: Optional[Path] = None
    full_pod5_path: Optional[Path] = None
    full_unaligned_bam_path: Optional[Path] = None
    dorado_model_dir: Optional[Path] = None
    def build_paths(self, common_paths: Paths):
        self.basecalled_output_dir = common_paths.data_dir / self.basecalled_output_dir_name
        self.demultiplexed_output_dir = common_paths.data_dir / self.demultiplexed_dir_name
        self.pod5_dir = common_paths.data_dir / self.pod5_dir_name
        self.full_pod5_path = self.pod5_dir / self.pod5_name
        self.full_unaligned_bam_path = self.basecalled_output_dir / self.unaligned_bam_name
        self.dorado_model_dir = common_paths.root / self.dorado_model_dir_name

class BasecallingStep(BaseModel):
    params: BasecallingParams
    paths: BasecallingPaths


class AlignmentPaths(BaseModel):
    # File names
    indexed_ref_fasta_name: str
    aligned_bam_name: str
    alignment_flagstat_name: str
    alignment_stats_name: str
    alignment_output_dir_name: str = "alignment_output"
    alignment_qc_dir_name: str = "alignment_qc"

    ref_genome_subdir: Optional[Path] = ""
    alignment_output_dir: Optional[Path] = None
    qc_dir: Optional[Path] = None
    full_indexed_genome_path: Optional[Path] = None
    full_aligned_bam_path: Optional[Path] = None
    full_flagstat_path: Optional[Path] = None
    full_stats_path: Optional[Path] = None

    def build_paths(self, common_paths: Paths):
        self.alignment_output_dir = common_paths.data_dir / self.alignment_output_dir_name
        self.qc_dir = common_paths.data_dir / self.alignment_qc_dir_name
        self.full_indexed_genome_path = common_paths.reference_genome_dir / self.ref_genome_subdir / self.indexed_ref_fasta_name
        self.full_aligned_bam_path = self.alignment_output_dir / self.aligned_bam_name
        self.full_flagstat_path = self.qc_dir / self.alignment_flagstat_name
        self.full_stats_path = self.qc_dir / self.alignment_stats_name


class AlignmentStep(BaseModel):
    paths: AlignmentPaths


class MethylationPaths(BaseModel):
    # File names
    methylation_bed_name: str
    methylation_log_name: str
    methylation_dir_name: str = "methylation"

    methylation_dir: Optional[Path] = None
    full_bed_path: Optional[Path] = None
    full_meth_log_path: Optional[Path] = None

    def build_paths(self, common_paths: Paths):
        self.methylation_dir = common_paths.data_dir / self.methylation_dir_name
        self.full_bed_path = self.methylation_dir / self.methylation_bed_name
        self.full_meth_log_path = self.methylation_dir / self.methylation_log_name


class MethylationStep(BaseModel):
    paths: MethylationPaths


class AnalysisParams(BaseModel):
    methylation_aggregation_chunksize: int
    deconv_algorithm: str


class AnalysisTools(BaseModel):
    uxm_dir_name: str = "UXM_deconv"
    wgbstools_dir_name: str = "wgbs_tools"
    methatlas_dir_name: str = "meth_atlas"

    uxm_dir: Optional[Path] = None
    wgbstools_dir: Optional[Path] = None
    meth_atlas_dir: Optional[Path] = None
    uxm_exe: Path
    wgbstools_exe: Path
    methatlas_exe: Path

    def build_paths(self, common_paths: Paths):
        self.uxm_dir = common_paths.externals_dir / self.uxm_dir_name
        self.wgbstools_dir = common_paths.externals_dir / self.wgbstools_dir_name
        self.meth_atlas_dir = common_paths.externals_dir / self.methatlas_dir_name


class AnalysisPaths(BaseModel):
    atlas_file_name: str
    manifest_name: str
    deconvolution_results_name: str
    atlas_dir_name: str = "atlas"
    deconvolution_dir_name: str = "deconvolution"
    file_to_deconvolute_name: str = "file_to_deconvolute.csv"

    analysis_dir: Optional[Path] = None
    atlas_dir: Optional[Path] = None
    deconvolution_dir: Optional[Path] = None
    full_path_atlas_file: Optional[Path] = None
    full_path_manifest: Optional[Path] = None
    full_path_deconvolution_results: Optional[Path] = None
    file_to_deconvolute: Optional[Path] = None

    def build_paths(self, common_paths: Paths):
        self.analysis_dir = common_paths.results_dir
        self.atlas_dir = common_paths.data_dir / self.atlas_dir_name
        self.deconvolution_dir = self.analysis_dir / self.deconvolution_dir_name
        self.full_path_atlas_file = self.atlas_dir / self.atlas_file_name
        self.full_path_manifest = self.atlas_dir / self.manifest_name
        self.full_path_deconvolution_results = self.analysis_dir / self.deconvolution_results_name
        self.file_to_deconvolute = self.analysis_dir / self.file_to_deconvolute_name

class AnalysisStep(BaseModel):
    params: AnalysisParams
    tools: AnalysisTools
    paths: AnalysisPaths


class PipelineSteps(BaseModel):
    setup: SetupStep
    basecalling: BasecallingStep
    alignment: AlignmentStep
    methylation: MethylationStep
    analysis: AnalysisStep


class Tools(BaseModel):
    dorado: str


class AppSettings(BaseModel):
    globals: Globals
    pipeline_control: PipelineControl
    pipeline_steps: PipelineSteps
    tools: Tools
    paths: Paths = Paths()  # Default to an empty instance


def deep_merge(d1: Dict[str, Any], d2: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merges 2 dictionaries. Values from d1 will overwrite d2.

    This is used to combine config.yaml and runtime_config.yaml. runtime_config.yaml
    overwrites config.yaml.
    """

    for k, v in d1.items():
        if k in d2 and isinstance(d2[k], dict) and isinstance(v, collections.abc.Mapping):
            deep_merge(v, d2[k])
        else:
            d2[k] = v

    return d2


def load_and_validate_configs(config_path: Path, runtime_config_path: Path) -> AppSettings:
    """
    Loads two YAML config files (runtime_config.yaml and config.yaml), deep-merges them, and validates them with
    Pydantic. The config in the second position will override values from the config in the first position.
    """

    # Load the YAML files into dictionaries
    try:
        with open(config_path, 'r') as f:
            user_config_data = yaml.safe_load(f) or {}
    except FileNotFoundError:
        raise FileNotFoundError(f"Base config file not found: {config_path}")

    try:
        with open(runtime_config_path, 'r') as f:
            runtime_config_data = yaml.safe_load(f) or {}
    except FileNotFoundError:
        raise FileNotFoundError(f"Base config file not found: {config_path}")

    # Do the deep merge (on a copy of the config data so we don't change in-place)
    config_data = deep_merge(runtime_config_data, user_config_data.copy())

    # Validate the final merged dictionary with Pydantic
    try:
        return AppSettings(**config_data)
    except ValidationError as e:
        raise ValueError(f"Configuration error after merging configs:\n{e}")

def print_config(part_to_print: AppSettings):
    config_dict = part_to_print.model_dump(mode='json')
    yaml_str = yaml.dump(config_dict, sort_keys=False, indent=2)
    print(f"\n--- Config {part_to_print} Configuration ---")
    print(yaml_str)