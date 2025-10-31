from __future__ import annotations

import collections

import yaml
from pathlib import Path
from typing import Optional, Any, Dict, Literal
from pydantic import BaseModel, ValidationError, AnyUrl, computed_field

from src.config.validators import validate_path, validate_pod5, resolve_path


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
    fast5_input_dir_name: str = "fast5_input"
    fast5_input_dir: Optional[Path] = None

    def _validate(self):
        pass

    def _build(self, common_paths: Paths):
        """Builds full paths for setup step """
        self.fast5_input_dir = resolve_path(common_paths.data_dir, self.fast5_input_dir_name)

    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()

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
    basecalled_bam_name: str
    pod5_input_path: str
    basecalled_output_dir_name: str = "basecalled_output"
    demultiplexed_dir_name: str = "demultiplexed"
    dorado_model_dir_name: str = "models"

    full_pod5_path: Optional[Path] = None
    full_unaligned_bam_path: Optional[Path] = None
    full_demultiplexed_output_dir: Optional[Path] = None
    full_dorado_model_dir: Optional[Path] = None

    def _validate(self):
        # Pod5 path must exist but can be either a file or a directory
        pod5_param_name = "Setup Input Path ('full_pod5_path')"

        validate_pod5(self.full_pod5_path, pod5_param_name)


    def _build(self, common_paths: Paths):
        root_dir = common_paths.root
        basecalled_dir = resolve_path(common_paths.data_dir, self.basecalled_output_dir_name)
        demux_dir = resolve_path(common_paths.data_dir, self.demultiplexed_dir_name)
        model_dir = resolve_path(root_dir, self.dorado_model_dir_name)

        self.full_pod5_path = resolve_path(root_dir, self.pod5_input_path)
        self.full_unaligned_bam_path = resolve_path(basecalled_dir, self.basecalled_bam_name)
        self.full_demultiplexed_output_dir = resolve_path(root_dir, demux_dir)
        self.full_dorado_model_dir = resolve_path(root_dir, model_dir)

    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()

class BasecallingStep(BaseModel):
    params: BasecallingParams
    paths: BasecallingPaths


class AlignmentPaths(BaseModel):
    reference_genome_dir_name: str = "reference_genomes"
    alignment_dir_name: str = "alignment"
    qc_dir_name: str = "qc"
    genome_id: Optional[str] = None
    custom_fasta_reference: Optional[str] = None
    aligned_bam_name: str
    alignment_flagstat_name: str
    alignment_stats_name: str

    full_ref_fasta_path: Optional[Path] = None
    full_aligned_bam_path: Optional[Path] = None
    full_flagstat_path: Optional[Path] = None
    full_stats_path: Optional[Path] = None
    alignment_output_dir: Optional[Path] = None
    qc_output_dir: Optional[Path] = None

    def _find_reference_fasta(self, common_paths: Paths):
        # Reference fasta depends on what the user defined in the config.yaml file.
        if self.custom_fasta_reference:
            # If the user did specify a custom fasta reference, just use that one.
            print(f"Using custom reference FASTA path: {self.custom_fasta_reference}")
            user_path = self.custom_fasta_reference
            return resolve_path(common_paths.root, user_path)

        elif self.genome_id:
            # If they just provided a genome id (e.g. hg38), look in the expected places.

            # Base reference dir is the reference_genomes directory.
            base_ref_dir = common_paths.root / self.reference_genome_dir_name

            # Define the possible places the genome will be found in.
            search_paths = [
                # Inside a dedicated folder (e.g. reference_genomes/hg38/hg38.fa)
                base_ref_dir / self.genome_id / f"{self.genome_id}.fa",
                base_ref_dir / self.genome_id / f"{self.genome_id}.fa.gz",
                # Inside a dedicated folder in iGenomes style (e.g. reference_genomes/hg38/genome.fa)
                base_ref_dir / self.genome_id / "genome.fa",
                base_ref_dir / self.genome_id / "genome.fa.gz",
                # Directly inside the reference genome folder (e.g. reference_genomes/hg38.fa)
                base_ref_dir / f"{self.genome_id}.fa",
                base_ref_dir / f"{self.genome_id}.fa.gz"
            ]

            return next((path for path in search_paths if path.is_file()), None)

        return None


    def _validate(self):
        if not (self.genome_id or self.custom_fasta_reference):
            raise ValueError("Configuration Error in alignment: Must provide either 'genome_id' or 'custom_fasta_reference'")
        if self.genome_id and self.custom_fasta_reference:
            raise ValueError("Configuration Error in alignment: Don't provide both a 'genome_id' and 'custom_fasta_reference'")

        validate_path(
            self.full_ref_fasta_path,
            must_exist=True,
            must_be_file=True,
            param_name="Reference Genome FASTA"
        )

    def _build(self, common_paths: Paths):
        self.alignment_output_dir = resolve_path(common_paths.data_dir, self.alignment_dir_name)
        self.qc_output_dir = resolve_path(self.alignment_output_dir, self.qc_dir_name)
        self.full_aligned_bam_path = resolve_path(self.alignment_output_dir, self.aligned_bam_name)
        self.full_flagstat_path = resolve_path(self.qc_output_dir, self.alignment_flagstat_name)
        self.full_stats_path = resolve_path(self.qc_output_dir, self.alignment_stats_name)
        self.full_ref_fasta_path = self._find_reference_fasta(common_paths)

    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()

class AlignmentStep(BaseModel):
    paths: AlignmentPaths


class MethylationPaths(BaseModel):
    methylation_output_dir_name: str = "methylation"
    methylation_bed_name: str
    methylation_log_name: str

    methylation_output_dir: Optional[Path] = None
    final_bed_file: Optional[Path] = None
    final_meth_log_file: Optional[Path] = None

    def _validate(self):
        pass
    def _build(self, common_paths: Paths):
        self.methylation_output_dir = resolve_path(common_paths.data_dir, self.methylation_output_dir_name)
        self.final_bed_file = resolve_path(self.methylation_output_dir, self.methylation_bed_name)
        self.final_meth_log_file = resolve_path(self.methylation_output_dir, self.methylation_log_name)
    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()

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

    def _validate(self):
        pass
    def _build(self, common_paths: Paths):
        self.uxm_dir = resolve_path(common_paths.externals_dir, self.uxm_dir_name)
        self.wgbstools_dir = resolve_path(common_paths.externals_dir, self.wgbstools_dir_name)
        self.meth_atlas_dir = resolve_path(common_paths.externals_dir, self.methatlas_dir_name)

    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()


class AnalysisPaths(BaseModel):
    atlas_file_name: str
    manifest_name: str
    deconvolution_results_name: str
    atlas_dir_name: str = "atlas"
    deconvolution_dir_name: str = "deconvolution"

    full_atlas_path: Optional[Path] = None
    full_manifest_path: Optional[Path] = None
    full_deconv_results_path: Optional[Path] = None

    def _validate(self):
        pass

    def _build(self, common_paths: Paths):
        analysis_dir = common_paths.results_dir
        atlas_dir = resolve_path(common_paths.data_dir, self.atlas_dir_name)
        deconvolution_dir = resolve_path(analysis_dir, self.deconvolution_dir_name)

        self.full_atlas_path = resolve_path(atlas_dir, self.atlas_file_name)
        self.full_manifest_path = resolve_path(atlas_dir, self.manifest_name)
        self.full_deconv_results_path = resolve_path(deconvolution_dir, self.deconvolution_results_name)

        print(f"Full atlas path - {self.full_atlas_path}")
        print(f"Full manifest path - {self.full_manifest_path}")
        print(f"Full deconv results path - {self.full_deconv_results_path}")

    def build_and_validate(self, common_paths: Paths):
        self._build(common_paths)
        self._validate()


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
