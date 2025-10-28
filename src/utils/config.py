from __future__ import annotations

import collections

import yaml
from pathlib import Path
from typing import Optional, Any, Dict, List, Literal
from pydantic import BaseModel, ValidationError, AnyUrl, model_validator


class GeneralParameters(BaseModel):
    threads: int = 4
    sort_memory_limit: str = "2G"


class SetupParameters(BaseModel):
    dorado_version: str
    num_fast5_files: int


class BasecallingMethodComplex(BaseModel):
    """Settings for the model complex selection method"""
    model_speed: str
    basecalling_modifications: Optional[str] = None
    kit_name: Optional[str] = None


class BasecallingMethodExplicit(BaseModel):
    base_model_name: str
    mod_model_name: Optional[str] = None


class BasecallingParameters(BaseModel):
    # Use EITHER the model complex OR the explicit base/mod model names
    method: Literal["complex", "explicit"]
    complex_settings: BasecallingMethodComplex
    explicit_settings: BasecallingMethodExplicit

    model_speed: Optional[str] = None
    kit_name: Optional[str] = None
    basecalling_modifications: Optional[str] = None

    base_model_name: Optional[str] = None
    mod_model_name: Optional[str] = None

    # Other settings
    batch_size: int


class AnalysisParameters(BaseModel):
    methylation_aggregation_chunksize: int
    deconv_algorithm: str


class Parameters(BaseModel):
    """ Container for all tool parameters """
    general: GeneralParameters
    setup: SetupParameters
    basecalling: BasecallingParameters
    analysis: AnalysisParameters


# ==================================================
# Pydantic tools for 'paths' section of config.yaml
# ==================================================

class Submodules(BaseModel):
    uxm_dir: Path
    wgbstools_dir: Path
    meth_atlas_dir: Path


class Paths(BaseModel):
    """ Container for all file and directory paths. """
    core_directories: List[Path]
    submodules: Submodules

    # External download urls
    fast5_download_url: AnyUrl
    reference_genome_url: AnyUrl
    uxm_atlas_url: AnyUrl
    manifest_url: AnyUrl
    full_atlas_url: AnyUrl

    # Raw data
    fast5_input_dir: Path
    pod5_dir: Path
    pod5_name: str

    # Step 1: Basecalling output
    basecalled_output_dir: Path
    demultiplexed_output_dir: Path

    # Step 2: Alignment output
    alignment_output_dir: Path
    qc_dir: Path
    indexed_ref_gen_fasta_name: str
    unaligned_bam_name: str
    aligned_bam_name: str
    alignment_flagstat_name: str
    alignment_stats_name: str
    reference_genome_dir: Path

    # Step 3: Methylation summary outputs
    methylation_dir: Path
    methylation_bed_name: str
    methylation_log_file: str

    # Step 4: Analysis and Deconvolution outputs
    analysis_dir: Path
    atlas_dir: Path
    deconvolution_dir: Path
    file_for_deconvolution_ilmn: str
    file_for_deconvolution_gc: str
    file_for_deconvolution_uxm: str
    file_to_deconvolute: str
    atlas_file_ilmn: str
    atlas_file_gc: str
    atlas_file_uxm: str
    atlas_file_name: str
    illumina_manifest: str
    processed_illumina_manifest: str
    deconvolution_results: str
    uxm_atlas_name: str


class AppSettings(BaseModel):
    parameters: Parameters
    paths: Paths


def load_config_from_yaml(config_file: Path) -> AppSettings:
    """
    Loads, validates and returns the application settings from a YAML file.
    """
    if not config_file.is_file():
        raise FileNotFoundError(f"Configuration file not found at {config_file}")

    try:
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        if not config_data:
            raise ValueError("Config file is empty or invalid.")
        return AppSettings(**config_data)
    except ValidationError as e:
        raise ValueError(f"Configuration error in '{config_file}':\n{e}")
    except Exception as e:
        raise ValueError(f"Failed to load or parse config '{config_file}': {e}")


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
