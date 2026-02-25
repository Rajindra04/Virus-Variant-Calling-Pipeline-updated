"""
Configuration loader for the Virus Variant Calling Pipeline.

Loads a YAML config file and provides validated access to all parameters.
Each virus gets its own config file (e.g., configs/denv1.yaml).
"""

import os
import yaml
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Default values -- used when a key is missing from the config file
DEFAULTS = {
    'virus_name': 'Unknown virus',
    'genome_type': 'Unknown',
    'ploidy': 1,
    'expected_genome_size': 0,
    'database_name': 'custom_db',
    'fastp': {
        'qualified_quality_phred': 20,
        'length_required': 50,
        'cut_front': True,
        'cut_tail': True,
        'cut_window_size': 4,
        'cut_mean_quality': 20,
        'detect_adapter_for_pe': True,
        'correction': True,
        'overlap_len_require': 30,
        'threads': 4,
    },
    'alignment_filtering': {
        'min_mapping_quality': 20,
        'exclude_flags': '0x904',
    },
    'deduplication': {
        'enabled': True,
        'remove_duplicates': True,
    },
    'primer_trimming': {
        'min_quality': 20,
        'min_length': 30,
        'sliding_window': 4,
        'include_reads_no_primer': True,
    },
    'coverage': {
        'min_base_quality': 20,
        'min_mapping_quality': 20,
    },
    'consensus': {
        'mpileup_max_depth': 10000,
        'mpileup_min_base_quality': 20,
        'mpileup_min_mapping_quality': 20,
        'ivar_min_quality': 20,
        'ivar_min_frequency': 0.5,
        'ivar_min_depth': 10,
        'ivar_ambiguous_char': 'N',
    },
    'variant_calling': {
        'ploidy': 1,
        'standard_min_confidence': 30,
        'min_base_quality_score': 20,
    },
    'vcf_filtering': {
        'filters': {
            'LowQD': 'QD < 2.0',
            'StrandBias': 'FS > 60.0',
            'LowMQ': 'MQ < 40.0',
            'LowDepth': 'DP < 10',
        },
        'select_pass_only': True,
    },
    'annotation': {
        'snpeff_memory': '4g',
    },
}


def _deep_merge(defaults, overrides):
    """Recursively merge overrides into defaults."""
    result = defaults.copy()
    for key, value in overrides.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(config_path):
    """Load and validate a virus config YAML file.

    Args:
        config_path: Path to the YAML config file.

    Returns:
        dict with all pipeline parameters, defaults filled in for missing keys.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        user_config = yaml.safe_load(f) or {}

    config = _deep_merge(DEFAULTS, user_config)

    # Ensure variant_calling.ploidy matches top-level ploidy
    config['variant_calling']['ploidy'] = config['ploidy']

    logging.info(f"Loaded config for: {config['virus_name']} "
                 f"(ploidy={config['ploidy']}, db={config['database_name']})")

    return config
