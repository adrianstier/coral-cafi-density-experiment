"""
Configuration for CAFI MRB Analysis Agents
==========================================
Central configuration file for all agents in the coral reef analysis pipeline.

Author: CAFI Analysis Team
Created: 2025-01-10
"""

from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import os

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Project root (parent of agents directory)
PROJECT_ROOT = Path(__file__).parent.parent.resolve()

# Key directories
DATA_DIR = PROJECT_ROOT / "data"
SCRIPTS_DIR = PROJECT_ROOT / "scripts" / "MRB"
OUTPUT_DIR = PROJECT_ROOT / "output" / "MRB"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

# Data subdirectories
RAW_DATA_DIR = DATA_DIR / "MRB Amount"

# Output subdirectories
FIGURES_DIR = OUTPUT_DIR / "figures"
TABLES_DIR = OUTPUT_DIR / "tables"
PUBLICATION_FIGURES_DIR = FIGURES_DIR / "publication-figures"


# =============================================================================
# SCRIPT CONFIGURATION
# =============================================================================

@dataclass
class ScriptConfig:
    """Configuration for an R analysis script."""
    name: str
    path: Path
    description: str
    dependencies: List[str] = field(default_factory=list)
    outputs: List[str] = field(default_factory=list)
    estimated_runtime_seconds: int = 30
    required_data_files: List[str] = field(default_factory=list)


# Define the analysis pipeline scripts in execution order
PIPELINE_SCRIPTS: Dict[str, ScriptConfig] = {
    "1.libraries": ScriptConfig(
        name="1.libraries",
        path=SCRIPTS_DIR / "1.libraries.R",
        description="Load all required R packages",
        dependencies=[],
        outputs=[],
        estimated_runtime_seconds=5
    ),
    "utils": ScriptConfig(
        name="utils",
        path=SCRIPTS_DIR / "utils.R",
        description="Utility functions for data loading and analysis",
        dependencies=["1.libraries"],
        outputs=[],
        estimated_runtime_seconds=2
    ),
    "mrb_figure_standards": ScriptConfig(
        name="mrb_figure_standards",
        path=SCRIPTS_DIR / "mrb_figure_standards.R",
        description="Publication figure styling standards",
        dependencies=["1.libraries"],
        outputs=[],
        estimated_runtime_seconds=2
    ),
    "3.abundance": ScriptConfig(
        name="3.abundance",
        path=SCRIPTS_DIR / "3.abundance.R",
        description="CAFI abundance analysis with negative binomial GLMs",
        dependencies=["1.libraries", "utils"],
        outputs=["abundance_figures", "abundance_tables"],
        estimated_runtime_seconds=60,
        required_data_files=["1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"]
    ),
    "4d.diversity": ScriptConfig(
        name="4d.diversity",
        path=SCRIPTS_DIR / "4d.diversity.R",
        description="Diversity metrics and rarefaction analysis",
        dependencies=["1.libraries", "utils", "3.abundance"],
        outputs=["diversity_figures", "rarefaction_curves"],
        estimated_runtime_seconds=90
    ),
    "5.fishes": ScriptConfig(
        name="5.fishes",
        path=SCRIPTS_DIR / "5.fishes.R",
        description="Fish community analysis",
        dependencies=["1.libraries", "utils"],
        outputs=["fish_figures", "fish_tables"],
        estimated_runtime_seconds=45
    ),
    "6.coral-growth": ScriptConfig(
        name="6.coral-growth",
        path=SCRIPTS_DIR / "6.coral-growth.R",
        description="Allometric coral growth analysis from photogrammetry",
        dependencies=["1.libraries", "utils"],
        outputs=["coral_growth.csv", "growth_figures"],
        estimated_runtime_seconds=60,
        required_data_files=[
            "MRB_2019_200K_mesh_measure.csv",
            "MRB_May_2021_200K_mesh_measure.csv"
        ]
    ),
    "7.coral-physiology": ScriptConfig(
        name="7.coral-physiology",
        path=SCRIPTS_DIR / "7.coral-physiology.R",
        description="Coral physiology integration and performance PCA",
        dependencies=["1.libraries", "utils", "6.coral-growth"],
        outputs=["physiology_figures", "performance_metrics"],
        estimated_runtime_seconds=45,
        required_data_files=["1. amount_master_phys_data_v5.csv"]
    ),
    "8.coral-caffi": ScriptConfig(
        name="8.coral-caffi",
        path=SCRIPTS_DIR / "8.coral-caffi.R",
        description="CAFI-coral community relationships analysis",
        dependencies=["1.libraries", "utils", "6.coral-growth", "7.coral-physiology"],
        outputs=["cafi_coral_figures", "species_lmm_results"],
        estimated_runtime_seconds=120
    ),
    "12.nmds_permanova_cafi": ScriptConfig(
        name="12.nmds_permanova_cafi",
        path=SCRIPTS_DIR / "12.nmds_permanova_cafi.R",
        description="NMDS ordination and PERMANOVA community composition",
        dependencies=["1.libraries", "utils"],
        outputs=["nmds_figures", "permanova_tables"],
        estimated_runtime_seconds=60
    ),
    "14.compile-manuscript-statistics": ScriptConfig(
        name="14.compile-manuscript-statistics",
        path=SCRIPTS_DIR / "14.compile-manuscript-statistics.R",
        description="Compile all statistics for manuscript",
        dependencies=["1.libraries", "utils"],
        outputs=["MANUSCRIPT_STATS_TABLE.csv", "statistics_html"],
        estimated_runtime_seconds=30
    ),
}

# Execution order for full pipeline
PIPELINE_ORDER = [
    "1.libraries",
    "utils",
    "mrb_figure_standards",
    "3.abundance",
    "4d.diversity",
    "5.fishes",
    "6.coral-growth",
    "7.coral-physiology",
    "8.coral-caffi",
    "12.nmds_permanova_cafi",
    "14.compile-manuscript-statistics",
]


# =============================================================================
# DATA FILE CONFIGURATION
# =============================================================================

@dataclass
class DataFileConfig:
    """Configuration for a data file."""
    name: str
    path: Path
    description: str
    file_type: str  # csv, xlsx, rds
    required: bool = True
    columns: List[str] = field(default_factory=list)


REQUIRED_DATA_FILES: Dict[str, DataFileConfig] = {
    "cafi_community": DataFileConfig(
        name="cafi_community",
        path=RAW_DATA_DIR / "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv",
        description="Main CAFI community abundance data",
        file_type="csv",
        required=True,
        columns=["coral_id", "site", "species", "count"]
    ),
    "treatment": DataFileConfig(
        name="treatment",
        path=RAW_DATA_DIR / "coral_id_position_treatment.csv",
        description="Coral treatment assignments",
        file_type="csv",
        required=True,
        columns=["coral_id", "position", "treatment"]
    ),
    "physiology": DataFileConfig(
        name="physiology",
        path=RAW_DATA_DIR / "1. amount_master_phys_data_v5.csv",
        description="Coral physiology measurements",
        file_type="csv",
        required=True,
        columns=["coral_id", "carbohydrate", "protein", "zooxanthellae", "afdw"]
    ),
    "growth_2019": DataFileConfig(
        name="growth_2019",
        path=RAW_DATA_DIR / "MRB_2019_200K_mesh_measure.csv",
        description="Initial 3D mesh measurements (2019)",
        file_type="csv",
        required=True
    ),
    "growth_2021": DataFileConfig(
        name="growth_2021",
        path=RAW_DATA_DIR / "MRB_May_2021_200K_mesh_measure.csv",
        description="Final 3D mesh measurements (2021)",
        file_type="csv",
        required=True
    ),
}


# =============================================================================
# STATISTICAL ANALYSIS CONFIGURATION
# =============================================================================

@dataclass
class StatisticalTestConfig:
    """Configuration for statistical tests."""
    name: str
    test_type: str
    variables: List[str]
    script_source: str
    description: str


KEY_STATISTICAL_TESTS: List[StatisticalTestConfig] = [
    StatisticalTestConfig(
        name="treatment_growth_ancova",
        test_type="ANCOVA (LMM)",
        variables=["treatment", "initial_volume", "final_volume"],
        script_source="6.coral-growth.R",
        description="Treatment effect on size-corrected coral growth"
    ),
    StatisticalTestConfig(
        name="treatment_carbohydrate",
        test_type="LMM",
        variables=["treatment", "carbohydrate_mg_cm2"],
        script_source="7.coral-physiology.R",
        description="Treatment effect on coral carbohydrate content"
    ),
    StatisticalTestConfig(
        name="treatment_performance_pc1",
        test_type="LMM",
        variables=["treatment", "performance_PC1"],
        script_source="7.coral-physiology.R",
        description="Treatment effect on integrated coral performance"
    ),
    StatisticalTestConfig(
        name="cafi_coral_relationship",
        test_type="LMM",
        variables=["cafi_PC1", "coral_performance_PC1"],
        script_source="8.coral-caffi.R",
        description="CAFI community composition predicting coral performance"
    ),
    StatisticalTestConfig(
        name="community_composition_permanova",
        test_type="PERMANOVA",
        variables=["treatment", "cafi_community_matrix"],
        script_source="12.nmds_permanova_cafi.R",
        description="Treatment effect on CAFI community composition"
    ),
]


# =============================================================================
# FIGURE CONFIGURATION
# =============================================================================

@dataclass
class FigureConfig:
    """Configuration for publication figures."""
    figure_number: int
    name: str
    script_source: str
    output_path: Path
    width_inches: float = 8
    height_inches: float = 6
    dpi: int = 300
    description: str = ""


PUBLICATION_FIGURES: List[FigureConfig] = [
    FigureConfig(
        figure_number=1,
        name="community_abundance_richness",
        script_source="3.abundance.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure1_community_abundance.png",
        width_inches=10,
        height_inches=8,
        description="CAFI community abundance and richness by treatment"
    ),
    FigureConfig(
        figure_number=2,
        name="species_abundance_scaling",
        script_source="3.abundance.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure2_species_scaling.png",
        description="Species abundance scaling relationships"
    ),
    FigureConfig(
        figure_number=3,
        name="nmds_ordination",
        script_source="12.nmds_permanova_cafi.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure3_nmds.png",
        description="NMDS ordination of CAFI community composition"
    ),
    FigureConfig(
        figure_number=4,
        name="cafi_pca",
        script_source="8.coral-caffi.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure4_cafi_pca.png",
        description="PCA of CAFI community composition"
    ),
    FigureConfig(
        figure_number=5,
        name="coral_physiology_treatment",
        script_source="7.coral-physiology.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure5_physiology.png",
        description="Coral physiology metrics by treatment"
    ),
    FigureConfig(
        figure_number=6,
        name="coral_condition_pc1",
        script_source="7.coral-physiology.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure6_condition_pc1.png",
        description="Coral integrated performance (PC1) by treatment"
    ),
    FigureConfig(
        figure_number=7,
        name="species_coral_relationships",
        script_source="8.coral-caffi.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure7_species_relationships.png",
        description="Individual species relationships with coral condition"
    ),
    FigureConfig(
        figure_number=8,
        name="performance_correlation_heatmap",
        script_source="8.coral-caffi.R",
        output_path=PUBLICATION_FIGURES_DIR / "Figure8_correlation_heatmap.png",
        description="Performance correlation heatmap"
    ),
]


# =============================================================================
# AGENT CONFIGURATION
# =============================================================================

@dataclass
class AgentConfig:
    """Configuration for an agent."""
    name: str
    description: str
    capabilities: List[str]
    default_timeout_seconds: int = 300
    max_retries: int = 3


AGENTS: Dict[str, AgentConfig] = {
    "pipeline": AgentConfig(
        name="Pipeline Runner",
        description="Orchestrates R script execution and manages dependencies",
        capabilities=["run_script", "run_pipeline", "check_dependencies", "validate_outputs"],
        default_timeout_seconds=600
    ),
    "stats": AgentConfig(
        name="Statistical Analyst",
        description="Performs statistical analysis and interprets model results",
        capabilities=["run_statistical_test", "extract_results", "format_for_manuscript", "check_assumptions"],
        default_timeout_seconds=300
    ),
    "figure": AgentConfig(
        name="Figure Generator",
        description="Creates and updates publication-quality figures",
        capabilities=["generate_figure", "update_styling", "batch_export", "check_quality"],
        default_timeout_seconds=180
    ),
    "validation": AgentConfig(
        name="Data Validator",
        description="Validates data integrity, checks for missing values, and ensures reproducibility",
        capabilities=["validate_data", "check_schema", "compare_outputs", "generate_report"],
        default_timeout_seconds=120
    ),
}


# =============================================================================
# ENVIRONMENT CONFIGURATION
# =============================================================================

# R executable path (can be overridden by environment variable)
R_EXECUTABLE = os.environ.get("R_EXECUTABLE", "Rscript")

# Logging configuration
LOG_LEVEL = os.environ.get("CAFI_LOG_LEVEL", "INFO")
LOG_DIR = OUTPUT_DIR / "logs"

# Ensure log directory exists
LOG_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_script_config(script_name: str) -> Optional[ScriptConfig]:
    """Get configuration for a specific script."""
    return PIPELINE_SCRIPTS.get(script_name)


def get_data_file_config(file_name: str) -> Optional[DataFileConfig]:
    """Get configuration for a specific data file."""
    return REQUIRED_DATA_FILES.get(file_name)


def get_figure_config(figure_number: int) -> Optional[FigureConfig]:
    """Get configuration for a specific figure by number."""
    for fig in PUBLICATION_FIGURES:
        if fig.figure_number == figure_number:
            return fig
    return None


def validate_environment() -> Dict[str, bool]:
    """Validate that all required paths and files exist."""
    results = {
        "project_root": PROJECT_ROOT.exists(),
        "data_dir": DATA_DIR.exists(),
        "scripts_dir": SCRIPTS_DIR.exists(),
        "output_dir": OUTPUT_DIR.exists(),
        "r_executable": os.system(f"which {R_EXECUTABLE} > /dev/null 2>&1") == 0,
    }

    # Check required data files
    for name, config in REQUIRED_DATA_FILES.items():
        results[f"data_{name}"] = config.path.exists()

    # Check scripts
    for name, config in PIPELINE_SCRIPTS.items():
        results[f"script_{name}"] = config.path.exists()

    return results


if __name__ == "__main__":
    # Print configuration summary when run directly
    print("=" * 60)
    print("CAFI MRB Analysis Agent Configuration")
    print("=" * 60)
    print(f"\nProject Root: {PROJECT_ROOT}")
    print(f"Data Directory: {DATA_DIR}")
    print(f"Scripts Directory: {SCRIPTS_DIR}")
    print(f"Output Directory: {OUTPUT_DIR}")
    print(f"\nPipeline Scripts: {len(PIPELINE_SCRIPTS)}")
    print(f"Required Data Files: {len(REQUIRED_DATA_FILES)}")
    print(f"Publication Figures: {len(PUBLICATION_FIGURES)}")
    print(f"Statistical Tests: {len(KEY_STATISTICAL_TESTS)}")

    print("\n" + "=" * 60)
    print("Environment Validation")
    print("=" * 60)
    validation = validate_environment()
    for check, passed in validation.items():
        status = "OK" if passed else "MISSING"
        print(f"  {check}: {status}")
