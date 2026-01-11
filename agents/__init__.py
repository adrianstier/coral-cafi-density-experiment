"""
CAFI MRB Analysis Agents
========================
A suite of Python agents for automating and orchestrating
the coral reef analysis pipeline.

Agents:
- Orchestrator: Central coordinator for all agents
- PipelineAgent: R script execution and dependency management
- StatsAgent: Statistical analysis and interpretation
- FigureAgent: Figure generation and quality control
- ValidationAgent: Data validation and quality assurance

Quick Start:
    from agents import Orchestrator

    orch = Orchestrator()
    status = orch.get_status()
    results = orch.run_full_analysis()

CLI Usage:
    python -m agents --help
    python -m agents --status
    python -m agents --run-full
"""

from .config import (
    PROJECT_ROOT,
    DATA_DIR,
    SCRIPTS_DIR,
    OUTPUT_DIR,
    PIPELINE_SCRIPTS,
    PIPELINE_ORDER,
)

from .orchestrator import Orchestrator
from .pipeline_agent import PipelineAgent
from .stats_agent import StatsAgent
from .figure_agent import FigureAgent
from .validation_agent import ValidationAgent

__version__ = "1.0.0"
__author__ = "CAFI Analysis Team"

__all__ = [
    "Orchestrator",
    "PipelineAgent",
    "StatsAgent",
    "FigureAgent",
    "ValidationAgent",
    "PROJECT_ROOT",
    "DATA_DIR",
    "SCRIPTS_DIR",
    "OUTPUT_DIR",
    "PIPELINE_SCRIPTS",
    "PIPELINE_ORDER",
]
