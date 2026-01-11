"""
Main Agent Orchestrator for CAFI MRB Analysis
==============================================
Central orchestrator that coordinates all specialized agents
for the coral reef analysis pipeline.

Author: CAFI Analysis Team
Created: 2025-01-10
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    PROJECT_ROOT, OUTPUT_DIR, PIPELINE_ORDER,
    validate_environment
)
from utils import setup_logger, write_json_report
from pipeline_agent import PipelineAgent, PipelineRun, PipelineStatus
from stats_agent import StatsAgent
from figure_agent import FigureAgent
from validation_agent import ValidationAgent, ValidationReport


# =============================================================================
# ORCHESTRATOR
# =============================================================================

class Orchestrator:
    """
    Central orchestrator for all CAFI MRB analysis agents.

    Coordinates:
    - PipelineAgent: R script execution and dependency management
    - StatsAgent: Statistical analysis and result interpretation
    - FigureAgent: Publication figure generation and quality control
    - ValidationAgent: Data validation and quality assurance

    Provides high-level workflows for common analysis tasks.
    """

    def __init__(self, log_file: Optional[Path] = None):
        """Initialize the orchestrator with all agents."""
        self.logger = setup_logger(
            "Orchestrator",
            log_file or (OUTPUT_DIR / "logs" / "orchestrator.log")
        )

        # Initialize all agents
        self.pipeline = PipelineAgent()
        self.stats = StatsAgent()
        self.figures = FigureAgent()
        self.validation = ValidationAgent()

        self.logger.info("Orchestrator initialized with all agents")

    # -------------------------------------------------------------------------
    # ENVIRONMENT & STATUS
    # -------------------------------------------------------------------------

    def check_environment(self) -> Dict[str, Any]:
        """
        Check the analysis environment.

        Returns:
            Dictionary with environment status
        """
        validation = validate_environment()

        status = {
            "timestamp": datetime.now().isoformat(),
            "project_root": str(PROJECT_ROOT),
            "environment_checks": validation,
            "all_ok": all(validation.values()),
            "missing": [k for k, v in validation.items() if not v]
        }

        if status["all_ok"]:
            self.logger.info("Environment check: All OK")
        else:
            self.logger.warning(f"Environment check: Missing items: {status['missing']}")

        return status

    def get_status(self) -> Dict[str, Any]:
        """
        Get overall analysis status.

        Returns:
            Comprehensive status dictionary
        """
        return {
            "timestamp": datetime.now().isoformat(),
            "environment": self.check_environment(),
            "pipeline": self.pipeline.get_pipeline_status(),
            "figures": {
                "statuses": self.figures.get_all_figure_status(),
                "total": 8,
                "generated": sum(1 for s in self.figures.get_all_figure_status() if s.get("exists"))
            },
            "agents": {
                "pipeline": "ready",
                "stats": "ready",
                "figures": "ready",
                "validation": "ready"
            }
        }

    # -------------------------------------------------------------------------
    # HIGH-LEVEL WORKFLOWS
    # -------------------------------------------------------------------------

    def run_full_analysis(
        self,
        validate_first: bool = True,
        generate_figures: bool = True,
        compile_stats: bool = True
    ) -> Dict[str, Any]:
        """
        Run the complete analysis workflow.

        Args:
            validate_first: Run data validation before analysis
            generate_figures: Generate publication figures
            compile_stats: Compile manuscript statistics

        Returns:
            Dictionary with results from all steps
        """
        results = {
            "timestamp": datetime.now().isoformat(),
            "workflow": "full_analysis",
            "steps": {}
        }

        self.logger.info("=" * 60)
        self.logger.info("STARTING FULL ANALYSIS WORKFLOW")
        self.logger.info("=" * 60)

        # Step 1: Validate data
        if validate_first:
            self.logger.info("Step 1: Data Validation")
            validation_report = self.validation.run_all_validations()
            results["steps"]["validation"] = {
                "status": "completed",
                "passed": validation_report.passed,
                "warnings": validation_report.warnings,
                "failed": validation_report.failed,
                "quality_score": validation_report.data_quality_score
            }

            if validation_report.failed > 0:
                self.logger.warning(f"Validation found {validation_report.failed} failures")
                # Continue anyway but note the issues

        # Step 2: Run pipeline
        self.logger.info("Step 2: Running R Pipeline")
        pipeline_run = self.pipeline.run_pipeline(stop_on_error=False)
        results["steps"]["pipeline"] = {
            "status": pipeline_run.status.value,
            "succeeded": pipeline_run.scripts_succeeded,
            "failed": pipeline_run.scripts_failed,
            "duration_seconds": pipeline_run.total_duration_seconds
        }

        # Step 3: Compile statistics
        if compile_stats:
            self.logger.info("Step 3: Compiling Statistics")
            stats_results = self.stats.extract_key_results()
            results["steps"]["statistics"] = {
                "status": "completed",
                "results_count": len(stats_results)
            }

            # Export statistics
            self.stats.export_results_csv()
            self.stats.export_results_markdown()

        # Step 4: Check/generate figures
        if generate_figures:
            self.logger.info("Step 4: Figure Generation")
            figure_statuses = self.figures.get_all_figure_status()
            missing_figures = [s for s in figure_statuses if not s.get("exists")]

            results["steps"]["figures"] = {
                "status": "completed",
                "total": 8,
                "existing": 8 - len(missing_figures),
                "missing": len(missing_figures)
            }

        # Generate final report
        self.logger.info("Generating Final Report")
        report_path = OUTPUT_DIR / "logs" / f"full_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        write_json_report(results, report_path)

        self.logger.info("=" * 60)
        self.logger.info("ANALYSIS WORKFLOW COMPLETE")
        self.logger.info(f"Report saved to: {report_path}")
        self.logger.info("=" * 60)

        return results

    def run_quick_analysis(self) -> Dict[str, Any]:
        """
        Run a quick analysis focusing on key results.

        Runs only the core scripts and extracts key statistics.
        """
        results = {
            "timestamp": datetime.now().isoformat(),
            "workflow": "quick_analysis",
            "steps": {}
        }

        self.logger.info("Running quick analysis workflow")

        # Run core scripts only
        core_scripts = [
            "6.coral-growth",
            "7.coral-physiology",
            "14.compile-manuscript-statistics"
        ]

        pipeline_run = self.pipeline.run_pipeline(scripts=core_scripts)
        results["steps"]["pipeline"] = {
            "status": pipeline_run.status.value,
            "scripts_run": len(core_scripts)
        }

        # Extract key results
        treatment_effects = self.stats.get_treatment_effects()
        results["steps"]["key_results"] = {
            "treatment_effects": len(treatment_effects),
            "summary": [
                {
                    "name": r.effect_name,
                    "p_value": r.p_value,
                    "significant": r.p_value < 0.05
                }
                for r in treatment_effects
            ]
        }

        return results

    def validate_and_report(self) -> ValidationReport:
        """
        Run comprehensive validation and generate reports.

        Returns:
            ValidationReport
        """
        self.logger.info("Running comprehensive data validation")

        report = self.validation.run_all_validations()

        # Generate both report formats
        self.validation.generate_report(report)
        self.validation.generate_markdown_report(report)

        return report

    def generate_manuscript_materials(self) -> Dict[str, Path]:
        """
        Generate all materials needed for manuscript.

        Returns:
            Dictionary mapping material type to file paths
        """
        materials = {}

        self.logger.info("Generating manuscript materials")

        # 1. Statistics summary
        stats_csv = self.stats.export_results_csv()
        stats_md = self.stats.export_results_markdown()
        materials["statistics_csv"] = stats_csv
        materials["statistics_markdown"] = stats_md

        # 2. Figure archive
        figure_archive = self.figures.create_figure_archive()
        materials["figures_archive"] = figure_archive

        # 3. Validation report
        validation_report = self.validation.run_all_validations()
        validation_path = self.validation.generate_report(validation_report)
        materials["validation_report"] = validation_path

        # 4. Figure quality check
        quality_results = self.figures.batch_quality_check()
        quality_path = OUTPUT_DIR / "figure_quality_check.json"
        write_json_report(quality_results, quality_path)
        materials["figure_quality"] = quality_path

        self.logger.info(f"Generated {len(materials)} manuscript materials")
        return materials

    # -------------------------------------------------------------------------
    # INTERACTIVE METHODS
    # -------------------------------------------------------------------------

    def query(self, question: str) -> str:
        """
        Answer a question about the analysis.

        Args:
            question: Natural language question

        Returns:
            Answer string
        """
        question_lower = question.lower()

        # Route to appropriate agent based on keywords
        if any(kw in question_lower for kw in ["p-value", "significant", "statistic", "effect"]):
            results = self.stats.extract_key_results()
            summary = self.stats.generate_results_summary()
            return summary

        elif any(kw in question_lower for kw in ["figure", "plot", "image"]):
            statuses = self.figures.get_all_figure_status()
            lines = ["Figure Status:"]
            for s in statuses:
                status = "exists" if s.get("exists") else "missing"
                lines.append(f"  Figure {s['figure_number']}: {status} - {s['name']}")
            return "\n".join(lines)

        elif any(kw in question_lower for kw in ["data", "validation", "quality", "missing"]):
            report = self.validation.run_all_validations()
            return f"Data Quality Score: {report.data_quality_score:.1f}/100\n" \
                   f"Checks Passed: {report.passed}/{report.total_checks}\n" \
                   f"Issues Found: {len(report.issues)}"

        elif any(kw in question_lower for kw in ["script", "pipeline", "run"]):
            return f"Pipeline Scripts: {len(PIPELINE_ORDER)}\n" \
                   f"Order: {', '.join(PIPELINE_ORDER[:5])}..."

        else:
            return "I can answer questions about:\n" \
                   "- Statistics (p-values, effects, significance)\n" \
                   "- Figures (status, generation)\n" \
                   "- Data quality (validation, missing values)\n" \
                   "- Pipeline (scripts, execution)"


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    """Command-line interface for the orchestrator."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CAFI MRB Analysis Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python orchestrator.py --status              Show overall status
  python orchestrator.py --run-full            Run full analysis workflow
  python orchestrator.py --run-quick           Run quick analysis
  python orchestrator.py --validate            Run data validation
  python orchestrator.py --manuscript          Generate manuscript materials
  python orchestrator.py --query "p-values"    Query the analysis

Workflows:
  Full Analysis: Validation -> Pipeline -> Statistics -> Figures
  Quick Analysis: Core scripts only (growth, physiology, stats)
  Manuscript: Generate all materials for publication
        """
    )

    parser.add_argument(
        "--status", action="store_true",
        help="Show overall analysis status"
    )
    parser.add_argument(
        "--check-env", action="store_true",
        help="Check the analysis environment"
    )
    parser.add_argument(
        "--run-full", action="store_true",
        help="Run full analysis workflow"
    )
    parser.add_argument(
        "--run-quick", action="store_true",
        help="Run quick analysis workflow"
    )
    parser.add_argument(
        "--validate", action="store_true",
        help="Run data validation"
    )
    parser.add_argument(
        "--manuscript", action="store_true",
        help="Generate manuscript materials"
    )
    parser.add_argument(
        "--query", type=str,
        help="Query the analysis (natural language)"
    )

    args = parser.parse_args()

    orchestrator = Orchestrator()

    if args.status:
        status = orchestrator.get_status()
        print("\n" + "=" * 60)
        print("CAFI MRB Analysis Status")
        print("=" * 60)
        print(f"\nEnvironment: {'OK' if status['environment']['all_ok'] else 'ISSUES FOUND'}")
        if status['environment'].get('missing'):
            print(f"  Missing: {status['environment']['missing']}")
        print(f"\nFigures: {status['figures']['generated']}/{status['figures']['total']} generated")
        print(f"\nAgents: All ready")
        return

    if args.check_env:
        env = orchestrator.check_environment()
        print("\nEnvironment Check:")
        print("=" * 50)
        for check, passed in env["environment_checks"].items():
            status = "OK" if passed else "MISSING"
            print(f"  [{status}] {check}")
        return

    if args.run_full:
        results = orchestrator.run_full_analysis()
        print("\nFull Analysis Complete!")
        print("=" * 50)
        for step, info in results["steps"].items():
            print(f"  {step}: {info.get('status', 'completed')}")
        return

    if args.run_quick:
        results = orchestrator.run_quick_analysis()
        print("\nQuick Analysis Complete!")
        print("=" * 50)
        print(f"  Pipeline: {results['steps']['pipeline']['status']}")
        print(f"  Key Results: {results['steps']['key_results']['treatment_effects']} treatment effects")
        return

    if args.validate:
        report = orchestrator.validate_and_report()
        print("\nValidation Complete!")
        print("=" * 50)
        print(f"  Quality Score: {report.data_quality_score:.1f}/100")
        print(f"  Passed: {report.passed}")
        print(f"  Warnings: {report.warnings}")
        print(f"  Failed: {report.failed}")
        return

    if args.manuscript:
        materials = orchestrator.generate_manuscript_materials()
        print("\nManuscript Materials Generated:")
        print("=" * 50)
        for name, path in materials.items():
            print(f"  {name}: {path}")
        return

    if args.query:
        answer = orchestrator.query(args.query)
        print(f"\n{answer}")
        return

    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
