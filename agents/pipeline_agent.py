"""
Pipeline Runner Agent for CAFI MRB Analysis
============================================
Orchestrates R script execution, manages dependencies, validates outputs,
and provides intelligent pipeline execution.

Author: CAFI Analysis Team
Created: 2025-01-10
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import json

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    PROJECT_ROOT, SCRIPTS_DIR, OUTPUT_DIR, PROCESSED_DATA_DIR,
    PIPELINE_SCRIPTS, PIPELINE_ORDER, ScriptConfig,
    get_script_config
)
from utils import (
    setup_logger, run_r_script, RScriptResult,
    write_json_report, get_file_info, ProgressTracker
)


# =============================================================================
# ENUMS AND DATA CLASSES
# =============================================================================

class PipelineStatus(Enum):
    """Status of a pipeline execution."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class ScriptExecution:
    """Record of a script execution."""
    script_name: str
    status: PipelineStatus
    result: Optional[RScriptResult] = None
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    skipped_reason: Optional[str] = None


@dataclass
class PipelineRun:
    """Record of a full pipeline run."""
    run_id: str
    start_time: datetime
    end_time: Optional[datetime] = None
    scripts: List[ScriptExecution] = field(default_factory=list)
    status: PipelineStatus = PipelineStatus.PENDING
    total_duration_seconds: float = 0
    scripts_succeeded: int = 0
    scripts_failed: int = 0
    scripts_skipped: int = 0


# =============================================================================
# PIPELINE AGENT
# =============================================================================

class PipelineAgent:
    """
    Agent for orchestrating the R analysis pipeline.

    Capabilities:
    - Run individual scripts or full pipeline
    - Manage script dependencies
    - Validate outputs
    - Generate execution reports
    - Handle errors and retries
    """

    def __init__(self, log_file: Optional[Path] = None):
        """Initialize the pipeline agent."""
        self.logger = setup_logger(
            "PipelineAgent",
            log_file or (OUTPUT_DIR / "logs" / "pipeline_agent.log")
        )
        self.current_run: Optional[PipelineRun] = None
        self.execution_history: List[PipelineRun] = []

    # -------------------------------------------------------------------------
    # DEPENDENCY MANAGEMENT
    # -------------------------------------------------------------------------

    def get_dependency_order(self, scripts: List[str]) -> List[str]:
        """
        Get execution order respecting dependencies.

        Args:
            scripts: List of script names to execute

        Returns:
            Ordered list of scripts to execute
        """
        # Build dependency graph
        all_deps = set()

        def collect_deps(script_name: str):
            config = get_script_config(script_name)
            if config:
                for dep in config.dependencies:
                    if dep not in all_deps:
                        all_deps.add(dep)
                        collect_deps(dep)

        for script in scripts:
            all_deps.add(script)
            collect_deps(script)

        # Order by pipeline order
        ordered = [s for s in PIPELINE_ORDER if s in all_deps]

        # Add any scripts not in standard order at the end
        for script in scripts:
            if script not in ordered:
                ordered.append(script)

        return ordered

    def check_dependencies(self, script_name: str) -> Dict[str, bool]:
        """
        Check if script dependencies are satisfied.

        Args:
            script_name: Name of the script to check

        Returns:
            Dictionary of dependency name -> satisfied status
        """
        config = get_script_config(script_name)
        if not config:
            return {"error": f"Unknown script: {script_name}"}

        results = {}

        # Check script dependencies exist
        for dep in config.dependencies:
            dep_config = get_script_config(dep)
            if dep_config:
                results[f"script_{dep}"] = dep_config.path.exists()
            else:
                results[f"script_{dep}"] = False

        # Check required data files
        from config import RAW_DATA_DIR
        for data_file in config.required_data_files:
            file_path = RAW_DATA_DIR / data_file
            results[f"data_{data_file}"] = file_path.exists()

        return results

    # -------------------------------------------------------------------------
    # SCRIPT EXECUTION
    # -------------------------------------------------------------------------

    def run_script(
        self,
        script_name: str,
        timeout: int = 600,
        check_deps: bool = True,
        dry_run: bool = False
    ) -> ScriptExecution:
        """
        Run a single R script.

        Args:
            script_name: Name of the script (without path)
            timeout: Maximum execution time in seconds
            check_deps: Whether to verify dependencies first
            dry_run: If True, only check without executing

        Returns:
            ScriptExecution record
        """
        execution = ScriptExecution(
            script_name=script_name,
            status=PipelineStatus.PENDING,
            start_time=datetime.now()
        )

        config = get_script_config(script_name)
        if not config:
            execution.status = PipelineStatus.FAILED
            execution.skipped_reason = f"Unknown script: {script_name}"
            self.logger.error(f"Unknown script: {script_name}")
            return execution

        # Check dependencies
        if check_deps:
            deps = self.check_dependencies(script_name)
            unsatisfied = [k for k, v in deps.items() if not v]
            if unsatisfied:
                execution.status = PipelineStatus.SKIPPED
                execution.skipped_reason = f"Unsatisfied dependencies: {unsatisfied}"
                self.logger.warning(f"Skipping {script_name}: {execution.skipped_reason}")
                return execution

        # Dry run - just check
        if dry_run:
            self.logger.info(f"[DRY RUN] Would execute: {config.path}")
            execution.status = PipelineStatus.SKIPPED
            execution.skipped_reason = "Dry run"
            return execution

        # Execute script
        self.logger.info(f"Executing: {script_name}")
        execution.status = PipelineStatus.RUNNING

        result = run_r_script(
            script_path=config.path,
            working_dir=PROJECT_ROOT,
            timeout=timeout
        )

        execution.result = result
        execution.end_time = datetime.now()

        if result.success:
            execution.status = PipelineStatus.COMPLETED
            self.logger.info(
                f"Completed: {script_name} ({result.duration_seconds:.1f}s)"
            )
            if result.warnings:
                self.logger.warning(f"Warnings: {len(result.warnings)}")
        else:
            execution.status = PipelineStatus.FAILED
            self.logger.error(f"Failed: {script_name}")
            for error in result.errors[:3]:  # Log first 3 errors
                self.logger.error(f"  {error[:200]}")

        return execution

    def run_pipeline(
        self,
        scripts: Optional[List[str]] = None,
        stop_on_error: bool = True,
        dry_run: bool = False,
        skip_completed: bool = False
    ) -> PipelineRun:
        """
        Run the full analysis pipeline or specified scripts.

        Args:
            scripts: List of scripts to run (default: all in PIPELINE_ORDER)
            stop_on_error: Stop pipeline if a script fails
            dry_run: Only check without executing
            skip_completed: Skip scripts that have recent output

        Returns:
            PipelineRun record with all execution details
        """
        scripts = scripts or PIPELINE_ORDER
        ordered_scripts = self.get_dependency_order(scripts)

        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.current_run = PipelineRun(
            run_id=run_id,
            start_time=datetime.now(),
            status=PipelineStatus.RUNNING
        )

        self.logger.info("=" * 60)
        self.logger.info(f"Starting pipeline run: {run_id}")
        self.logger.info(f"Scripts to execute: {len(ordered_scripts)}")
        self.logger.info("=" * 60)

        progress = ProgressTracker(len(ordered_scripts), "Pipeline")

        for script_name in ordered_scripts:
            progress.update(script_name)

            execution = self.run_script(
                script_name,
                check_deps=True,
                dry_run=dry_run
            )
            self.current_run.scripts.append(execution)

            # Update counters
            if execution.status == PipelineStatus.COMPLETED:
                self.current_run.scripts_succeeded += 1
            elif execution.status == PipelineStatus.FAILED:
                self.current_run.scripts_failed += 1
                if stop_on_error:
                    self.logger.error(f"Pipeline stopped due to error in {script_name}")
                    break
            elif execution.status == PipelineStatus.SKIPPED:
                self.current_run.scripts_skipped += 1

        # Finalize run
        self.current_run.end_time = datetime.now()
        self.current_run.total_duration_seconds = (
            self.current_run.end_time - self.current_run.start_time
        ).total_seconds()

        if self.current_run.scripts_failed > 0:
            self.current_run.status = PipelineStatus.FAILED
        else:
            self.current_run.status = PipelineStatus.COMPLETED

        self.execution_history.append(self.current_run)

        # Log summary
        self.logger.info("=" * 60)
        self.logger.info("Pipeline Summary")
        self.logger.info(f"  Status: {self.current_run.status.value}")
        self.logger.info(f"  Succeeded: {self.current_run.scripts_succeeded}")
        self.logger.info(f"  Failed: {self.current_run.scripts_failed}")
        self.logger.info(f"  Skipped: {self.current_run.scripts_skipped}")
        self.logger.info(f"  Duration: {self.current_run.total_duration_seconds:.1f}s")
        self.logger.info("=" * 60)

        return self.current_run

    # -------------------------------------------------------------------------
    # OUTPUT VALIDATION
    # -------------------------------------------------------------------------

    def validate_outputs(self, script_name: str) -> Dict[str, Any]:
        """
        Validate that a script produced expected outputs.

        Args:
            script_name: Name of the script

        Returns:
            Validation results
        """
        config = get_script_config(script_name)
        if not config:
            return {"error": f"Unknown script: {script_name}"}

        results = {
            "script": script_name,
            "expected_outputs": config.outputs,
            "validated": [],
            "missing": [],
            "all_valid": True
        }

        # Check standard output locations based on script
        output_checks = self._get_output_checks(script_name)

        for output_name, output_path in output_checks.items():
            info = get_file_info(output_path)
            if info["exists"]:
                results["validated"].append({
                    "name": output_name,
                    "path": str(output_path),
                    "size": info["size_human"],
                    "modified": info["modified"]
                })
            else:
                results["missing"].append(output_name)
                results["all_valid"] = False

        return results

    def _get_output_checks(self, script_name: str) -> Dict[str, Path]:
        """Get expected output files for a script."""
        checks = {}

        if script_name == "6.coral-growth":
            checks["coral_growth_data"] = PROCESSED_DATA_DIR / "coral_growth.csv"
            checks["growth_figures"] = OUTPUT_DIR / "figures" / "coral"

        elif script_name == "7.coral-physiology":
            checks["physiology_figures"] = OUTPUT_DIR / "figures" / "coral" / "physio"

        elif script_name == "14.compile-manuscript-statistics":
            checks["stats_table"] = OUTPUT_DIR / "MANUSCRIPT_STATS_TABLE.csv"
            checks["stats_html"] = OUTPUT_DIR / "tables" / "MANUSCRIPT_STATISTICAL_TESTS.html"

        return checks

    # -------------------------------------------------------------------------
    # REPORTING
    # -------------------------------------------------------------------------

    def generate_report(self, run: Optional[PipelineRun] = None) -> Path:
        """
        Generate a detailed execution report.

        Args:
            run: Pipeline run to report on (default: most recent)

        Returns:
            Path to the generated report
        """
        run = run or self.current_run
        if not run:
            self.logger.error("No pipeline run to report on")
            return None

        report = {
            "run_id": run.run_id,
            "status": run.status.value,
            "start_time": run.start_time.isoformat(),
            "end_time": run.end_time.isoformat() if run.end_time else None,
            "duration_seconds": run.total_duration_seconds,
            "summary": {
                "total_scripts": len(run.scripts),
                "succeeded": run.scripts_succeeded,
                "failed": run.scripts_failed,
                "skipped": run.scripts_skipped
            },
            "scripts": []
        }

        for execution in run.scripts:
            script_report = {
                "name": execution.script_name,
                "status": execution.status.value,
                "start_time": execution.start_time.isoformat() if execution.start_time else None,
                "end_time": execution.end_time.isoformat() if execution.end_time else None,
            }

            if execution.result:
                script_report["duration_seconds"] = execution.result.duration_seconds
                script_report["return_code"] = execution.result.return_code
                script_report["warnings_count"] = len(execution.result.warnings)
                script_report["errors_count"] = len(execution.result.errors)
                if execution.result.errors:
                    script_report["errors"] = execution.result.errors[:5]

            if execution.skipped_reason:
                script_report["skipped_reason"] = execution.skipped_reason

            report["scripts"].append(script_report)

        # Save report
        report_path = OUTPUT_DIR / "logs" / f"pipeline_report_{run.run_id}.json"
        write_json_report(report, report_path)
        self.logger.info(f"Report saved to: {report_path}")

        return report_path

    # -------------------------------------------------------------------------
    # CONVENIENCE METHODS
    # -------------------------------------------------------------------------

    def run_growth_analysis(self) -> ScriptExecution:
        """Run the coral growth analysis script."""
        return self.run_script("6.coral-growth")

    def run_physiology_analysis(self) -> ScriptExecution:
        """Run the coral physiology analysis script."""
        # Ensure growth analysis is done first
        self.run_script("6.coral-growth")
        return self.run_script("7.coral-physiology")

    def run_community_analysis(self) -> List[ScriptExecution]:
        """Run all community-related analyses."""
        scripts = ["3.abundance", "4d.diversity", "5.fishes", "12.nmds_permanova_cafi"]
        results = []
        for script in scripts:
            results.append(self.run_script(script))
        return results

    def run_cafi_coral_analysis(self) -> PipelineRun:
        """Run the full CAFI-coral relationship analysis."""
        scripts = [
            "6.coral-growth",
            "7.coral-physiology",
            "8.coral-caffi"
        ]
        return self.run_pipeline(scripts)

    def compile_statistics(self) -> ScriptExecution:
        """Compile all manuscript statistics."""
        return self.run_script("14.compile-manuscript-statistics")

    def get_pipeline_status(self) -> Dict[str, Any]:
        """Get current pipeline status and recent history."""
        return {
            "current_run": {
                "run_id": self.current_run.run_id if self.current_run else None,
                "status": self.current_run.status.value if self.current_run else None,
            },
            "history_count": len(self.execution_history),
            "recent_runs": [
                {
                    "run_id": run.run_id,
                    "status": run.status.value,
                    "duration": run.total_duration_seconds
                }
                for run in self.execution_history[-5:]
            ]
        }


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    """Command-line interface for the pipeline agent."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CAFI MRB Analysis Pipeline Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pipeline_agent.py --run-all              Run full pipeline
  python pipeline_agent.py --script 6.coral-growth   Run specific script
  python pipeline_agent.py --check-deps 8.coral-caffi  Check dependencies
  python pipeline_agent.py --validate 6.coral-growth   Validate outputs
  python pipeline_agent.py --list                 List available scripts
        """
    )

    parser.add_argument(
        "--run-all", action="store_true",
        help="Run the complete analysis pipeline"
    )
    parser.add_argument(
        "--script", type=str,
        help="Run a specific script by name"
    )
    parser.add_argument(
        "--check-deps", type=str,
        help="Check dependencies for a script"
    )
    parser.add_argument(
        "--validate", type=str,
        help="Validate outputs for a script"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List all available scripts"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Check without executing"
    )
    parser.add_argument(
        "--no-stop-on-error", action="store_true",
        help="Continue pipeline even if a script fails"
    )

    args = parser.parse_args()

    agent = PipelineAgent()

    if args.list:
        print("\nAvailable Pipeline Scripts:")
        print("=" * 60)
        for name, config in PIPELINE_SCRIPTS.items():
            status = "OK" if config.path.exists() else "MISSING"
            print(f"  [{status}] {name}")
            print(f"         {config.description}")
        return

    if args.check_deps:
        deps = agent.check_dependencies(args.check_deps)
        print(f"\nDependencies for {args.check_deps}:")
        for name, satisfied in deps.items():
            status = "OK" if satisfied else "MISSING"
            print(f"  [{status}] {name}")
        return

    if args.validate:
        results = agent.validate_outputs(args.validate)
        print(f"\nOutput validation for {args.validate}:")
        print(f"  Validated: {len(results['validated'])}")
        print(f"  Missing: {len(results['missing'])}")
        for item in results["validated"]:
            print(f"    OK: {item['name']} ({item['size']})")
        for item in results["missing"]:
            print(f"    MISSING: {item}")
        return

    if args.script:
        execution = agent.run_script(
            args.script,
            dry_run=args.dry_run
        )
        print(f"\nScript: {execution.script_name}")
        print(f"Status: {execution.status.value}")
        if execution.result:
            print(f"Duration: {execution.result.duration_seconds:.1f}s")
        return

    if args.run_all:
        run = agent.run_pipeline(
            stop_on_error=not args.no_stop_on_error,
            dry_run=args.dry_run
        )
        agent.generate_report(run)
        print(f"\nPipeline completed: {run.status.value}")
        return

    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
