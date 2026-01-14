"""
Utility Functions for CAFI MRB Analysis Agents
===============================================
Shared utilities for logging, subprocess management, R integration, and common operations.

Author: CAFI Analysis Team
Created: 2025-01-10
"""

import subprocess
import logging
import json
import re
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
import csv

from config import (
    PROJECT_ROOT, SCRIPTS_DIR, OUTPUT_DIR, LOG_DIR,
    R_EXECUTABLE, LOG_LEVEL
)


# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logger(name: str, log_file: Optional[Path] = None) -> logging.Logger:
    """
    Set up a logger with console and optional file handlers.

    Args:
        name: Logger name
        log_file: Optional path to log file

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter(
        '%(asctime)s | %(name)s | %(levelname)s | %(message)s',
        datefmt='%H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            '%(asctime)s | %(name)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)

    return logger


# =============================================================================
# R SCRIPT EXECUTION
# =============================================================================

@dataclass
class RScriptResult:
    """Result of an R script execution."""
    success: bool
    script_name: str
    return_code: int
    stdout: str
    stderr: str
    duration_seconds: float
    timestamp: str
    warnings: List[str]
    errors: List[str]


def run_r_script(
    script_path: Path,
    working_dir: Optional[Path] = None,
    timeout: int = 600,
    capture_output: bool = True,
    env_vars: Optional[Dict[str, str]] = None
) -> RScriptResult:
    """
    Execute an R script and capture results.

    Args:
        script_path: Path to the R script
        working_dir: Working directory for execution (default: project root)
        timeout: Maximum execution time in seconds
        capture_output: Whether to capture stdout/stderr
        env_vars: Additional environment variables

    Returns:
        RScriptResult with execution details
    """
    import os
    import time

    if not script_path.exists():
        return RScriptResult(
            success=False,
            script_name=script_path.name,
            return_code=-1,
            stdout="",
            stderr=f"Script not found: {script_path}",
            duration_seconds=0,
            timestamp=datetime.now().isoformat(),
            warnings=[],
            errors=[f"Script not found: {script_path}"]
        )

    working_dir = working_dir or PROJECT_ROOT
    start_time = time.time()
    timestamp = datetime.now().isoformat()

    # Prepare environment
    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)

    try:
        result = subprocess.run(
            [R_EXECUTABLE, str(script_path)],
            cwd=str(working_dir),
            capture_output=capture_output,
            text=True,
            timeout=timeout,
            env=env
        )

        duration = time.time() - start_time
        stdout = result.stdout if capture_output else ""
        stderr = result.stderr if capture_output else ""

        # Parse warnings and errors from output
        warnings = extract_r_warnings(stdout + stderr)
        errors = extract_r_errors(stderr)

        return RScriptResult(
            success=result.returncode == 0,
            script_name=script_path.name,
            return_code=result.returncode,
            stdout=stdout,
            stderr=stderr,
            duration_seconds=duration,
            timestamp=timestamp,
            warnings=warnings,
            errors=errors
        )

    except subprocess.TimeoutExpired:
        duration = time.time() - start_time
        return RScriptResult(
            success=False,
            script_name=script_path.name,
            return_code=-1,
            stdout="",
            stderr=f"Script timed out after {timeout} seconds",
            duration_seconds=duration,
            timestamp=timestamp,
            warnings=[],
            errors=[f"Timeout after {timeout} seconds"]
        )

    except Exception as e:
        duration = time.time() - start_time
        return RScriptResult(
            success=False,
            script_name=script_path.name,
            return_code=-1,
            stdout="",
            stderr=str(e),
            duration_seconds=duration,
            timestamp=timestamp,
            warnings=[],
            errors=[str(e)]
        )


def extract_r_warnings(output: str) -> List[str]:
    """Extract warning messages from R output."""
    warnings = []
    pattern = r"Warning[^:]*:(.+?)(?=Warning|Error|$)"
    matches = re.findall(pattern, output, re.DOTALL | re.IGNORECASE)
    for match in matches:
        warning = match.strip()[:500]  # Limit length
        if warning:
            warnings.append(warning)
    return warnings


def extract_r_errors(output: str) -> List[str]:
    """Extract error messages from R output."""
    errors = []
    pattern = r"Error[^:]*:(.+?)(?=Error|Warning|$)"
    matches = re.findall(pattern, output, re.DOTALL | re.IGNORECASE)
    for match in matches:
        error = match.strip()[:500]
        if error:
            errors.append(error)
    return errors


# =============================================================================
# STATISTICAL RESULTS PARSING
# =============================================================================

@dataclass
class StatisticalResult:
    """Parsed statistical test result."""
    test_name: str
    test_type: str
    statistic_name: str
    statistic_value: float
    p_value: float
    df: Optional[str] = None
    effect_size: Optional[float] = None
    confidence_interval: Optional[Tuple[float, float]] = None
    interpretation: Optional[str] = None
    source_script: Optional[str] = None


def parse_lmm_anova(output: str) -> List[StatisticalResult]:
    """
    Parse ANOVA table from lmerTest output.

    Example output to parse:
    Analysis of Deviance Table (Type III Wald chisquare tests)
    Response: growth_vol_b
                     Chisq Df Pr(>Chisq)
    (Intercept)    147.93  1  < 2.2e-16 ***
    treatment        2.64  2     0.2674
    """
    results = []
    lines = output.split('\n')

    in_table = False
    for line in lines:
        if 'Chisq' in line and 'Df' in line:
            in_table = True
            continue

        if in_table and line.strip():
            # Parse table row
            parts = line.split()
            if len(parts) >= 4:
                try:
                    term = parts[0]
                    chi_sq = float(parts[1])
                    df = parts[2]
                    p_val_str = parts[3]

                    # Handle p-value parsing
                    if '<' in p_val_str:
                        p_value = 0.0001  # Approximate very small p-values
                    else:
                        p_value = float(p_val_str)

                    results.append(StatisticalResult(
                        test_name=term,
                        test_type="Type III Wald Chi-square",
                        statistic_name="Chi-square",
                        statistic_value=chi_sq,
                        p_value=p_value,
                        df=df
                    ))
                except (ValueError, IndexError):
                    continue

    return results


def parse_permanova_result(output: str) -> Optional[StatisticalResult]:
    """
    Parse PERMANOVA result from vegan::adonis2 output.

    Example output:
    Permutation test for adonis under reduced model
    Terms added sequentially (first to last)
    Permutation: free
    Number of permutations: 999

    adonis2(formula = comm ~ treatment, data = meta, method = "bray")
             Df SumOfSqs      R2      F Pr(>F)
    treatment  2   0.7834 0.16742 2.0102  0.015 *
    Residual  20   3.8940 0.83258
    Total     22   4.6774 1.00000
    """
    lines = output.split('\n')

    for line in lines:
        if 'treatment' in line.lower() and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    df = parts[1]
                    f_stat = float(parts[4])
                    p_value = float(parts[5])

                    # Extract R2 if available
                    r2 = float(parts[3]) if len(parts) > 3 else None

                    return StatisticalResult(
                        test_name="PERMANOVA",
                        test_type="PERMANOVA (adonis2)",
                        statistic_name="F",
                        statistic_value=f_stat,
                        p_value=p_value,
                        df=df,
                        effect_size=r2,
                        interpretation=f"R² = {r2:.3f}" if r2 else None
                    )
                except (ValueError, IndexError):
                    continue

    return None


# =============================================================================
# FILE OPERATIONS
# =============================================================================

def read_csv_safely(file_path: Path) -> Tuple[Optional[List[Dict]], Optional[str]]:
    """
    Read a CSV file with error handling.

    Returns:
        Tuple of (data as list of dicts, error message or None)
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            data = list(reader)
        return data, None
    except FileNotFoundError:
        return None, f"File not found: {file_path}"
    except Exception as e:
        return None, f"Error reading {file_path}: {str(e)}"


def write_json_report(data: Dict[str, Any], output_path: Path) -> bool:
    """Write a JSON report with proper formatting."""
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str)
        return True
    except Exception as e:
        print(f"Error writing JSON: {e}")
        return False


def get_file_info(file_path: Path) -> Dict[str, Any]:
    """Get information about a file."""
    if not file_path.exists():
        return {"exists": False, "path": str(file_path)}

    stat = file_path.stat()
    return {
        "exists": True,
        "path": str(file_path),
        "size_bytes": stat.st_size,
        "size_human": format_bytes(stat.st_size),
        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "created": datetime.fromtimestamp(stat.st_ctime).isoformat(),
    }


def format_bytes(size: int) -> str:
    """Format bytes to human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


# =============================================================================
# DATA VALIDATION
# =============================================================================

def validate_csv_schema(
    file_path: Path,
    required_columns: List[str]
) -> Tuple[bool, List[str]]:
    """
    Validate that a CSV file has required columns.

    Returns:
        Tuple of (is_valid, list of missing columns)
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader)

        # Clean header (remove BOM, whitespace)
        header = [col.strip().replace('\ufeff', '') for col in header]

        missing = [col for col in required_columns if col not in header]
        return len(missing) == 0, missing

    except Exception as e:
        return False, [f"Error reading file: {str(e)}"]


def check_data_quality(file_path: Path) -> Dict[str, Any]:
    """
    Check data quality metrics for a CSV file.

    Returns:
        Dictionary with quality metrics
    """
    data, error = read_csv_safely(file_path)
    if error:
        return {"error": error}

    if not data:
        return {"error": "Empty file", "row_count": 0}

    # Calculate metrics
    row_count = len(data)
    columns = list(data[0].keys())

    # Check for missing values
    missing_by_column = {}
    for col in columns:
        missing_count = sum(1 for row in data if not row.get(col) or row[col] == '')
        missing_by_column[col] = {
            "count": missing_count,
            "percentage": (missing_count / row_count) * 100
        }

    # Check for duplicate rows
    unique_rows = set(tuple(row.values()) for row in data)
    duplicate_count = row_count - len(unique_rows)

    return {
        "file": str(file_path),
        "row_count": row_count,
        "column_count": len(columns),
        "columns": columns,
        "missing_values": missing_by_column,
        "duplicate_rows": duplicate_count,
        "quality_score": calculate_quality_score(missing_by_column, duplicate_count, row_count)
    }


def calculate_quality_score(
    missing_by_column: Dict,
    duplicate_count: int,
    row_count: int
) -> float:
    """Calculate a 0-100 quality score."""
    if row_count == 0:
        return 0.0

    # Penalize for missing values
    avg_missing_pct = sum(v["percentage"] for v in missing_by_column.values()) / len(missing_by_column)

    # Penalize for duplicates
    duplicate_pct = (duplicate_count / row_count) * 100

    # Calculate score
    score = 100 - avg_missing_pct - (duplicate_pct * 0.5)
    return max(0, min(100, score))


# =============================================================================
# RESULT FORMATTING
# =============================================================================

def format_p_value(p: float, threshold: float = 0.001) -> str:
    """Format p-value for display."""
    if p < threshold:
        return f"< {threshold}"
    elif p < 0.01:
        return f"{p:.4f}"
    elif p < 0.05:
        return f"{p:.3f}"
    else:
        return f"{p:.2f}"


def format_statistic(value: float, name: str = "") -> str:
    """Format a statistical value for display."""
    if abs(value) >= 100:
        return f"{value:.1f}"
    elif abs(value) >= 10:
        return f"{value:.2f}"
    else:
        return f"{value:.3f}"


def generate_result_summary(results: List[StatisticalResult]) -> str:
    """Generate a markdown summary of statistical results."""
    lines = ["# Statistical Results Summary\n"]

    for result in results:
        lines.append(f"## {result.test_name}")
        lines.append(f"- **Test Type:** {result.test_type}")
        lines.append(f"- **{result.statistic_name}:** {format_statistic(result.statistic_value)}")
        lines.append(f"- **p-value:** {format_p_value(result.p_value)}")

        if result.df:
            lines.append(f"- **df:** {result.df}")
        if result.effect_size:
            lines.append(f"- **Effect Size:** {result.effect_size:.3f}")
        if result.interpretation:
            lines.append(f"- **Interpretation:** {result.interpretation}")

        lines.append("")

    return "\n".join(lines)


# =============================================================================
# PROGRESS TRACKING
# =============================================================================

class ProgressTracker:
    """Track progress of multi-step operations."""

    def __init__(self, total_steps: int, description: str = "Progress"):
        self.total_steps = total_steps
        self.current_step = 0
        self.description = description
        self.start_time = datetime.now()
        self.step_times: List[float] = []

    def update(self, step_name: str = "") -> None:
        """Update progress by one step."""
        import time
        self.current_step += 1
        self.step_times.append(time.time())

        pct = (self.current_step / self.total_steps) * 100
        print(f"[{self.current_step}/{self.total_steps}] {pct:.0f}% - {step_name}")

    def get_eta(self) -> Optional[str]:
        """Estimate time remaining."""
        if self.current_step == 0:
            return None

        elapsed = (datetime.now() - self.start_time).total_seconds()
        avg_time_per_step = elapsed / self.current_step
        remaining_steps = self.total_steps - self.current_step
        eta_seconds = avg_time_per_step * remaining_steps

        if eta_seconds < 60:
            return f"{eta_seconds:.0f} seconds"
        elif eta_seconds < 3600:
            return f"{eta_seconds/60:.1f} minutes"
        else:
            return f"{eta_seconds/3600:.1f} hours"

    def summary(self) -> Dict[str, Any]:
        """Get progress summary."""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        return {
            "description": self.description,
            "completed": self.current_step,
            "total": self.total_steps,
            "percentage": (self.current_step / self.total_steps) * 100,
            "elapsed_seconds": elapsed,
            "eta": self.get_eta()
        }


if __name__ == "__main__":
    # Test utilities
    print("Testing CAFI Agent Utilities")
    print("=" * 50)

    # Test logger
    logger = setup_logger("test", LOG_DIR / "test.log")
    logger.info("Logger initialized successfully")

    # Test file operations
    print("\nTesting file operations...")
    from config import REQUIRED_DATA_FILES
    for name, config in list(REQUIRED_DATA_FILES.items())[:2]:
        info = get_file_info(config.path)
        print(f"  {name}: {'exists' if info['exists'] else 'MISSING'}")

    print("\nUtilities test complete!")
