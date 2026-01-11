"""
Data Validation Agent for CAFI MRB Analysis
============================================
Validates data integrity, checks for missing values, ensures reproducibility,
and generates data quality reports.

Author: CAFI Analysis Team
Created: 2025-01-10
"""

import sys
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import json
import csv

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    PROJECT_ROOT, DATA_DIR, OUTPUT_DIR, PROCESSED_DATA_DIR,
    RAW_DATA_DIR, REQUIRED_DATA_FILES, DataFileConfig,
    get_data_file_config
)
from utils import (
    setup_logger, read_csv_safely, write_json_report,
    get_file_info, check_data_quality, validate_csv_schema
)


# =============================================================================
# ENUMS AND DATA CLASSES
# =============================================================================

class ValidationStatus(Enum):
    """Status of a validation check."""
    PASSED = "passed"
    WARNING = "warning"
    FAILED = "failed"
    SKIPPED = "skipped"


class DataIssueType(Enum):
    """Types of data quality issues."""
    MISSING_FILE = "missing_file"
    MISSING_COLUMN = "missing_column"
    MISSING_VALUES = "missing_values"
    INVALID_TYPE = "invalid_type"
    OUT_OF_RANGE = "out_of_range"
    DUPLICATE_ROWS = "duplicate_rows"
    INCONSISTENT_FORMAT = "inconsistent_format"
    REFERENTIAL_INTEGRITY = "referential_integrity"


@dataclass
class ValidationCheck:
    """Result of a single validation check."""
    check_name: str
    status: ValidationStatus
    message: str
    details: Dict[str, Any] = field(default_factory=dict)
    severity: str = "info"  # info, warning, error


@dataclass
class DataIssue:
    """A data quality issue."""
    issue_type: DataIssueType
    file_name: str
    description: str
    affected_rows: Optional[int] = None
    affected_columns: Optional[List[str]] = None
    suggestion: Optional[str] = None


@dataclass
class ValidationReport:
    """Complete validation report."""
    timestamp: datetime
    total_checks: int
    passed: int
    warnings: int
    failed: int
    checks: List[ValidationCheck]
    issues: List[DataIssue]
    data_quality_score: float


# =============================================================================
# DATA VALIDATION AGENT
# =============================================================================

class ValidationAgent:
    """
    Agent for validating data quality and integrity.

    Capabilities:
    - Validate data file existence and structure
    - Check for missing values and duplicates
    - Validate data types and ranges
    - Check referential integrity between files
    - Generate comprehensive validation reports
    - Compute data quality scores
    """

    def __init__(self, log_file: Optional[Path] = None):
        """Initialize the validation agent."""
        self.logger = setup_logger(
            "ValidationAgent",
            log_file or (OUTPUT_DIR / "logs" / "validation_agent.log")
        )
        self.issues: List[DataIssue] = []
        self.checks: List[ValidationCheck] = []

    # -------------------------------------------------------------------------
    # FILE VALIDATION
    # -------------------------------------------------------------------------

    def validate_required_files(self) -> List[ValidationCheck]:
        """Validate that all required data files exist."""
        checks = []

        for name, config in REQUIRED_DATA_FILES.items():
            if config.path.exists():
                info = get_file_info(config.path)
                checks.append(ValidationCheck(
                    check_name=f"file_exists_{name}",
                    status=ValidationStatus.PASSED,
                    message=f"File exists: {config.path.name}",
                    details={
                        "path": str(config.path),
                        "size": info["size_human"],
                        "modified": info["modified"]
                    }
                ))
            else:
                checks.append(ValidationCheck(
                    check_name=f"file_exists_{name}",
                    status=ValidationStatus.FAILED,
                    message=f"Missing required file: {config.path.name}",
                    details={"expected_path": str(config.path)},
                    severity="error"
                ))
                self.issues.append(DataIssue(
                    issue_type=DataIssueType.MISSING_FILE,
                    file_name=config.path.name,
                    description=f"Required file not found: {config.path}",
                    suggestion=f"Ensure {config.path.name} is in the data directory"
                ))

        self.checks.extend(checks)
        return checks

    def validate_file_schema(self, file_name: str) -> List[ValidationCheck]:
        """
        Validate that a file has the expected columns.

        Args:
            file_name: Name of the file configuration

        Returns:
            List of validation checks
        """
        checks = []
        config = get_data_file_config(file_name)

        if not config:
            checks.append(ValidationCheck(
                check_name=f"schema_{file_name}",
                status=ValidationStatus.SKIPPED,
                message=f"Unknown file: {file_name}"
            ))
            return checks

        if not config.path.exists():
            checks.append(ValidationCheck(
                check_name=f"schema_{file_name}",
                status=ValidationStatus.SKIPPED,
                message=f"File not found: {config.path.name}"
            ))
            return checks

        if not config.columns:
            checks.append(ValidationCheck(
                check_name=f"schema_{file_name}",
                status=ValidationStatus.SKIPPED,
                message=f"No schema defined for {file_name}"
            ))
            return checks

        # Check schema
        is_valid, missing = validate_csv_schema(config.path, config.columns)

        if is_valid:
            checks.append(ValidationCheck(
                check_name=f"schema_{file_name}",
                status=ValidationStatus.PASSED,
                message=f"Schema valid: {config.path.name}",
                details={"expected_columns": config.columns}
            ))
        else:
            checks.append(ValidationCheck(
                check_name=f"schema_{file_name}",
                status=ValidationStatus.FAILED,
                message=f"Missing columns in {config.path.name}",
                details={"missing_columns": missing},
                severity="error"
            ))
            self.issues.append(DataIssue(
                issue_type=DataIssueType.MISSING_COLUMN,
                file_name=config.path.name,
                description=f"Missing columns: {', '.join(missing)}",
                affected_columns=missing,
                suggestion="Check if column names match expected format"
            ))

        self.checks.extend(checks)
        return checks

    # -------------------------------------------------------------------------
    # DATA QUALITY VALIDATION
    # -------------------------------------------------------------------------

    def validate_data_quality(self, file_name: str) -> List[ValidationCheck]:
        """
        Run data quality checks on a file.

        Args:
            file_name: Name of the file configuration

        Returns:
            List of validation checks
        """
        checks = []
        config = get_data_file_config(file_name)

        if not config or not config.path.exists():
            return checks

        quality = check_data_quality(config.path)

        if "error" in quality:
            checks.append(ValidationCheck(
                check_name=f"quality_{file_name}",
                status=ValidationStatus.FAILED,
                message=quality["error"]
            ))
            return checks

        # Check row count
        if quality["row_count"] > 0:
            checks.append(ValidationCheck(
                check_name=f"rows_{file_name}",
                status=ValidationStatus.PASSED,
                message=f"File has {quality['row_count']} rows",
                details={"row_count": quality["row_count"]}
            ))

        # Check for missing values
        for col, stats in quality.get("missing_values", {}).items():
            if stats["percentage"] > 20:
                checks.append(ValidationCheck(
                    check_name=f"missing_{file_name}_{col}",
                    status=ValidationStatus.WARNING,
                    message=f"High missing values in {col}: {stats['percentage']:.1f}%",
                    details=stats,
                    severity="warning"
                ))
                self.issues.append(DataIssue(
                    issue_type=DataIssueType.MISSING_VALUES,
                    file_name=config.path.name,
                    description=f"Column '{col}' has {stats['percentage']:.1f}% missing values",
                    affected_columns=[col],
                    affected_rows=stats["count"]
                ))
            elif stats["percentage"] > 0:
                checks.append(ValidationCheck(
                    check_name=f"missing_{file_name}_{col}",
                    status=ValidationStatus.PASSED,
                    message=f"Missing values in {col}: {stats['percentage']:.1f}%",
                    details=stats
                ))

        # Check for duplicates
        if quality.get("duplicate_rows", 0) > 0:
            checks.append(ValidationCheck(
                check_name=f"duplicates_{file_name}",
                status=ValidationStatus.WARNING,
                message=f"Found {quality['duplicate_rows']} duplicate rows",
                details={"duplicate_count": quality["duplicate_rows"]},
                severity="warning"
            ))
            self.issues.append(DataIssue(
                issue_type=DataIssueType.DUPLICATE_ROWS,
                file_name=config.path.name,
                description=f"File contains {quality['duplicate_rows']} duplicate rows",
                affected_rows=quality["duplicate_rows"]
            ))

        # Quality score
        checks.append(ValidationCheck(
            check_name=f"quality_score_{file_name}",
            status=ValidationStatus.PASSED if quality["quality_score"] >= 80 else ValidationStatus.WARNING,
            message=f"Data quality score: {quality['quality_score']:.1f}/100",
            details={"score": quality["quality_score"]}
        ))

        self.checks.extend(checks)
        return checks

    # -------------------------------------------------------------------------
    # SPECIFIC VALIDATIONS FOR CAFI DATA
    # -------------------------------------------------------------------------

    def validate_coral_ids(self) -> List[ValidationCheck]:
        """Validate coral ID consistency across files."""
        checks = []

        # Load coral IDs from different sources
        coral_ids_by_source = {}

        # From treatment file
        treatment_config = get_data_file_config("treatment")
        if treatment_config and treatment_config.path.exists():
            data, _ = read_csv_safely(treatment_config.path)
            if data:
                coral_ids_by_source["treatment"] = set(
                    row.get("coral_id", "") for row in data if row.get("coral_id")
                )

        # From physiology file
        physio_config = get_data_file_config("physiology")
        if physio_config and physio_config.path.exists():
            data, _ = read_csv_safely(physio_config.path)
            if data:
                coral_ids_by_source["physiology"] = set(
                    row.get("coral_id", "") for row in data if row.get("coral_id")
                )

        # Check consistency
        if len(coral_ids_by_source) >= 2:
            sources = list(coral_ids_by_source.keys())
            for i, source1 in enumerate(sources):
                for source2 in sources[i+1:]:
                    ids1 = coral_ids_by_source[source1]
                    ids2 = coral_ids_by_source[source2]

                    only_in_1 = ids1 - ids2
                    only_in_2 = ids2 - ids1

                    if only_in_1 or only_in_2:
                        checks.append(ValidationCheck(
                            check_name=f"coral_id_consistency_{source1}_{source2}",
                            status=ValidationStatus.WARNING,
                            message=f"Coral ID mismatch between {source1} and {source2}",
                            details={
                                f"only_in_{source1}": list(only_in_1)[:10],
                                f"only_in_{source2}": list(only_in_2)[:10],
                            },
                            severity="warning"
                        ))
                        self.issues.append(DataIssue(
                            issue_type=DataIssueType.REFERENTIAL_INTEGRITY,
                            file_name=f"{source1} vs {source2}",
                            description=f"Coral IDs don't match between files",
                            suggestion="Check for typos or missing records"
                        ))
                    else:
                        checks.append(ValidationCheck(
                            check_name=f"coral_id_consistency_{source1}_{source2}",
                            status=ValidationStatus.PASSED,
                            message=f"Coral IDs consistent between {source1} and {source2}"
                        ))

        self.checks.extend(checks)
        return checks

    def validate_treatment_values(self) -> List[ValidationCheck]:
        """Validate treatment values are valid (1, 3, or 6)."""
        checks = []
        valid_treatments = {"1", "3", "6"}

        treatment_config = get_data_file_config("treatment")
        if not treatment_config or not treatment_config.path.exists():
            return checks

        data, _ = read_csv_safely(treatment_config.path)
        if not data:
            return checks

        treatments = set(str(row.get("treatment", "")).strip() for row in data)
        invalid = treatments - valid_treatments - {""}

        if invalid:
            checks.append(ValidationCheck(
                check_name="treatment_values",
                status=ValidationStatus.FAILED,
                message=f"Invalid treatment values found",
                details={"invalid_values": list(invalid)},
                severity="error"
            ))
            self.issues.append(DataIssue(
                issue_type=DataIssueType.OUT_OF_RANGE,
                file_name="treatment",
                description=f"Invalid treatment values: {invalid}",
                suggestion="Treatment should be 1, 3, or 6"
            ))
        else:
            checks.append(ValidationCheck(
                check_name="treatment_values",
                status=ValidationStatus.PASSED,
                message="All treatment values are valid (1, 3, or 6)",
                details={"treatments": list(treatments)}
            ))

        self.checks.extend(checks)
        return checks

    def validate_numeric_ranges(self) -> List[ValidationCheck]:
        """Validate that numeric values are within expected ranges."""
        checks = []

        # Define expected ranges
        ranges = {
            "physiology": {
                "carbohydrate": (0, 100),
                "protein": (0, 100),
                "zooxanthellae": (0, 1e9),
                "afdw": (0, 100),
            }
        }

        physio_config = get_data_file_config("physiology")
        if physio_config and physio_config.path.exists():
            data, _ = read_csv_safely(physio_config.path)
            if data:
                for col, (min_val, max_val) in ranges.get("physiology", {}).items():
                    # Get values, handling various column name formats
                    values = []
                    for row in data:
                        for key in row.keys():
                            if col.lower() in key.lower():
                                try:
                                    val = float(row[key])
                                    values.append(val)
                                except (ValueError, TypeError):
                                    pass

                    if values:
                        out_of_range = [v for v in values if v < min_val or v > max_val]
                        if out_of_range:
                            checks.append(ValidationCheck(
                                check_name=f"range_{col}",
                                status=ValidationStatus.WARNING,
                                message=f"{len(out_of_range)} values out of range for {col}",
                                details={
                                    "expected_range": [min_val, max_val],
                                    "out_of_range_count": len(out_of_range)
                                }
                            ))
                        else:
                            checks.append(ValidationCheck(
                                check_name=f"range_{col}",
                                status=ValidationStatus.PASSED,
                                message=f"All {col} values within expected range"
                            ))

        self.checks.extend(checks)
        return checks

    # -------------------------------------------------------------------------
    # REPRODUCIBILITY VALIDATION
    # -------------------------------------------------------------------------

    def compute_file_checksums(self) -> Dict[str, str]:
        """Compute MD5 checksums for all data files."""
        checksums = {}

        for name, config in REQUIRED_DATA_FILES.items():
            if config.path.exists():
                with open(config.path, 'rb') as f:
                    checksum = hashlib.md5(f.read()).hexdigest()
                checksums[name] = checksum
            else:
                checksums[name] = "FILE_NOT_FOUND"

        return checksums

    def validate_reproducibility(self, expected_checksums: Optional[Dict[str, str]] = None) -> List[ValidationCheck]:
        """
        Validate data files against expected checksums.

        Args:
            expected_checksums: Dictionary of expected checksums

        Returns:
            List of validation checks
        """
        checks = []
        current_checksums = self.compute_file_checksums()

        if not expected_checksums:
            # Just report current checksums
            checks.append(ValidationCheck(
                check_name="checksums",
                status=ValidationStatus.PASSED,
                message="Computed checksums for reproducibility",
                details={"checksums": current_checksums}
            ))
        else:
            # Compare to expected
            for name, expected in expected_checksums.items():
                current = current_checksums.get(name)
                if current == expected:
                    checks.append(ValidationCheck(
                        check_name=f"checksum_{name}",
                        status=ValidationStatus.PASSED,
                        message=f"Checksum matches for {name}"
                    ))
                else:
                    checks.append(ValidationCheck(
                        check_name=f"checksum_{name}",
                        status=ValidationStatus.WARNING,
                        message=f"Checksum mismatch for {name}",
                        details={
                            "expected": expected,
                            "current": current
                        }
                    ))

        self.checks.extend(checks)
        return checks

    # -------------------------------------------------------------------------
    # COMPREHENSIVE VALIDATION
    # -------------------------------------------------------------------------

    def run_all_validations(self) -> ValidationReport:
        """
        Run all validation checks and generate a report.

        Returns:
            ValidationReport with all results
        """
        self.issues = []
        self.checks = []

        self.logger.info("Starting comprehensive data validation")

        # Run all validation checks
        self.validate_required_files()

        for file_name in REQUIRED_DATA_FILES.keys():
            self.validate_file_schema(file_name)
            self.validate_data_quality(file_name)

        self.validate_coral_ids()
        self.validate_treatment_values()
        self.validate_numeric_ranges()
        self.validate_reproducibility()

        # Compile report
        passed = sum(1 for c in self.checks if c.status == ValidationStatus.PASSED)
        warnings = sum(1 for c in self.checks if c.status == ValidationStatus.WARNING)
        failed = sum(1 for c in self.checks if c.status == ValidationStatus.FAILED)

        # Calculate quality score
        total = len(self.checks)
        quality_score = ((passed * 1.0) + (warnings * 0.5)) / total * 100 if total > 0 else 0

        report = ValidationReport(
            timestamp=datetime.now(),
            total_checks=total,
            passed=passed,
            warnings=warnings,
            failed=failed,
            checks=self.checks,
            issues=self.issues,
            data_quality_score=quality_score
        )

        self.logger.info(f"Validation complete: {passed} passed, {warnings} warnings, {failed} failed")
        self.logger.info(f"Data quality score: {quality_score:.1f}/100")

        return report

    # -------------------------------------------------------------------------
    # REPORTING
    # -------------------------------------------------------------------------

    def generate_report(self, report: ValidationReport, output_path: Optional[Path] = None) -> Path:
        """
        Generate a validation report file.

        Args:
            report: ValidationReport to save
            output_path: Output file path

        Returns:
            Path to the saved report
        """
        output_path = output_path or (OUTPUT_DIR / "validation_report.json")

        report_dict = {
            "timestamp": report.timestamp.isoformat(),
            "summary": {
                "total_checks": report.total_checks,
                "passed": report.passed,
                "warnings": report.warnings,
                "failed": report.failed,
                "quality_score": report.data_quality_score
            },
            "checks": [
                {
                    "name": c.check_name,
                    "status": c.status.value,
                    "message": c.message,
                    "severity": c.severity,
                    "details": c.details
                }
                for c in report.checks
            ],
            "issues": [
                {
                    "type": i.issue_type.value,
                    "file": i.file_name,
                    "description": i.description,
                    "affected_rows": i.affected_rows,
                    "affected_columns": i.affected_columns,
                    "suggestion": i.suggestion
                }
                for i in report.issues
            ]
        }

        write_json_report(report_dict, output_path)
        self.logger.info(f"Report saved to: {output_path}")

        return output_path

    def generate_markdown_report(self, report: ValidationReport, output_path: Optional[Path] = None) -> Path:
        """Generate a markdown validation report."""
        output_path = output_path or (OUTPUT_DIR / "VALIDATION_REPORT.md")

        lines = [
            "# Data Validation Report",
            "",
            f"*Generated: {report.timestamp.strftime('%Y-%m-%d %H:%M')}*",
            "",
            "## Summary",
            "",
            f"| Metric | Value |",
            f"|--------|-------|",
            f"| Total Checks | {report.total_checks} |",
            f"| Passed | {report.passed} |",
            f"| Warnings | {report.warnings} |",
            f"| Failed | {report.failed} |",
            f"| Quality Score | {report.data_quality_score:.1f}/100 |",
            "",
        ]

        # Issues section
        if report.issues:
            lines.extend([
                "## Issues Found",
                ""
            ])
            for issue in report.issues:
                lines.append(f"### {issue.issue_type.value}")
                lines.append(f"- **File:** {issue.file_name}")
                lines.append(f"- **Description:** {issue.description}")
                if issue.suggestion:
                    lines.append(f"- **Suggestion:** {issue.suggestion}")
                lines.append("")

        # Detailed checks
        lines.extend([
            "## Detailed Check Results",
            ""
        ])

        for check in report.checks:
            status_icon = {"passed": "OK", "warning": "WARN", "failed": "FAIL", "skipped": "SKIP"}
            icon = status_icon.get(check.status.value, "?")
            lines.append(f"- **[{icon}]** {check.check_name}: {check.message}")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write("\n".join(lines))

        return output_path


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    """Command-line interface for the validation agent."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CAFI MRB Data Validation Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validation_agent.py --validate-all      Run all validations
  python validation_agent.py --check-files       Check required files
  python validation_agent.py --check-quality     Check data quality
  python validation_agent.py --checksums         Compute file checksums
  python validation_agent.py --report            Generate full report
        """
    )

    parser.add_argument(
        "--validate-all", action="store_true",
        help="Run all validation checks"
    )
    parser.add_argument(
        "--check-files", action="store_true",
        help="Check for required files"
    )
    parser.add_argument(
        "--check-quality", action="store_true",
        help="Run data quality checks"
    )
    parser.add_argument(
        "--checksums", action="store_true",
        help="Compute and display file checksums"
    )
    parser.add_argument(
        "--report", action="store_true",
        help="Generate validation report"
    )
    parser.add_argument(
        "--output", type=str,
        help="Output path for report"
    )

    args = parser.parse_args()

    agent = ValidationAgent()

    if args.check_files:
        checks = agent.validate_required_files()
        print("\nRequired Files Check:")
        print("=" * 50)
        for check in checks:
            status = "OK" if check.status == ValidationStatus.PASSED else "MISSING"
            print(f"  [{status}] {check.message}")
        return

    if args.check_quality:
        print("\nData Quality Checks:")
        print("=" * 50)
        for file_name in REQUIRED_DATA_FILES.keys():
            checks = agent.validate_data_quality(file_name)
            for check in checks:
                if "quality_score" in check.check_name:
                    print(f"\n  {file_name}: {check.message}")
        return

    if args.checksums:
        checksums = agent.compute_file_checksums()
        print("\nFile Checksums (MD5):")
        print("=" * 50)
        for name, checksum in checksums.items():
            print(f"  {name}: {checksum}")
        return

    if args.validate_all or args.report:
        report = agent.run_all_validations()

        print("\n" + "=" * 60)
        print("VALIDATION REPORT")
        print("=" * 60)
        print(f"Total Checks: {report.total_checks}")
        print(f"  Passed:   {report.passed}")
        print(f"  Warnings: {report.warnings}")
        print(f"  Failed:   {report.failed}")
        print(f"\nData Quality Score: {report.data_quality_score:.1f}/100")

        if report.issues:
            print(f"\nIssues Found: {len(report.issues)}")
            for issue in report.issues[:5]:
                print(f"  - {issue.issue_type.value}: {issue.description}")

        if args.report:
            output_path = Path(args.output) if args.output else None
            json_path = agent.generate_report(report, output_path)
            md_path = agent.generate_markdown_report(report)
            print(f"\nReports saved:")
            print(f"  JSON: {json_path}")
            print(f"  Markdown: {md_path}")

        return

    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
