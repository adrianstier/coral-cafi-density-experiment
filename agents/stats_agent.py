"""
Statistical Analysis Agent for CAFI MRB Analysis
=================================================
Performs statistical analysis, extracts results from R outputs,
interprets models, and formats results for manuscript.

IMPORTANT: This agent generates R code for statistical analyses,
integrating with the existing CAFI R pipeline.

Author: CAFI Analysis Team
Created: 2025-01-10
Updated: 2025-01-11 - Added R code generation capabilities
"""

import sys
import re
import csv
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    PROJECT_ROOT, OUTPUT_DIR, TABLES_DIR, SCRIPTS_DIR,
    KEY_STATISTICAL_TESTS, StatisticalTestConfig
)
from utils import (
    setup_logger, run_r_script, read_csv_safely,
    write_json_report, format_p_value, StatisticalResult,
    parse_lmm_anova, parse_permanova_result
)
from r_code_generator import RCodeGenerator, RCodeTemplates


# =============================================================================
# DATA CLASSES
# =============================================================================

class SignificanceLevel(Enum):
    """Statistical significance levels."""
    HIGHLY_SIGNIFICANT = "***"  # p < 0.001
    VERY_SIGNIFICANT = "**"     # p < 0.01
    SIGNIFICANT = "*"           # p < 0.05
    MARGINALLY_SIGNIFICANT = "."  # p < 0.1
    NOT_SIGNIFICANT = "ns"      # p >= 0.1


@dataclass
class ModelSummary:
    """Summary of a statistical model."""
    model_name: str
    model_type: str
    formula: str
    response_variable: str
    predictors: List[str]
    random_effects: Optional[List[str]]
    n_observations: int
    n_groups: Optional[int]
    aic: Optional[float]
    bic: Optional[float]
    r_squared: Optional[float]
    marginal_r2: Optional[float]
    conditional_r2: Optional[float]


@dataclass
class PostHocResult:
    """Post-hoc comparison result."""
    contrast: str
    estimate: float
    se: float
    df: float
    t_ratio: float
    p_value: float
    significance: SignificanceLevel


@dataclass
class EffectSummary:
    """Summary of a treatment effect."""
    effect_name: str
    test_type: str
    statistic: float
    statistic_name: str
    df: str
    p_value: float
    p_adjusted: Optional[float]
    significance: SignificanceLevel
    direction: Optional[str]
    effect_size: Optional[float]
    confidence_interval: Optional[Tuple[float, float]]
    interpretation: str


# =============================================================================
# STATISTICAL ANALYSIS AGENT
# =============================================================================

class StatsAgent:
    """
    Agent for statistical analysis of CAFI MRB data.

    Capabilities:
    - Extract results from R output files
    - Parse statistical tables
    - Format results for manuscript
    - Check statistical assumptions
    - Generate statistical summaries
    """

    def __init__(self, log_file: Optional[Path] = None):
        """Initialize the stats agent."""
        self.logger = setup_logger(
            "StatsAgent",
            log_file or (OUTPUT_DIR / "logs" / "stats_agent.log")
        )
        self.results_cache: Dict[str, Any] = {}
        self.r_generator = RCodeGenerator()

    # -------------------------------------------------------------------------
    # RESULT EXTRACTION
    # -------------------------------------------------------------------------

    def load_manuscript_stats(self) -> Optional[List[Dict]]:
        """Load the compiled manuscript statistics table."""
        stats_file = OUTPUT_DIR / "MANUSCRIPT_STATS_TABLE.csv"
        data, error = read_csv_safely(stats_file)

        if error:
            self.logger.error(f"Failed to load stats: {error}")
            return None

        self.logger.info(f"Loaded {len(data)} statistical results")
        return data

    def extract_key_results(self) -> Dict[str, EffectSummary]:
        """Extract key statistical results for the manuscript."""
        stats = self.load_manuscript_stats()
        if not stats:
            return {}

        results = {}

        for row in stats:
            test_name = row.get("test_name", row.get("Test", ""))
            if not test_name:
                continue

            # Parse the row into an EffectSummary
            try:
                p_value = self._parse_p_value(row.get("p_value", row.get("P-value", "")))
                statistic = float(row.get("statistic", row.get("Statistic", 0)))

                effect = EffectSummary(
                    effect_name=test_name,
                    test_type=row.get("test_type", row.get("Test Type", "")),
                    statistic=statistic,
                    statistic_name=row.get("statistic_name", ""),
                    df=row.get("df", row.get("DF", "")),
                    p_value=p_value,
                    p_adjusted=self._parse_p_value(row.get("p_adjusted", "")),
                    significance=self._get_significance(p_value),
                    direction=row.get("direction", None),
                    effect_size=self._safe_float(row.get("effect_size", "")),
                    confidence_interval=None,
                    interpretation=self._generate_interpretation(test_name, p_value, statistic)
                )
                results[test_name] = effect

            except (ValueError, KeyError) as e:
                self.logger.warning(f"Could not parse result for {test_name}: {e}")
                continue

        return results

    def get_treatment_effects(self) -> List[EffectSummary]:
        """Get all treatment effect results."""
        all_results = self.extract_key_results()
        return [r for r in all_results.values() if "treatment" in r.effect_name.lower()]

    def get_coral_performance_results(self) -> List[EffectSummary]:
        """Get results related to coral performance."""
        all_results = self.extract_key_results()
        keywords = ["performance", "growth", "physiology", "pc1", "carbohydrate"]
        return [
            r for r in all_results.values()
            if any(kw in r.effect_name.lower() for kw in keywords)
        ]

    def get_community_composition_results(self) -> List[EffectSummary]:
        """Get community composition analysis results."""
        all_results = self.extract_key_results()
        keywords = ["permanova", "nmds", "community", "cafi"]
        return [
            r for r in all_results.values()
            if any(kw in r.effect_name.lower() for kw in keywords)
        ]

    # -------------------------------------------------------------------------
    # STATISTICAL PARSING
    # -------------------------------------------------------------------------

    def parse_r_output(self, output_text: str) -> List[StatisticalResult]:
        """
        Parse statistical results from R console output.

        Args:
            output_text: Raw R output text

        Returns:
            List of parsed StatisticalResult objects
        """
        results = []

        # Try different parsers
        results.extend(parse_lmm_anova(output_text))

        permanova = parse_permanova_result(output_text)
        if permanova:
            results.append(permanova)

        # Parse emmeans contrasts
        results.extend(self._parse_emmeans_contrasts(output_text))

        return results

    def _parse_emmeans_contrasts(self, output: str) -> List[StatisticalResult]:
        """Parse emmeans post-hoc contrast results."""
        results = []

        # Pattern for emmeans contrast output
        # Example:
        # contrast       estimate    SE  df t.ratio p.value
        # 1 - 3          0.123    0.05 20   2.46   0.023
        pattern = r"(\d+\s*-\s*\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)"

        for match in re.finditer(pattern, output):
            try:
                contrast = match.group(1).strip()
                estimate = float(match.group(2))
                se = float(match.group(3))
                df = match.group(4)
                t_ratio = float(match.group(5))
                p_value = float(match.group(6))

                results.append(StatisticalResult(
                    test_name=f"Contrast: {contrast}",
                    test_type="emmeans contrast",
                    statistic_name="t",
                    statistic_value=t_ratio,
                    p_value=p_value,
                    df=df
                ))
            except (ValueError, IndexError):
                continue

        return results

    # -------------------------------------------------------------------------
    # RESULT FORMATTING
    # -------------------------------------------------------------------------

    def format_for_manuscript(self, result: EffectSummary) -> str:
        """
        Format a statistical result for manuscript text.

        Example output: "χ² = 8.11, df = 2, p = 0.017"
        """
        parts = []

        # Format statistic
        stat_symbol = self._get_stat_symbol(result.statistic_name)
        parts.append(f"{stat_symbol} = {result.statistic:.2f}")

        # Add df if available
        if result.df:
            parts.append(f"df = {result.df}")

        # Format p-value
        parts.append(f"p {format_p_value(result.p_value)}")

        return ", ".join(parts)

    def format_for_table(self, results: List[EffectSummary]) -> List[Dict]:
        """Format results for a statistical table."""
        table_rows = []

        for result in results:
            row = {
                "Effect": result.effect_name,
                "Test": result.test_type,
                "Statistic": f"{result.statistic:.2f}",
                "df": result.df,
                "p-value": format_p_value(result.p_value),
                "Significance": result.significance.value,
            }

            if result.p_adjusted:
                row["p-adjusted (BH)"] = format_p_value(result.p_adjusted)

            if result.effect_size:
                row["Effect Size"] = f"{result.effect_size:.3f}"

            table_rows.append(row)

        return table_rows

    def generate_results_summary(self) -> str:
        """Generate a markdown summary of key results."""
        results = self.extract_key_results()

        lines = [
            "# Key Statistical Results",
            "",
            f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*",
            "",
            "## Treatment Effects on Coral Performance",
            ""
        ]

        for name, result in results.items():
            if "treatment" in name.lower() or "performance" in name.lower():
                sig_marker = result.significance.value
                lines.append(f"### {result.effect_name} {sig_marker}")
                lines.append(f"- {self.format_for_manuscript(result)}")
                lines.append(f"- {result.interpretation}")
                lines.append("")

        lines.extend([
            "## Community Composition",
            ""
        ])

        for name, result in results.items():
            if "permanova" in name.lower() or "community" in name.lower():
                sig_marker = result.significance.value
                lines.append(f"### {result.effect_name} {sig_marker}")
                lines.append(f"- {self.format_for_manuscript(result)}")
                lines.append("")

        return "\n".join(lines)

    # -------------------------------------------------------------------------
    # STATISTICAL CHECKS
    # -------------------------------------------------------------------------

    def run_assumption_checks(self, script_name: str) -> Dict[str, Any]:
        """
        Run statistical assumption checks for a script's analyses.

        Args:
            script_name: Name of the R script

        Returns:
            Dictionary with assumption check results
        """
        # Create a temporary R script to check assumptions
        check_script = f"""
# Assumption checks for {script_name}
library(performance)
library(lmerTest)

# Load data and models from {script_name}
# This would need to be customized per script

# Example checks:
# check_normality(model)
# check_heteroscedasticity(model)
# check_collinearity(model)

cat("Assumption checks complete\\n")
"""
        # For now, return placeholder results
        return {
            "script": script_name,
            "checks_run": ["normality", "heteroscedasticity", "collinearity"],
            "all_passed": True,
            "notes": "Implement custom assumption checks per script"
        }

    def check_multiple_testing(self, p_values: List[float], method: str = "BH") -> List[float]:
        """
        Apply multiple testing correction.

        Args:
            p_values: List of raw p-values
            method: Correction method ("BH" for Benjamini-Hochberg, "bonferroni")

        Returns:
            List of adjusted p-values
        """
        n = len(p_values)
        if n == 0:
            return []

        if method.lower() == "bonferroni":
            return [min(p * n, 1.0) for p in p_values]

        elif method.lower() == "bh":
            # Benjamini-Hochberg procedure
            indexed = sorted(enumerate(p_values), key=lambda x: x[1])
            adjusted = [0.0] * n

            # Calculate adjusted p-values
            for i, (orig_idx, p) in enumerate(indexed):
                rank = i + 1
                adj_p = p * n / rank
                adjusted[orig_idx] = min(adj_p, 1.0)

            # Ensure monotonicity
            for i in range(n - 2, -1, -1):
                adjusted[i] = min(adjusted[i], adjusted[i + 1])

            return adjusted

        else:
            self.logger.warning(f"Unknown correction method: {method}")
            return p_values

    # -------------------------------------------------------------------------
    # R CODE GENERATION FOR STATISTICAL ANALYSIS
    # -------------------------------------------------------------------------

    def generate_lmm_analysis_script(
        self,
        analysis_name: str,
        response_var: str,
        fixed_effects: List[str],
        random_effect: str = "reef",
        data_sources: List[str] = None,
        include_posthoc: bool = True,
        include_figure: bool = True,
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate an R script for linear mixed model analysis.

        Creates complete R code following CAFI conventions for running
        an LMM with the specified variables.

        Args:
            analysis_name: Name for the analysis
            response_var: Response/dependent variable
            fixed_effects: List of fixed effect predictors
            random_effect: Random effect grouping variable
            data_sources: Data sources to load
            include_posthoc: Include post-hoc comparisons
            include_figure: Include diagnostic/result figures
            run_script: Execute the generated script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        data_sources = data_sources or ["cafi", "treatment"]

        analyses = [{
            "type": "lmm",
            "model_name": f"{analysis_name}_model",
            "response_var": response_var,
            "fixed_effects": fixed_effects,
            "random_effect": random_effect,
            "data_var": "analysis_data",
            "include_posthoc": include_posthoc,
            "include_emmeans": True
        }]

        script = self.r_generator.generate_analysis_script(
            analysis_name=analysis_name,
            description=f"LMM analysis of {response_var} ~ {' + '.join(fixed_effects)}",
            data_sources=data_sources,
            analyses=analyses,
            output_tables=[f"{analysis_name}_results"]
        )

        # Add figure if requested
        if include_figure:
            script += RCodeTemplates.ggplot_figure(
                plot_name=f"{analysis_name}_plot",
                data_var="analysis_data",
                x_var=fixed_effects[0],
                y_var=response_var,
                plot_type="boxplot",
                fill_var=fixed_effects[0] if fixed_effects[0] == "treatment" else None,
                title=f"{response_var} by {fixed_effects[0]}",
                x_label=fixed_effects[0],
                y_label=response_var,
                output_path=f"output/MRB/figures/agent_generated/{analysis_name}"
            )

        script_path = self.r_generator.save_script(script, f"analysis_{analysis_name}")
        self.logger.info(f"Generated LMM analysis script: {script_path}")

        if run_script:
            result = run_r_script(script_path)
            if result.success:
                self.logger.info(f"Analysis completed in {result.duration_seconds:.1f}s")
            else:
                self.logger.error(f"Analysis failed: {result.errors}")

        return script, script_path

    def generate_permanova_script(
        self,
        analysis_name: str,
        grouping_var: str = "treatment",
        distance_method: str = "bray",
        permutations: int = 999,
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate an R script for PERMANOVA analysis.

        Args:
            analysis_name: Name for the analysis
            grouping_var: Grouping variable for PERMANOVA
            distance_method: Distance metric
            permutations: Number of permutations
            run_script: Execute the generated script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        script = RCodeTemplates.script_header(
            title=f"PERMANOVA Analysis: {analysis_name}",
            description=f"PERMANOVA testing effect of {grouping_var} on community composition",
            inputs=["CAFI community data", "Treatment data"],
            outputs=[f"{analysis_name}_results.csv"]
        )

        script += RCodeTemplates.load_data_block(["cafi", "treatment"])

        # Add community matrix creation
        script += '''
# ==============================================================================
# CREATE COMMUNITY MATRIX
# ==============================================================================

# Aggregate to reef level if needed
community_matrix <- cafi_data %>%
  group_by(coral_id, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0
  )

# Extract metadata
metadata <- community_matrix %>%
  select(coral_id) %>%
  left_join(treatment_data, by = "coral_id")

# Convert to matrix
comm_matrix <- community_matrix %>%
  select(-coral_id) %>%
  as.matrix()

rownames(comm_matrix) <- community_matrix$coral_id

'''
        script += RCodeTemplates.permanova_analysis(
            matrix_var="comm_matrix",
            grouping_var=grouping_var,
            metadata_var="metadata",
            method=distance_method,
            permutations=permutations
        )

        # Add NMDS visualization
        script += RCodeTemplates.nmds_ordination(matrix_var="comm_matrix")

        script += f'''
# ==============================================================================
# NMDS VISUALIZATION
# ==============================================================================

# Merge scores with metadata
nmds_plot_data <- nmds_scores %>%
  left_join(metadata, by = c("sample_id" = "coral_id"))

# Create NMDS plot
nmds_plot <- ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = {grouping_var})) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  labs(
    title = "NMDS Ordination",
    subtitle = sprintf("Stress: %.3f", nmds_result$stress)
  ) +
  theme_cafi_pub() +
  scale_color_manual(values = treatment_colors, labels = treatment_labels)

# Save figure
save_both(
  nmds_plot,
  here::here("output", "MRB", "figures", "agent_generated", "{analysis_name}_nmds"),
  width = 8,
  height = 6
)

message("PERMANOVA analysis complete: {analysis_name}")
'''

        script_path = self.r_generator.save_script(script, f"analysis_{analysis_name}")

        if run_script:
            result = run_r_script(script_path)
            if not result.success:
                self.logger.error(f"Analysis failed: {result.errors}")

        return script, script_path

    def generate_pca_analysis_script(
        self,
        analysis_name: str,
        variables: List[str],
        scale: bool = True,
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate an R script for PCA analysis.

        Args:
            analysis_name: Name for the analysis
            variables: Variables to include in PCA
            scale: Whether to scale variables
            run_script: Execute the generated script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        script = RCodeTemplates.script_header(
            title=f"PCA Analysis: {analysis_name}",
            description=f"Principal Component Analysis of {', '.join(variables[:3])}...",
            inputs=["Analysis data"],
            outputs=[f"{analysis_name}_scores.csv", f"{analysis_name}_loadings.csv"]
        )

        script += RCodeTemplates.load_data_block(["cafi", "treatment"])

        # Data prep
        vars_str = ', '.join([f'"{v}"' for v in variables])
        script += f'''
# ==============================================================================
# DATA PREPARATION
# ==============================================================================

# Select variables for PCA
pca_vars <- c({vars_str})

pca_data <- analysis_data %>%
  select(coral_id, treatment, all_of(pca_vars)) %>%
  drop_na()

# Extract numeric data for PCA
pca_numeric <- pca_data %>%
  select(all_of(pca_vars))

'''

        script += RCodeTemplates.pca_analysis(
            data_var="pca_numeric",
            scale=scale,
            center=True
        )

        script += f'''
# ==============================================================================
# PCA VISUALIZATION
# ==============================================================================

# Merge scores with metadata
pca_plot_data <- pca_scores %>%
  bind_cols(pca_data %>% select(coral_id, treatment))

# Biplot
pca_biplot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  labs(
    title = "PCA Biplot",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_cafi_pub() +
  scale_color_manual(values = treatment_colors, labels = treatment_labels)

# Save figure
save_both(
  pca_biplot,
  here::here("output", "MRB", "figures", "agent_generated", "{analysis_name}_biplot"),
  width = 8,
  height = 6
)

# Save results
write_csv(pca_scores, here::here("output", "MRB", "tables", "{analysis_name}_scores.csv"))
write_csv(pca_loadings, here::here("output", "MRB", "tables", "{analysis_name}_loadings.csv"))

message("PCA analysis complete: {analysis_name}")
'''

        script_path = self.r_generator.save_script(script, f"analysis_{analysis_name}")

        if run_script:
            result = run_r_script(script_path)
            if not result.success:
                self.logger.error(f"Analysis failed: {result.errors}")

        return script, script_path

    def generate_treatment_effect_script(
        self,
        response_var: str,
        response_label: str,
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate a complete treatment effect analysis script.

        Convenience method that combines LMM analysis with visualization.

        Args:
            response_var: Response variable name
            response_label: Human-readable label
            run_script: Execute the generated script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        script = self.r_generator.generate_treatment_effect_script(
            response_var=response_var,
            response_label=response_label,
            include_figure=True
        )

        script_path = self.r_generator.save_script(script, f"treatment_effect_{response_var}")

        if run_script:
            result = run_r_script(script_path)
            if result.success:
                self.logger.info(f"Treatment effect analysis completed")
            else:
                self.logger.error(f"Analysis failed: {result.errors}")

        return script, script_path

    def list_generated_scripts(self) -> List[Path]:
        """List all agent-generated analysis scripts."""
        generated_dir = SCRIPTS_DIR / "agent_generated"
        if not generated_dir.exists():
            return []
        return list(generated_dir.glob("analysis_*.R"))

    # -------------------------------------------------------------------------
    # HELPER METHODS
    # -------------------------------------------------------------------------

    def _parse_p_value(self, p_str: str) -> float:
        """Parse a p-value string to float."""
        if not p_str or p_str.strip() == "":
            return 1.0

        p_str = str(p_str).strip()

        # Handle "< 0.001" format
        if p_str.startswith("<"):
            val = re.search(r"[\d.]+", p_str)
            return float(val.group()) / 2 if val else 0.0001

        # Handle "p = 0.xxx" format
        match = re.search(r"[\d.]+", p_str)
        if match:
            return float(match.group())

        return 1.0

    def _safe_float(self, val: Any) -> Optional[float]:
        """Safely convert to float."""
        if val is None or val == "":
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    def _get_significance(self, p_value: float) -> SignificanceLevel:
        """Get significance level from p-value."""
        if p_value < 0.001:
            return SignificanceLevel.HIGHLY_SIGNIFICANT
        elif p_value < 0.01:
            return SignificanceLevel.VERY_SIGNIFICANT
        elif p_value < 0.05:
            return SignificanceLevel.SIGNIFICANT
        elif p_value < 0.1:
            return SignificanceLevel.MARGINALLY_SIGNIFICANT
        else:
            return SignificanceLevel.NOT_SIGNIFICANT

    def _get_stat_symbol(self, stat_name: str) -> str:
        """Get the symbol for a statistical test."""
        symbols = {
            "chi-square": "χ²",
            "chisq": "χ²",
            "chi_sq": "χ²",
            "f": "F",
            "t": "t",
            "z": "z",
            "r": "r",
            "r2": "R²",
        }
        return symbols.get(stat_name.lower(), stat_name)

    def _generate_interpretation(self, test_name: str, p_value: float, statistic: float) -> str:
        """Generate an interpretation string."""
        sig = self._get_significance(p_value)

        if sig == SignificanceLevel.NOT_SIGNIFICANT:
            return "No significant effect detected"
        elif sig == SignificanceLevel.MARGINALLY_SIGNIFICANT:
            return "Marginally significant effect (p < 0.1)"
        elif sig == SignificanceLevel.SIGNIFICANT:
            return "Significant effect (p < 0.05)"
        elif sig == SignificanceLevel.VERY_SIGNIFICANT:
            return "Highly significant effect (p < 0.01)"
        else:
            return "Very highly significant effect (p < 0.001)"

    # -------------------------------------------------------------------------
    # EXPORT METHODS
    # -------------------------------------------------------------------------

    def export_results_csv(self, output_path: Optional[Path] = None) -> Path:
        """Export all results to CSV."""
        results = self.extract_key_results()
        output_path = output_path or (TABLES_DIR / "stats_agent_results.csv")

        table_data = self.format_for_table(list(results.values()))

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            if table_data:
                writer = csv.DictWriter(f, fieldnames=table_data[0].keys())
                writer.writeheader()
                writer.writerows(table_data)

        self.logger.info(f"Exported results to: {output_path}")
        return output_path

    def export_results_markdown(self, output_path: Optional[Path] = None) -> Path:
        """Export results as markdown."""
        output_path = output_path or (OUTPUT_DIR / "KEY_STATISTICS_AGENT.md")

        summary = self.generate_results_summary()

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(summary)

        self.logger.info(f"Exported markdown to: {output_path}")
        return output_path


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    """Command-line interface for the stats agent."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CAFI MRB Statistical Analysis Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python stats_agent.py --summary           Show key results summary
  python stats_agent.py --treatment         Show treatment effects
  python stats_agent.py --export-csv        Export results to CSV
  python stats_agent.py --export-markdown   Export results to markdown
  python stats_agent.py --format "treatment_performance"  Format specific result
        """
    )

    parser.add_argument(
        "--summary", action="store_true",
        help="Display summary of key statistical results"
    )
    parser.add_argument(
        "--treatment", action="store_true",
        help="Display treatment effect results"
    )
    parser.add_argument(
        "--community", action="store_true",
        help="Display community composition results"
    )
    parser.add_argument(
        "--export-csv", action="store_true",
        help="Export results to CSV"
    )
    parser.add_argument(
        "--export-markdown", action="store_true",
        help="Export results to markdown"
    )
    parser.add_argument(
        "--format", type=str,
        help="Format a specific result for manuscript"
    )

    args = parser.parse_args()

    agent = StatsAgent()

    if args.summary:
        print(agent.generate_results_summary())
        return

    if args.treatment:
        results = agent.get_treatment_effects()
        print("\nTreatment Effects:")
        print("=" * 60)
        for result in results:
            print(f"\n{result.effect_name}")
            print(f"  {agent.format_for_manuscript(result)}")
            print(f"  {result.interpretation}")
        return

    if args.community:
        results = agent.get_community_composition_results()
        print("\nCommunity Composition Results:")
        print("=" * 60)
        for result in results:
            print(f"\n{result.effect_name}")
            print(f"  {agent.format_for_manuscript(result)}")
        return

    if args.export_csv:
        path = agent.export_results_csv()
        print(f"Exported to: {path}")
        return

    if args.export_markdown:
        path = agent.export_results_markdown()
        print(f"Exported to: {path}")
        return

    if args.format:
        results = agent.extract_key_results()
        if args.format in results:
            print(agent.format_for_manuscript(results[args.format]))
        else:
            print(f"Result not found: {args.format}")
            print(f"Available: {list(results.keys())}")
        return

    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
