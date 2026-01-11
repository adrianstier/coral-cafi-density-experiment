"""
Figure Generation Agent for CAFI MRB Analysis
==============================================
Creates, updates, and manages publication-quality figures
for the coral reef analysis manuscript.

IMPORTANT: This agent generates R code for figures, which integrates
with the existing CAFI R analysis pipeline.

Author: CAFI Analysis Team
Created: 2025-01-10
Updated: 2025-01-11 - Added R code generation capabilities
"""

import sys
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import json

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    PROJECT_ROOT, OUTPUT_DIR, FIGURES_DIR, PUBLICATION_FIGURES_DIR,
    SCRIPTS_DIR, PUBLICATION_FIGURES, FigureConfig, get_figure_config
)
from utils import (
    setup_logger, run_r_script, get_file_info,
    write_json_report, ProgressTracker
)
from r_code_generator import RCodeGenerator, RCodeTemplates


# =============================================================================
# ENUMS AND DATA CLASSES
# =============================================================================

class FigureStatus(Enum):
    """Status of a figure."""
    NOT_FOUND = "not_found"
    OUTDATED = "outdated"
    CURRENT = "current"
    GENERATING = "generating"
    ERROR = "error"


class FigureFormat(Enum):
    """Supported figure formats."""
    PNG = "png"
    PDF = "pdf"
    SVG = "svg"
    TIFF = "tiff"


@dataclass
class FigureInfo:
    """Information about a figure file."""
    path: Path
    exists: bool
    format: FigureFormat
    size_bytes: int = 0
    modified: Optional[datetime] = None
    width_px: Optional[int] = None
    height_px: Optional[int] = None
    dpi: Optional[int] = None


@dataclass
class FigureGenerationResult:
    """Result of figure generation."""
    figure_name: str
    success: bool
    output_paths: List[Path]
    duration_seconds: float
    errors: List[str]
    warnings: List[str]


# =============================================================================
# FIGURE GENERATION AGENT
# =============================================================================

class FigureAgent:
    """
    Agent for generating and managing publication figures.

    Capabilities:
    - Generate individual or batch figures
    - Update figure styling
    - Check figure quality
    - Export in multiple formats
    - Track figure versions
    """

    def __init__(self, log_file: Optional[Path] = None):
        """Initialize the figure agent."""
        self.logger = setup_logger(
            "FigureAgent",
            log_file or (OUTPUT_DIR / "logs" / "figure_agent.log")
        )

        # Initialize R code generator
        self.r_generator = RCodeGenerator()

        # Standard figure dimensions (in inches)
        self.default_dimensions = {
            "single_column": (3.5, 3.5),
            "double_column": (7.0, 5.0),
            "full_page": (7.0, 9.0),
            "square": (6.0, 6.0),
            "wide": (10.0, 6.0),
        }

        # Standard DPI settings
        self.dpi_settings = {
            "screen": 72,
            "presentation": 150,
            "publication": 300,
            "high_quality": 600,
        }

    # -------------------------------------------------------------------------
    # FIGURE INVENTORY
    # -------------------------------------------------------------------------

    def scan_figures(self) -> Dict[str, List[FigureInfo]]:
        """
        Scan all figure directories and inventory existing figures.

        Returns:
            Dictionary mapping directory names to lists of FigureInfo
        """
        inventory = {}

        # Scan publication figures
        pub_figures = []
        if PUBLICATION_FIGURES_DIR.exists():
            for fig_path in PUBLICATION_FIGURES_DIR.glob("*"):
                if fig_path.suffix.lower() in [".png", ".pdf", ".svg", ".tiff"]:
                    info = self._get_figure_info(fig_path)
                    pub_figures.append(info)

        inventory["publication"] = pub_figures

        # Scan other figure directories
        figure_subdirs = ["coral", "community", "cafi", "fish", "nmds"]
        for subdir in figure_subdirs:
            subdir_path = FIGURES_DIR / subdir
            if subdir_path.exists():
                figures = []
                for fig_path in subdir_path.rglob("*.png"):
                    figures.append(self._get_figure_info(fig_path))
                inventory[subdir] = figures

        return inventory

    def get_figure_status(self, figure_number: int) -> Dict[str, Any]:
        """
        Get status of a publication figure.

        Args:
            figure_number: Figure number (1-8)

        Returns:
            Status dictionary
        """
        config = get_figure_config(figure_number)
        if not config:
            return {"error": f"Unknown figure: {figure_number}"}

        info = self._get_figure_info(config.output_path)

        return {
            "figure_number": figure_number,
            "name": config.name,
            "description": config.description,
            "source_script": config.script_source,
            "exists": info.exists,
            "status": FigureStatus.CURRENT.value if info.exists else FigureStatus.NOT_FOUND.value,
            "path": str(config.output_path),
            "size": f"{info.size_bytes / 1024:.1f} KB" if info.exists else None,
            "modified": info.modified.isoformat() if info.modified else None,
        }

    def get_all_figure_status(self) -> List[Dict[str, Any]]:
        """Get status of all publication figures."""
        return [self.get_figure_status(i) for i in range(1, 9)]

    # -------------------------------------------------------------------------
    # FIGURE GENERATION
    # -------------------------------------------------------------------------

    def generate_figure(
        self,
        figure_number: int,
        formats: List[FigureFormat] = None,
        dpi: int = 300,
        timeout: int = 180
    ) -> FigureGenerationResult:
        """
        Generate a specific publication figure.

        Args:
            figure_number: Figure number to generate
            formats: Output formats (default: PNG and PDF)
            dpi: Resolution for raster formats
            timeout: Maximum execution time

        Returns:
            FigureGenerationResult
        """
        formats = formats or [FigureFormat.PNG, FigureFormat.PDF]
        config = get_figure_config(figure_number)

        if not config:
            return FigureGenerationResult(
                figure_name=f"Figure_{figure_number}",
                success=False,
                output_paths=[],
                duration_seconds=0,
                errors=[f"Unknown figure number: {figure_number}"],
                warnings=[]
            )

        self.logger.info(f"Generating Figure {figure_number}: {config.name}")

        # Find the source script
        script_path = SCRIPTS_DIR / config.script_source
        if not script_path.exists():
            return FigureGenerationResult(
                figure_name=config.name,
                success=False,
                output_paths=[],
                duration_seconds=0,
                errors=[f"Source script not found: {script_path}"],
                warnings=[]
            )

        # Run the R script
        result = run_r_script(script_path, timeout=timeout)

        output_paths = []
        if result.success:
            # Check for generated files
            for fmt in formats:
                expected_path = config.output_path.with_suffix(f".{fmt.value}")
                if expected_path.exists():
                    output_paths.append(expected_path)

        return FigureGenerationResult(
            figure_name=config.name,
            success=result.success and len(output_paths) > 0,
            output_paths=output_paths,
            duration_seconds=result.duration_seconds,
            errors=result.errors,
            warnings=result.warnings
        )

    def generate_all_figures(
        self,
        figure_numbers: List[int] = None,
        parallel: bool = False
    ) -> List[FigureGenerationResult]:
        """
        Generate multiple figures.

        Args:
            figure_numbers: Figures to generate (default: all)
            parallel: Whether to run in parallel (not implemented)

        Returns:
            List of generation results
        """
        figure_numbers = figure_numbers or list(range(1, 9))
        results = []

        progress = ProgressTracker(len(figure_numbers), "Figure Generation")

        for fig_num in figure_numbers:
            progress.update(f"Figure {fig_num}")
            result = self.generate_figure(fig_num)
            results.append(result)

            if result.success:
                self.logger.info(f"Figure {fig_num} generated successfully")
            else:
                self.logger.error(f"Figure {fig_num} failed: {result.errors}")

        return results

    # -------------------------------------------------------------------------
    # FIGURE STYLING
    # -------------------------------------------------------------------------

    def create_style_script(self, output_path: Optional[Path] = None) -> Path:
        """
        Create an R script with standard figure styling.

        Returns:
            Path to the created style script
        """
        output_path = output_path or (SCRIPTS_DIR / "agent_figure_style.R")

        style_content = '''
# ============================================================================
# CAFI Figure Styling Standards (Agent Generated)
# ============================================================================
# Generated by FigureAgent on ''' + datetime.now().strftime("%Y-%m-%d") + '''
#
# Usage: source("agent_figure_style.R") at the start of figure scripts
# ============================================================================

library(ggplot2)
library(scales)

# Color palette for treatments
treatment_colors <- c(
  "1" = "#4DAF4A",   # Green for single coral
  "3" = "#377EB8",   # Blue for three corals
  "6" = "#E41A1C"    # Red for six corals
)

treatment_labels <- c(
  "1" = "1 Coral",
  "3" = "3 Corals",
  "6" = "6 Corals"
)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),

      # Text
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 2),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),

      # Legend
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = NA),

      # Strips (for facets)
      strip.background = element_rect(fill = "grey92", color = "black"),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

# Standard scale functions
scale_color_treatment <- function(...) {
  scale_color_manual(values = treatment_colors, labels = treatment_labels, ...)
}

scale_fill_treatment <- function(...) {
  scale_fill_manual(values = treatment_colors, labels = treatment_labels, ...)
}

# Save function with standard settings
save_pub_figure <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  # Ensure directory exists
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  # Save PNG
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  # Save PDF
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  message("Saved: ", filename, " (PNG + PDF)")
}

message("CAFI figure styling loaded successfully")
'''

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(style_content)

        self.logger.info(f"Created style script: {output_path}")
        return output_path

    # -------------------------------------------------------------------------
    # R CODE GENERATION FOR FIGURES
    # -------------------------------------------------------------------------

    def generate_figure_r_script(
        self,
        figure_name: str,
        description: str,
        data_sources: List[str],
        x_var: str,
        y_var: str,
        plot_type: str = "boxplot",
        color_var: Optional[str] = None,
        fill_var: Optional[str] = None,
        facet_var: Optional[str] = None,
        title: str = "",
        x_label: str = "",
        y_label: str = "",
        width: float = 8,
        height: float = 6,
        run_script: bool = False
    ) -> Tuple[str, Optional[Path]]:
        """
        Generate an R script for creating a figure.

        This method generates R code that follows CAFI pipeline conventions
        and can be run as part of the existing R workflow.

        Args:
            figure_name: Name for the figure (used in filenames)
            description: Description of what the figure shows
            data_sources: List of data sources (e.g., ["cafi", "treatment"])
            x_var: Variable for x-axis
            y_var: Variable for y-axis
            plot_type: Type of plot ("boxplot", "point", "line", "violin", "bar")
            color_var: Variable for color aesthetic
            fill_var: Variable for fill aesthetic
            facet_var: Variable for faceting
            title: Plot title
            x_label: X-axis label
            y_label: Y-axis label
            width: Figure width in inches
            height: Figure height in inches
            run_script: Whether to execute the generated script

        Returns:
            Tuple of (R script content, path to saved script or None)
        """
        plot_config = {
            "plot_name": f"{figure_name}_plot",
            "data_var": "analysis_data",
            "x_var": x_var,
            "y_var": y_var,
            "plot_type": plot_type,
            "title": title or description,
            "x_label": x_label or x_var,
            "y_label": y_label or y_var,
            "width": width,
            "height": height,
            "output_path": f"output/MRB/figures/agent_generated/{figure_name}"
        }

        if color_var:
            plot_config["color_var"] = color_var
        if fill_var:
            plot_config["fill_var"] = fill_var
        if facet_var:
            plot_config["facet_var"] = facet_var

        # Generate the R script
        script = self.r_generator.generate_figure_script(
            figure_name=figure_name,
            description=description,
            data_sources=data_sources,
            plots=[plot_config]
        )

        # Save the script
        script_path = self.r_generator.save_script(script, f"figure_{figure_name}")
        self.logger.info(f"Generated R script: {script_path}")

        # Optionally run the script
        if run_script:
            self.logger.info(f"Executing generated R script...")
            result = run_r_script(script_path)
            if result.success:
                self.logger.info(f"Script executed successfully in {result.duration_seconds:.1f}s")
            else:
                self.logger.error(f"Script execution failed: {result.errors}")

        return script, script_path

    def generate_treatment_comparison_figure(
        self,
        response_var: str,
        response_label: str,
        data_source: str = "cafi",
        plot_type: str = "boxplot",
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate an R script for a treatment comparison figure.

        This is a convenience method for the common pattern of comparing
        response variables across treatment groups.

        Args:
            response_var: Name of the response variable
            response_label: Human-readable label for the y-axis
            data_source: Data source to use
            plot_type: Type of plot
            run_script: Whether to execute the script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        return self.generate_figure_r_script(
            figure_name=f"{response_var}_by_treatment",
            description=f"Treatment comparison of {response_label}",
            data_sources=[data_source, "treatment"],
            x_var="treatment",
            y_var=response_var,
            plot_type=plot_type,
            fill_var="treatment",
            title=f"{response_label} by Treatment",
            x_label="Treatment (# Corals)",
            y_label=response_label,
            run_script=run_script
        )

    def generate_multi_panel_figure(
        self,
        figure_name: str,
        panels: List[Dict[str, Any]],
        ncol: int = 2,
        run_script: bool = False
    ) -> Tuple[str, Path]:
        """
        Generate an R script for a multi-panel figure using patchwork.

        Args:
            figure_name: Name for the figure
            panels: List of panel configurations, each with keys:
                   - name, x_var, y_var, plot_type, title, etc.
            ncol: Number of columns in the layout
            run_script: Whether to execute the script

        Returns:
            Tuple of (R script content, path to saved script)
        """
        from textwrap import dedent

        # Build header
        script = RCodeTemplates.script_header(
            title=f"Multi-panel Figure: {figure_name}",
            description=f"Multi-panel figure with {len(panels)} panels",
            inputs=["CAFI data", "Treatment data"],
            outputs=[f"output/MRB/figures/agent_generated/{figure_name}.png/pdf"]
        )

        script += RCodeTemplates.load_data_block(["cafi", "treatment"])

        # Data prep
        script += dedent('''
        # ==============================================================================
        # DATA PREPARATION
        # ==============================================================================

        analysis_data <- cafi_data %>%
          left_join(treatment_data, by = "coral_id") %>%
          filter(!is.na(treatment))

        ''')

        # Generate each panel
        script += "# ==============================================================================\n"
        script += "# GENERATE PANELS\n"
        script += "# ==============================================================================\n\n"

        panel_names = []
        for i, panel in enumerate(panels):
            panel_name = panel.get("name", f"panel_{i+1}")
            panel_names.append(panel_name)

            script += RCodeTemplates.ggplot_figure(
                plot_name=panel_name,
                data_var="analysis_data",
                x_var=panel.get("x_var", "treatment"),
                y_var=panel.get("y_var", "value"),
                plot_type=panel.get("plot_type", "boxplot"),
                color_var=panel.get("color_var"),
                fill_var=panel.get("fill_var", "treatment"),
                title=panel.get("title", ""),
                x_label=panel.get("x_label", ""),
                y_label=panel.get("y_label", "")
            )

        # Combine with patchwork
        script += dedent(f'''
        # ==============================================================================
        # COMBINE PANELS WITH PATCHWORK
        # ==============================================================================

        library(patchwork)

        combined_figure <- {" + ".join(panel_names)} +
          plot_layout(ncol = {ncol}) +
          plot_annotation(
            title = "{figure_name}",
            theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
          )

        # Save combined figure
        save_both(
          combined_figure,
          here::here("output", "MRB", "figures", "agent_generated", "{figure_name}"),
          width = {ncol * 4},
          height = {(len(panels) // ncol + 1) * 4}
        )

        message("Multi-panel figure complete: {figure_name}")
        ''')

        # Save and optionally run
        script_path = self.r_generator.save_script(script, f"figure_{figure_name}")

        if run_script:
            result = run_r_script(script_path)
            if not result.success:
                self.logger.error(f"Script execution failed: {result.errors}")

        return script, script_path

    def list_generated_scripts(self) -> List[Path]:
        """List all agent-generated R scripts."""
        generated_dir = SCRIPTS_DIR / "agent_generated"
        if not generated_dir.exists():
            return []
        return list(generated_dir.glob("*.R"))

    # -------------------------------------------------------------------------
    # FIGURE QUALITY CHECKS
    # -------------------------------------------------------------------------

    def check_figure_quality(self, figure_path: Path) -> Dict[str, Any]:
        """
        Check quality metrics of a figure.

        Args:
            figure_path: Path to the figure file

        Returns:
            Quality check results
        """
        if not figure_path.exists():
            return {"error": f"File not found: {figure_path}"}

        info = self._get_figure_info(figure_path)

        # Basic checks
        checks = {
            "file_exists": info.exists,
            "file_size_ok": info.size_bytes > 10000,  # At least 10KB
            "format": info.format.value if info.format else None,
        }

        # Try to get image dimensions (requires PIL)
        try:
            from PIL import Image
            with Image.open(figure_path) as img:
                width, height = img.size
                checks["width_px"] = width
                checks["height_px"] = height
                checks["resolution_ok"] = width >= 1000 and height >= 800

                # Check if likely publication quality
                if info.format == FigureFormat.PNG:
                    # Estimate DPI from typical print size (7 inches)
                    estimated_dpi = width / 7
                    checks["estimated_dpi"] = int(estimated_dpi)
                    checks["publication_ready"] = estimated_dpi >= 300

        except ImportError:
            checks["note"] = "Install Pillow for image dimension checks"
        except Exception as e:
            checks["dimension_error"] = str(e)

        return checks

    def batch_quality_check(self) -> Dict[str, Any]:
        """Check quality of all publication figures."""
        results = {
            "timestamp": datetime.now().isoformat(),
            "figures": [],
            "summary": {
                "total": 0,
                "passed": 0,
                "failed": 0,
                "missing": 0
            }
        }

        for config in PUBLICATION_FIGURES:
            check = self.check_figure_quality(config.output_path)
            check["figure_number"] = config.figure_number
            check["name"] = config.name

            results["figures"].append(check)
            results["summary"]["total"] += 1

            if "error" in check:
                results["summary"]["missing"] += 1
            elif check.get("publication_ready", True):
                results["summary"]["passed"] += 1
            else:
                results["summary"]["failed"] += 1

        return results

    # -------------------------------------------------------------------------
    # FIGURE EXPORT
    # -------------------------------------------------------------------------

    def export_figure(
        self,
        source_path: Path,
        output_dir: Path,
        formats: List[FigureFormat],
        resize: Optional[Tuple[int, int]] = None
    ) -> List[Path]:
        """
        Export a figure to multiple formats.

        Args:
            source_path: Path to source figure
            output_dir: Output directory
            formats: List of output formats
            resize: Optional (width, height) to resize

        Returns:
            List of exported file paths
        """
        if not source_path.exists():
            self.logger.error(f"Source not found: {source_path}")
            return []

        output_dir.mkdir(parents=True, exist_ok=True)
        exported = []

        base_name = source_path.stem

        for fmt in formats:
            output_path = output_dir / f"{base_name}.{fmt.value}"

            if fmt == FigureFormat.PNG and source_path.suffix.lower() == ".png":
                # Direct copy
                shutil.copy2(source_path, output_path)
                exported.append(output_path)

            elif fmt == FigureFormat.PDF and source_path.suffix.lower() == ".pdf":
                shutil.copy2(source_path, output_path)
                exported.append(output_path)

            else:
                # Try conversion with PIL
                try:
                    from PIL import Image
                    with Image.open(source_path) as img:
                        if resize:
                            img = img.resize(resize, Image.Resampling.LANCZOS)

                        if fmt == FigureFormat.PNG:
                            img.save(output_path, "PNG", dpi=(300, 300))
                        elif fmt == FigureFormat.TIFF:
                            img.save(output_path, "TIFF", dpi=(300, 300))
                        # PDF requires additional handling

                        exported.append(output_path)

                except ImportError:
                    self.logger.warning("Install Pillow for format conversion")
                except Exception as e:
                    self.logger.error(f"Conversion error: {e}")

        return exported

    def create_figure_archive(self, output_path: Optional[Path] = None) -> Path:
        """
        Create a ZIP archive of all publication figures.

        Returns:
            Path to the created archive
        """
        import zipfile

        output_path = output_path or (OUTPUT_DIR / "publication_figures.zip")

        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zf:
            for config in PUBLICATION_FIGURES:
                # Add PNG
                png_path = config.output_path.with_suffix(".png")
                if png_path.exists():
                    zf.write(png_path, png_path.name)

                # Add PDF
                pdf_path = config.output_path.with_suffix(".pdf")
                if pdf_path.exists():
                    zf.write(pdf_path, pdf_path.name)

        self.logger.info(f"Created archive: {output_path}")
        return output_path

    # -------------------------------------------------------------------------
    # HELPER METHODS
    # -------------------------------------------------------------------------

    def _get_figure_info(self, path: Path) -> FigureInfo:
        """Get information about a figure file."""
        if not path.exists():
            return FigureInfo(
                path=path,
                exists=False,
                format=self._get_format(path)
            )

        stat = path.stat()

        return FigureInfo(
            path=path,
            exists=True,
            format=self._get_format(path),
            size_bytes=stat.st_size,
            modified=datetime.fromtimestamp(stat.st_mtime)
        )

    def _get_format(self, path: Path) -> Optional[FigureFormat]:
        """Get figure format from path."""
        suffix = path.suffix.lower().lstrip(".")
        try:
            return FigureFormat(suffix)
        except ValueError:
            return None


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    """Command-line interface for the figure agent."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CAFI MRB Figure Generation Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python figure_agent.py --status           Show status of all figures
  python figure_agent.py --generate 1       Generate Figure 1
  python figure_agent.py --generate-all     Generate all figures
  python figure_agent.py --quality-check    Check quality of all figures
  python figure_agent.py --create-archive   Create ZIP of all figures
  python figure_agent.py --create-style     Create R styling script
        """
    )

    parser.add_argument(
        "--status", action="store_true",
        help="Show status of all publication figures"
    )
    parser.add_argument(
        "--generate", type=int,
        help="Generate a specific figure by number"
    )
    parser.add_argument(
        "--generate-all", action="store_true",
        help="Generate all publication figures"
    )
    parser.add_argument(
        "--quality-check", action="store_true",
        help="Run quality checks on all figures"
    )
    parser.add_argument(
        "--create-archive", action="store_true",
        help="Create ZIP archive of all figures"
    )
    parser.add_argument(
        "--create-style", action="store_true",
        help="Create R styling script"
    )
    parser.add_argument(
        "--scan", action="store_true",
        help="Scan and inventory all figures"
    )

    args = parser.parse_args()

    agent = FigureAgent()

    if args.status:
        print("\nPublication Figure Status:")
        print("=" * 70)
        statuses = agent.get_all_figure_status()
        for status in statuses:
            exists_mark = "OK" if status.get("exists") else "MISSING"
            print(f"  Figure {status['figure_number']}: [{exists_mark}] {status['name']}")
            if status.get("exists"):
                print(f"           Size: {status['size']}, Modified: {status['modified']}")
        return

    if args.generate:
        result = agent.generate_figure(args.generate)
        print(f"\nFigure {args.generate}: {'SUCCESS' if result.success else 'FAILED'}")
        if result.output_paths:
            for path in result.output_paths:
                print(f"  Output: {path}")
        if result.errors:
            for error in result.errors:
                print(f"  Error: {error}")
        return

    if args.generate_all:
        results = agent.generate_all_figures()
        print("\nFigure Generation Results:")
        print("=" * 50)
        for result in results:
            status = "OK" if result.success else "FAIL"
            print(f"  [{status}] {result.figure_name}")
        return

    if args.quality_check:
        results = agent.batch_quality_check()
        print("\nFigure Quality Check:")
        print("=" * 50)
        print(f"  Total: {results['summary']['total']}")
        print(f"  Passed: {results['summary']['passed']}")
        print(f"  Failed: {results['summary']['failed']}")
        print(f"  Missing: {results['summary']['missing']}")
        return

    if args.create_archive:
        path = agent.create_figure_archive()
        print(f"\nArchive created: {path}")
        return

    if args.create_style:
        path = agent.create_style_script()
        print(f"\nStyle script created: {path}")
        return

    if args.scan:
        inventory = agent.scan_figures()
        print("\nFigure Inventory:")
        print("=" * 50)
        for category, figures in inventory.items():
            print(f"\n{category}: {len(figures)} figures")
            for fig in figures[:5]:  # Show first 5
                print(f"  - {fig.path.name}")
            if len(figures) > 5:
                print(f"  ... and {len(figures) - 5} more")
        return

    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
