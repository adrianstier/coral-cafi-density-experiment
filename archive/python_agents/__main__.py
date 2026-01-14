"""
CLI entry point for CAFI MRB Analysis Agents.

Usage:
    python -m agents --help
    python -m agents --status
    python -m agents --run-full
"""

from orchestrator import main

if __name__ == "__main__":
    main()
