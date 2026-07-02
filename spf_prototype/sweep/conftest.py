"""Put the sweep package dir (and the existing spf_prototype/python dir) on
sys.path so tests can `import coherence_metrics`, `import quasr_loader`, and reuse
the existing pipeline modules (quasr_geom, run_conformal, ...)."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent          # spf_prototype/sweep
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "python"))  # spf_prototype/python
