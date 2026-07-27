"""Quick inspection helper."""
import inspect

from bioforge.omics import batch, cluster, normalize, qc, trajectory

for m in [
    qc.run_qc,
    normalize.run_normalize,
    cluster.run_cluster,
    trajectory.run_trajectory,
    batch.run_harmony,
]:
    print(m.__module__ + "." + m.__name__, inspect.signature(m))
