"""Thread-level parallelism for independent work items.

numpy and scipy release the GIL inside FFTs, sorts/partitions and vector ufuncs, so
independent hypotheses (carrier offsets, raster channels, receiver packages) overlap
well on threads without the memory cost of processes. Results come back in input
order, so a parallel run gives exactly what the sequential run gives.

SPECTRUM_INTEL_THREADS sets the thread count (1 = strictly sequential); the default is
min(4, CPU count), which leaves headroom on an 8 GB / 4-big-core edge box.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor


def threads() -> int:
    v = os.environ.get("SPECTRUM_INTEL_THREADS")
    return max(1, int(v)) if v else min(4, os.cpu_count() or 1)


def pmap(fn, items, n_threads: int | None = None) -> list:
    """[fn(x) for x in items], on up to n_threads threads, in input order."""
    items = list(items)
    n = min(n_threads or threads(), len(items))
    if n <= 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=n) as ex:
        return list(ex.map(fn, items))
