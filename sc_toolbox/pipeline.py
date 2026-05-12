from dataclasses import dataclass, field
from typing import Protocol

import anndata as ad


class RunFn(Protocol):
    def __call__(self, adata: ad.AnnData, params: dict) -> ad.AnnData: ...


class RenderFn(Protocol):
    def __call__(self, adata: ad.AnnData, defaults: dict) -> dict: ...


@dataclass
class Step:
    name: str
    label: str
    auto: bool
    run: RunFn
    render: RenderFn | None = None
    needs: list[str] = field(default_factory=lambda: ["counts_matrix"])


@dataclass
class Pipeline:
    name: str
    steps: list[Step]
    enabled: dict[str, bool] = field(default_factory=dict)

    def active_steps(self) -> list[Step]:
        return [s for s in self.steps if self.enabled.get(s.name, True)]

    def aggregated_needs(self) -> set[str]:
        out: set[str] = set()
        for s in self.active_steps():
            out.update(s.needs)
        return out
