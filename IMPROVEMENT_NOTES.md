# Improvement Notes

These are follow-up improvements identified after the initial code review fixes.
The package is already fast enough for current use, so these are mainly about
preventing future performance regressions and keeping the graph/solver layer
manageable as the model grows.

## 1. Add a Benchmark Suite

Status: implemented in `benchmarks/runbenchmarks.jl`.

Add steady-state benchmarks for representative workflows:

- simple open Brayton cycle
- recuperated closed-loop design solve
- off-design map solve
- ForwardDiff sensitivity/gradient solve

The benchmark suite should avoid measuring first-run compilation as much as
practical, and should be easy to run before and after graph/solver changes.

## 2. Refactor `FlowNetwork.one_pass!`

Status: first pass implemented. `FlowNetwork` now caches an internal flow plan
with required inlets, back-edges, and forward edges grouped by source element.
`one_pass!` uses a ready queue instead of repeatedly scanning all elements and
all edges.

`one_pass!` currently works, but it does a lot of mutable bookkeeping:

- `Dict{Tuple{Int,Symbol}, Any}` for port availability
- repeated scans over all elements until traversal stalls
- repeated scans over all edges after each element compute
- concrete-element special cases in central traversal code

A future cleanup should precompute a compiled network plan:

- element execution order
- required inlet ports
- downstream edges grouped by source element
- back-edge list
- stable port indices instead of symbol/dict lookups inside each residual call

This is the most likely place for future speed and readability gains.

## 3. Introduce an Element Port Interface

Status: implemented as a private network-port interface in `FlowNetwork.jl`.
The central traversal now asks each element for required inlets, outlet ports,
inlet aliases, outlet values, and its network compute hook instead of branching
on each special element type in `_compute_element!`.

`FlowNetwork` currently knows too much about `HeatExchanger`, `Splitter`, and
`Mixer` internals. Prefer element-local methods such as:

```julia
required_inlets(el)
get_outlet(el, port)
set_inlet!(el, port, value)
compute_network_element!(el)
outlet_ports(el)
```

That would make new elements easier to add without editing central traversal
logic.

## 4. Split `Solver.jl`

Status: implemented. `Solver.jl` is now the public coordinator and includes
focused files for result handling, residual assembly, design solve,
off-design solve, and cycle metrics.

`Solver.jl` currently handles independent-variable collection, residual
assembly, design solving, off-design solving, result display, and cycle metrics.
Consider splitting it into smaller files:

- `SolveResult.jl`
- `ResidualAssembly.jl`
- `DesignSolve.jl`
- `OffDesignSolve.jl`
- `CycleMetrics.jl`

## 5. Make Map Bounds Behavior Explicit

Status: implemented. `PerformanceMap` now accepts `bounds=:error` by default,
with `:warn` and `:clamp` available for intentional clamping. Scaled maps
inherit the base map bounds policy unless overridden.

`FPTFluid` now has explicit bounds behavior, but `PerformanceMap.query` still
silently clamps corrected speed and corrected flow. Consider the same style:

- `bounds=:error`
- `bounds=:warn`
- `bounds=:clamp`

Silent map clamping can hide surge/choke/off-map behavior.

## 6. Add Network Validation

Add a `validate(net)` pass before solving to catch common setup mistakes:

- duplicated element names
- unreachable elements
- missing required inlet connections
- unknown/residual count mismatches
- unsupported fluid mixing
- heat-exchanger hot/cold registration issues
- explicit map/table bounds policy

This would make solver failures easier to diagnose.
