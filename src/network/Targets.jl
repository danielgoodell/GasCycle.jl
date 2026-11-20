# ── User-declared free variables and targets ─────────────────────────────────
#
# NPSS-style solver pairs: declare model inputs the solver may vary
# (`add_free!`) and model outputs it must drive to a value (`add_target!`).
# They join the element unknowns/residuals and back-edge states in one Newton
# vector, in both design and off-design solves. The counts need not pair up
# one-to-one, but the total system must be square.

_auto_name(net, kind, n) = "$kind $(n + 1)"

function _check_bounds(lower, upper, name)
    lower < upper || error("add_free! \"$name\": lower ($lower) must be < upper ($upper)")
end

# Write a free-variable value into an element field. Fields typed `Real` (the
# solver-owned off-design fields) accept anything; concretely typed design
# fields (e.g. HeatSource{Float64}.TtExit) accept plain numbers but cannot hold
# the Dual numbers of an outer gradient.
function _set_field!(el, field::Symbol, x, name::String)
    FT = fieldtype(typeof(el), field)
    if x isa FT
        setfield!(el, field, x)
    elseif !(x isa ForwardDiff.Dual)
        setfield!(el, field, convert(FT, x))
    else
        error("""free variable "$name": cannot store a ForwardDiff Dual in field \
               `$field::$FT` of element "$(el.name)". When differentiating through a \
               solve, construct the element with this parameter in the gradient's number \
               type, e.g. `$(nameof(typeof(el)))(...; $field = 1.0 * one(eltype(p)))`.""")
    end
    nothing
end

"""
    add_free!(net, el, field; name, ref, lower, upper)

Declare element field `el.field` as a free variable: the solver varies it to
satisfy the targets declared with [`add_target!`](@ref). `el` may be any element
in the network, a `Shaft` added with `add_shaft!`, or an electrical node.

    add_free!(net, field; ...)          # seed-state field: :W, :Tt or :Pt
    add_free!(net, name; get, set!, ...) # anything else, via accessor closures

Keywords:
- `name`  — label used in messages (default `"<element>.<field>"`).
- `ref`   — normalization scale; default `|initial value|` (or 1 if zero).
  The solver works in `x / ref`, so pick a typical magnitude of the variable.
- `lower`, `upper` — box bounds; trial values are clamped into `[lower, upper]`.

The current value is the initial guess. Parameters the element recomputes
itself (e.g. `HeatSource.Q` in `:fixed_TtExit` mode) are not meaningful free
variables. When a gradient flows through the solve, a concretely typed field
must already hold the gradient's number type (see the error message).
"""
function add_free!(net::FlowNetwork, el::Union{AbstractElement,AbstractElectrical}, field::Symbol;
                   name::AbstractString = "$(el.name).$field",
                   ref::Real = NaN, lower::Real = -Inf, upper::Real = Inf)
    (any(e -> e === el, net.elements) || any(s -> s === el, net.shafts) ||
     any(e -> e === el, net.electrical)) ||
        error("add_free!: element \"$(el.name)\" is not in the network — call add! / add_shaft! / connect_electrical! first")
    hasfield(typeof(el), field) ||
        error("add_free!: $(nameof(typeof(el))) \"$(el.name)\" has no field `$field`")
    getfield(el, field) isa Real ||
        error("add_free!: field `$field` of \"$(el.name)\" is not numeric")
    _check_bounds(lower, upper, name)
    nm = String(name)
    push!(net.free_vars, FreeVar(nm, () -> getfield(el, field),
                                 x -> _set_field!(el, field, x, nm),
                                 Float64(ref), Float64(lower), Float64(upper)))
    net
end

function add_free!(net::FlowNetwork, field::Symbol;
                   name::AbstractString = "seed.$field",
                   ref::Real = NaN, lower::Real = -Inf, upper::Real = Inf)
    field in (:W, :Tt, :Pt) ||
        error("add_free!: seed-state free variable must be :W, :Tt or :Pt, got :$field")
    _check_bounds(lower, upper, name)
    get = () -> begin
        isnothing(net.seed_state) && error("add_free! \"$name\": call set_state! before solving")
        getfield(net.seed_state, field)
    end
    set! = x -> (net.seed_state = update(net.seed_state; (field => x,)...); nothing)
    push!(net.free_vars, FreeVar(String(name), get, set!,
                                 Float64(ref), Float64(lower), Float64(upper)))
    net
end

function add_free!(net::FlowNetwork, name::AbstractString;
                   get::Function, set!::Function,
                   ref::Real = NaN, lower::Real = -Inf, upper::Real = Inf)
    _check_bounds(lower, upper, name)
    push!(net.free_vars, FreeVar(String(name), get, x -> (set!(x); nothing),
                                 Float64(ref), Float64(lower), Float64(upper)))
    net
end

"""
    add_target!(net, f, value; name, scale)

Declare a target: the solver drives `f(net)` to `value`. `f` is any function of
the network evaluated after each pass, e.g.

    add_target!(net, net_power, 10e3)                      # 10 kW net shaft power
    add_target!(net, n -> comp.inlet[].Tt, 400.0)          # compressor inlet 400 K

The residual is `(f(net) − value) / scale`; `scale` defaults to `|value|`
(or 1 when `value` is zero) — set it explicitly for targets near zero.
`value` may carry ForwardDiff Duals (differentiate w.r.t. a requirement).
"""
function add_target!(net::FlowNetwork, f::Function, value::Real;
                     name::AbstractString = _auto_name(net, "target", length(net.targets)),
                     scale::Real = NaN)
    sc = isnan(scale) ? (iszero(value) ? 1.0 : abs(Float64(_primal(value)))) : Float64(scale)
    sc > 0 || error("add_target! \"$name\": scale must be positive")
    push!(net.targets, Target(String(name), f, value, sc))
    net
end

"""
    clear_targets!(net)

Remove all free variables and targets declared with `add_free!`/`add_target!`.
"""
function clear_targets!(net::FlowNetwork)
    empty!(net.free_vars)
    empty!(net.targets)
    net
end

# Strip all Dual levels down to the plain number.
_primal(x::Real) = x
_primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))
