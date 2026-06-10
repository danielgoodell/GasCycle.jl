"""
Abstract base type for all cycle elements.

Every concrete element must implement:

  compute!(el, inlet::Port) -> Port
      Propagate thermodynamic state: read inlet, write and return outlet.
      Called once per solver iteration.

  residuals(el) -> Vector{Float64}
      Return the element's constraint violations (should be ≈ 0 at convergence).
      Called once per solver iteration to build the global residual vector.

  n_residuals(el) -> Int
      Number of residual equations this element contributes.
      Must be constant for a given element instance.

  indep_vars(el) -> Vector{Float64}
      Current values of the element's independent (free) variables.
      The solver adjusts these during iteration.

  set_indep_vars!(el, x::AbstractVector)
      Write a new set of independent variable values into the element.
"""
abstract type AbstractElement end

"""
    _polytropic_outlet(fp, Tt_in, Pt_in, Pt_out, η_p; N=20) -> Tt_out

Numerically integrate the polytropic process from `Pt_in` to `Pt_out` in `N`
equal log-pressure steps.  Each small stage applies the isentropic formula
within that step, scaled by the polytropic efficiency:

  Compression (Pt_out > Pt_in):  Δh_actual = Δh_is / η_p   (more work than ideal)
  Expansion   (Pt_out < Pt_in):  Δh_actual = Δh_is × η_p   (less work than ideal)

20 steps gives better than 0.01 % accuracy for He-Xe across typical CBC
pressure ratios.  All fluid calls are generic so ForwardDiff Dual numbers
propagate through correctly when the fluid backend supports it (IdealGasFluid
has closed-form inversions; FPTFluid uses bisection and does not support AD).
"""
function _polytropic_outlet(fp::FluidProperties, Tt_in, Pt_in, Pt_out, η_p; N=20)
    compress = Pt_out > Pt_in   # true → compressor, false → turbine/expander
    T = Tt_in
    P = Pt_in
    ln_step = log(Pt_out / Pt_in) / N   # positive for compression, negative for expansion
    for _ in 1:N
        P_next = P * exp(ln_step)
        s_here = entropy(fp, T, P)
        h_here = enthalpy(fp, T, P)
        T_is   = T_from_s(fp, s_here, P_next; T_guess = T)
        h_is   = enthalpy(fp, T_is, P_next)
        Δh_is  = h_is - h_here
        Δh     = compress ? Δh_is / η_p : Δh_is * η_p
        T      = T_from_h(fp, h_here + Δh, P_next; T_guess = T_is)
        P      = P_next
    end
    T
end

"""
    _isentropic_outlet(fp, Tt_in, Pt_in, Pt_out, η_ad) -> Tt_out

Outlet temperature using adiabatic (isentropic) efficiency semantics — the
single-step counterpart of `_polytropic_outlet`, matching NPSS `effDes`:

  Compression:  Δh_actual = Δh_is / η_ad
  Expansion:    Δh_actual = Δh_is × η_ad

with Δh_is evaluated for the full pressure change at constant inlet entropy.
"""
function _isentropic_outlet(fp::FluidProperties, Tt_in, Pt_in, Pt_out, η_ad)
    compress = Pt_out > Pt_in
    s_in  = entropy(fp, Tt_in, Pt_in)
    h_in  = enthalpy(fp, Tt_in, Pt_in)
    T_is  = T_from_s(fp, s_in, Pt_out; T_guess = Tt_in)
    Δh_is = enthalpy(fp, T_is, Pt_out) - h_in
    Δh    = compress ? Δh_is / η_ad : Δh_is * η_ad
    T_from_h(fp, h_in + Δh, Pt_out; T_guess = T_is)
end

"""Dispatch on an element's efficiency semantics (`:polytropic` or `:isentropic`)."""
_efficiency_outlet(η_type::Symbol, fp::FluidProperties, Tt_in, Pt_in, Pt_out, η) =
    η_type == :isentropic ? _isentropic_outlet(fp, Tt_in, Pt_in, Pt_out, η) :
                            _polytropic_outlet(fp, Tt_in, Pt_in, Pt_out, η)

function compute!(el::AbstractElement, inlet::Port)
    error("compute! not implemented for $(typeof(el))")
end

function residuals(el::AbstractElement)
    error("residuals not implemented for $(typeof(el))")
end

function n_residuals(el::AbstractElement)::Int
    error("n_residuals not implemented for $(typeof(el))")
end

function indep_vars(el::AbstractElement)
    error("indep_vars not implemented for $(typeof(el))")
end

function set_indep_vars!(el::AbstractElement, x::AbstractVector)
    error("set_indep_vars! not implemented for $(typeof(el))")
end
