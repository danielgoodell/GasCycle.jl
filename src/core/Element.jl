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
