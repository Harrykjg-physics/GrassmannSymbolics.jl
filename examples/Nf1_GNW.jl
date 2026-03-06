
@grassmann η1 ξ1 η2 ξ2 η1b ξ1b η2b ξ2b
@grassmann ψ1 ψ2 ψ1b ψ2b

@variables μ

term1 = exp(0.5 * sqrt(2) * (ψ1b - ψ2b) * η1)
term2 = exp(0.5 * sqrt(2) * η1b * (ψ1 - ψ2))
term3 = exp(0.5 * sqrt(2) * exp(0.5 * μ) * (ψ1b - im * ψ2b)) * η2

