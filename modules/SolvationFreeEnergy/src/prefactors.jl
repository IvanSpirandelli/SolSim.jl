struct WhiteBearParameters
    rs::Float64
    η::Float64
end

function get_prefactors_from_white_bear(rs::Float64, η::Float64)
    @SVector[
    pressure(rs, η), # Volume
    sigma(rs, η), # Area
    kappa(rs, η), # Mean curvature 
    kappa_bar(η), # Gaussian curvature
    ]
end

""" Prefactors.
    Given in equation.24 of Density functional theory for hard-sphere mixtures:
    the White Bear version mark II.
"""

function pressure(rs::Float64, η::Float64)
    fraction = (1+η+η^2-η^3)/(1-η)^3
    scaling = η*3/4/pi/rs^3
    scaling * fraction
end

function sigma(rs::Float64, η::Float64)
    fraction_1 = (1+2*η + 8*η^2 - 5*η^3)/(3*(1-η)^3)
    fraction_2 = (log(1-η)/(3*η))
    scaling = -η*3/4/pi/rs^2 #Note the MINUS
    scaling * (fraction_1 + fraction_2)
end

function kappa(rs::Float64, η::Float64)
    fraction_1 = (4-10*η+20*η^2-8*η^3)/(3*(1-η)^3)
    fraction_2 = (4*log(1-η)/(3*η))
    scaling = η*3/4/pi/rs
    scaling * (fraction_1 + fraction_2)
end

function kappa_bar(η::Float64)
    fraction_1 = (-4+11*η-13*η^2+4*η^3)/(3*(1-η)^3)
    fraction_2 = (4*log(1-η))/(3*η)
    scaling = η*3/4/pi
    scaling * (fraction_1 - fraction_2)
end
