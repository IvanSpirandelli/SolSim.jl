struct WhiteBearParameters
    rs::Float64
    η::Float64
end

function get_prefactors(rs, eta; id = "wb", r = 1.0)
    if id == "wb"
        return get_wb_prefactors(rs, eta)
    elseif id == "virial"
        return get_virial_prefactors(rs, eta)
    elseif id == "amp"
        return get_amped_prefactors(r, rs, eta)
    end
    @assert false
end

# These prefactors are equivalent to doubling the negative contributions to the mean curvature 
# for a union of hard spheres of radius r. (This was a bug in a previous version of AlphaMol) 
function get_amped_prefactors(r::Float64, rs::Float64, η::Float64)
    pf_wb = get_prefactors(rs, η)
    @SVector[pf_wb[1], pf_wb[2] - (pf_wb[3]/(r + rs)), 2.0*pf_wb[3], pf_wb[4]]
end

""" Prefactors.
    Given in equation.24 of Density functional theory for hard-sphere mixtures:
    the White Bear version mark II.
"""
function get_wb_prefactors(rs::Float64, η::Float64)
    @SVector[ 
    pressure(rs, η), # Volume
    sigma(rs, η), # Area
    kappa(rs, η), # Mean curvature 
    kappa_bar(η), # Gaussian curvature
    ]
end

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

function get_virial_prefactors(rs::Float64, η::Float64) 
    @SVector[
        pressure(rs, η), # Volume
        virial_sigma(rs, η), # Area
        virial_kappa(rs, η), # Mean curvature 
        virial_kappa_bar(η), # Gaussian curvature
    ]
end

function virial_sigma(rs, η)
    scaling = η*3/4/pi/rs^2 
    n1 = 18*η*((1+η+η^2-η^3)/(1-η)^3)
    n2 = sqrt(3)*pi*η*(1-16*η-4*η^2+7η^3)/(1-η)^3
    n3 = (12-sqrt(3)*pi)*log(1-η)
    n4 = 12*log((2-η)/(2*(1-η)^3))
    nominator = n1 + n2 - n3 - n4
    denominator = 9*(sqrt(3)*pi-4)*η
    scaling * nominator/denominator
end

function virial_kappa(rs, η)
    scaling = -η*3/4/pi/rs #Note the MINUS
    n1 = 2*η * (11-23*η+η^2+5*η^3)/(1-η)^3
    n2 = 8*log(1-η)
    n3 = 12*log((2-η)/(2*(1-η)^3))
    nominator = n1 - n2 - n3
    denominator = 3*(sqrt(3)*pi-4)*η
    scaling * nominator/denominator
end

function virial_kappa_bar(η)
    scaling = η*3/4/pi
    n1 = 12*η*(5-12*η+3*η^3)/(1-η)^3
    n2 = sqrt(3)*pi*η*(4-13*η-η^2+4*η^3)/(1-η)^3
    n3 = 4*sqrt(3)*pi*log(1-η)
    n4 = 24*log((2-η)/(2*(1-η)^3))
    nominator = n1 - n2 - n3 - n4
    denominator = 9*(sqrt(3)*pi-4)*η
    scaling * nominator/denominator
end