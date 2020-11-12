
#########################################################################

"""
$(TYPEDSIGNATURES)

Vector addition of magnetic (remnant + induced) components.
"""
function magcomp(modJind::Real,Iind::Real,Dind::Real,modJrem::Real,Irem::Real,
                 Drem::Real,C::Real)
    
    ## Induced magnetization components
    Jix = modJind*cos(Iind)*cos(C-Dind)
    Jiy = modJind*cos(Iind)*sin(C-Dind)
    Jiz = modJind*sin(Iind)

    ## Remnant magnetization components
    Jrx = modJrem*cos(Irem)*cos(C-Drem)
    Jry = modJrem*cos(Irem)*sin(C-Drem)
    Jrz = modJrem*sin(Irem)

    ## Vector addition    
    Jtotx = Jix+Jrx
    Jtoty = Jiy+Jry
    Jtotz = Jiz+Jrz
   
    return Jtotx,Jtoty,Jtotz
end

##############################################

"""
$(TYPEDSIGNATURES)

Convert from the field H (A/m) to B (nT).
"""
function convert_H_to_B_nT( H_Am::Real ) 
    ## permeabilita' del vuoto 
    ## muzero = 4.0 * pi * 10.0^-7
    ## B nanoTesla
    ## B_nT = ( muzero * H_Am ) * 10.0^9
    B_nT =  pi * 400.0 * H_Am 
    return B_nT
end

###############################################

"""
$(TYPEDSIGNATURES)

Convert from the field B (nT) to H (A/m).
"""
function convert_B_nT_to_H( B_nT::Real ) 
    H_Am = B_nT / (pi * 400.0)
    return H_Am
end

###############################################
