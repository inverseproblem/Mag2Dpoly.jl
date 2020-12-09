

########################################################################

"""
$(TYPEDSIGNATURES)

Total magnetic field (2D) for a set of polygonal bodies defined by their corners. Takes into account both induced and remnant magnetization.
 Based on Talwani & Heitzler (1964), the default algorithm in Mag2Dpoly package. 
"""
function tmagpolybodies2D(xzobs::Array{<:Real,2},Jinds::Vector{MagnetizVector},Jrems::Vector{MagnetizVector},
                    northxax::Real,bodies::MagPolygBodies2D)

    nbody = length(bodies.bo)
    tmag = zeros(eltype(Jinds[1].mod),size(xzobs,1))
    for i=1:nbody
        tmag .+= tmagpoly2D(xzobs,Jinds[i],Jrems[i],northxax,bodies.bo[i]) 
    end
    return tmag
end
###################################################################################

"""
$(TYPEDSIGNATURES)

Total magnetic field (2D) for a set of polygonal bodies defined by their corners. Takes into account both induced and remnant magnetization.
Generic version containing four different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani & Heitzler (1964)
  - "talwani_red"  --> Talwani & Heitzler (1964) rederived from Kravchinsky et al. 2019
  - "krav"         --> Kravchinsky et al. (2019) rectified by Ghirotto et al. (2020)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tmagpolybodies2Dgen(xzobs::Array{<:Real,2},Jinds::Vector{MagnetizVector},Jrems::Vector{MagnetizVector},
                       northxax::Real,bodies::MagPolygBodies2D,forwardtype::String)

    nbody = length(bodies.bo)
    tmag = zeros(eltype(Jinds[1].mod),size(xzobs,1))
    for i=1:nbody
        tmag .+= tmagpoly2Dgen(xzobs,Jinds[i],Jrems[i],northxax,bodies.bo[i],forwardtype)
    end
    return tmag
end


###################################################################################

"""
$(TYPEDSIGNATURES)

Total magnetic field (2D) for a polygon defined by its corners. Takes into account both induced and remnant magnetization.
Generic version containing four different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani & Heitzler (1964)
  - "talwani_red"  --> Talwani & Heitzler (1964) rederived from Kravchinsky et al. 2019
  - "krav"         --> Kravchinsky et al. (2019) rectified by Ghirotto et al. (2020)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tmagpoly2Dgen(xzobs::Array{<:Real,2},Jind::MagnetizVector,Jrem::MagnetizVector,
                       northxax::Real,body::BodySegments2D,forwardtype::String)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)

    if aclockw==false
        error("tmagforwardpoly2D(): vertices *not* ordered anticlockwise. Aborting.")
    end
   
    ##---------------------
    ## Get the angles
   
    ## check modules
    @assert Jind.mod >= 0.0
    @assert Jrem.mod >= 0.0

    ## `northxax` is the angle between geographic north and the positive x axis
    @assert 0.0 <= northxax <= 360.0
    Cnorth = deg2rad(northxax)
    
    ## check angles
    @assert -90.0 <= Jind.Ideg <= 90.0
    @assert -90.0 <= Jrem.Ideg <= 90.0
    
    @assert -180.0 <= Jind.Ddeg <= 180.0
    @assert -180.0 <= Jrem.Ddeg <= 180.0

    # check right forwardtype
    if forwardtype != "talwani" && forwardtype != "talwani_red" && forwardtype != "krav" && forwardtype != "wonbev"
        error("tmagforwardpoly2D(): [forwardtype] must be 'talwani' or 'talwani_red' or 'krav' or 'wonbev'")
    end 
    
    # deg to rad
    Iind = deg2rad(Jind.Ideg)
    Dind = deg2rad(Jind.Ddeg)
    Irem = deg2rad(Jrem.Ideg)
    Drem = deg2rad(Jrem.Ddeg)

    # Calculation of Jx and Jz only for case != wonbev
    if forwardtype != "wonbev"
        Jtotx,_,Jtotz = magcomp(Jind.mod,Iind,Dind,Jrem.mod,Irem,Drem,Cnorth)
    end
    
    ##-------------------------------------
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    totfield = Vector{typeof(Jind.mod)}(undef,nobs)

    ## loop on observation points
    for iob=1:nobs      
        
        xo = xzobs[iob,1]
        zo = xzobs[iob,2]

        ## loop on segments
        tsum = 0.0
        for ise=1:body.nsegm

            x1 = body.ver1[ise,1]-xo
            z1 = body.ver1[ise,2]-zo
            x2 = body.ver2[ise,1]-xo
            z2 = body.ver2[ise,2]-zo

            if forwardtype == "talwani"
                tsum += tmagtalwani(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "talwani_red"
                tsum += tmagtalwanired(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "krav"
                tsum += tmagkrav(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "wonbev"
                tsum += tmagwonbev(x1,z1,x2,z2,Jind.mod,Jrem.mod,Iind,Dind,Irem,Drem,Cnorth)
                
            end            
         end                   
                          
        # convert field from A/m to nT
        totfield[iob] = convert_H_to_B_nT( tsum )
    end

    return totfield
end

########################################################################

"""
$(TYPEDSIGNATURES)

Total magnetic field (2D) for a polygon defined by its corners. Takes into account both induced and remnant magnetization.
 Based on Talwani & Heitzler (1964), the default algorithm in Mag2Dpoly package. 
"""
function tmagpoly2D(xzobs::Array{<:Real,2},Jind::MagnetizVector,Jrem::MagnetizVector,
                    northxax::Real,body::BodySegments2D)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)

    if aclockw==false
        error("tmagforwardpoly2D(): vertices *not* ordered anticlockwise. Aborting.")
    end
   
    ##---------------------
    ## Get the angles
   
    ## check modules
    @assert Jind.mod >= 0.0
    @assert Jrem.mod >= 0.0

    ## `northxax` is the angle between geographic north and the positive x axis
    @assert 0.0 <= northxax <= 360.0
    Cnorth = deg2rad(northxax)
    
    ## check angles
    @assert -90.0 <= Jind.Ideg <= 90.0
    @assert -90.0 <= Jrem.Ideg <= 90.0
    
    @assert -180.0 <= Jind.Ddeg <= 180.0
    @assert -180.0 <= Jrem.Ddeg <= 180.0

    # deg to rad
    Iind = deg2rad(Jind.Ideg)
    Dind = deg2rad(Jind.Ddeg)
    Irem = deg2rad(Jrem.Ideg)
    Drem = deg2rad(Jrem.Ddeg)

    # Calculation of Jx and Jz only for case != wonbev
    Jtotx,_,Jtotz = magcomp(Jind.mod,Iind,Dind,Jrem.mod,Irem,Drem,Cnorth)

    ##-------------------------------------
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    totfield = Vector{typeof(Jind.mod)}(undef,nobs)

    ## loop on observation points
    for iob=1:nobs      
        
        xo = xzobs[iob,1]
        zo = xzobs[iob,2]

        ## loop on segments
        tsum = 0.0
        for ise=1:body.nsegm

            x1 = body.ver1[ise,1]-xo
            z1 = body.ver1[ise,2]-zo
            x2 = body.ver2[ise,1]-xo
            z2 = body.ver2[ise,2]-zo

            tsum += tmagtalwani(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)
         end                   
                          
        # convert field from A/m to nT
        totfield[iob] = convert_H_to_B_nT( tsum )
    end

    return totfield
end

#######################################################################################
"""
$(TYPEDSIGNATURES)

Check whether the polygonal body has segments ordered anticlockwise.
"""
function checkanticlockwiseorder(body::BodySegments2D)
    ## https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
    ## https://www.element84.com/blog/determining-the-winding-of-a-polygon-given-as-a-set-of-ordered-points
    #
    # Check direction (anti)clockwise for a reference
    #   system like the following:
    #
    #   z
    #  /\ 
    #  |      2
    #  |   1     3
    #  |      4
    #  |
    #  -------------> x
    #
    encarea2=0.0
    for ise=1:body.nsegm
        x1 = body.ver1[ise,1]
        z1 = body.ver1[ise,2] 
        x2 = body.ver2[ise,1]
        z2 = body.ver2[ise,2]
        encarea2 += (x2-x1)*(z2+z1)
    end
    #
    # anticlockwise -> encarea2 < 0.0
    # clockwise -> encarea2 > 0.0
    if encarea2<0.0
        anticlockw=true
    else
        anticlockw=false
    end
    #
    # The reference system for the magnetic anomaly functions
    #   is reversed in z:
    #
    #  -------------> x
    #  |
    #  |      4
    #  |   1     3
    #  |      2
    #  \/
    #  z
    #
    # so, consequently, we flip the direction of
    # clockwise/anticlockwise:   !(anticlockw)   
    return !anticlockw    
end

#####################################################################################

"""
$(TYPEDSIGNATURES)
 
 Total magnetic field (2D) for a line segment. Formulas from Talwani & Heitzler (1964).
"""
function tmagtalwani(x1::Real,z1::Real,x2::Real,z2::Real,
                     Jx::Real,Jz::Real,Iind::Real,Dind::Real,C::Real)
    

    # Quantities for error definitions
    small = 1e4*eps(typeof(Jx))
    anglelim = 0.995*π
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1
    s = sqrt(x21^2+z21^2)
    
    # Return 0 if two corners are too close
    if s < small
        return 0.0
    end

    #-----------------------
    
    # Get the angles
    θ1 = atan(z1,x1)
    θ2 = atan(z2,x2)
           
    # If z21 = 0.0 no contribution    
    if z21 != 0.0
        g = -x21/z21
    else
        return 0.0
    end

    ϕ = acot(g) 

    θdiff = θ2-θ1
    # In the case polygon sides cross the x axis    
    if θdiff < -π
        θdiff = θdiff + 2.0*π
    elseif θdiff > π
        θdiff = θdiff - 2.0*π
    end    

    #--------------------------------------

    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if x1 < small && z1 < small
        x1 = small
        z1 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if x2 < small && z2 < small
        x2 = small
        z2 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    r1 = sqrt(x1^2+z1^2)
    r2 = sqrt(x2^2+z2^2)
    
    flog = log(r2)-log(r1)
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #-------------------------------------   

    # vertical component
    V = 2.0*sin(ϕ) * (Jx*( (θdiff)*cos(ϕ) + sin(ϕ)*flog) -
                      Jz*( (θdiff)*sin(ϕ) - cos(ϕ)*flog) )

    # horizonatal component
    H = 2.0*sin(ϕ) * (Jx*( (θdiff)*sin(ϕ) - cos(ϕ)*flog)+
                      Jz*( (θdiff)*cos(ϕ) + sin(ϕ)*flog) )

    # Divided by 4π to take into account algorithm formulation in emu units 
    totfield = (1.0/(4.0*π)) * (H*cos(Iind)*cos(C-Dind) + V*sin(Iind))
    
    return totfield 
end

###########################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D) for a line segment. Formulas from Kravchinsky et al (2019) rectified by Ghirotto et al. (2021). 
"""
function tmagkrav(x1::Real,z1::Real,x2::Real,z2::Real,
                  Jtotx::Real,Jtotz::Real,Iind::Real,Dind::Real,Cnorth::Real)


    # Quantities for errors definitions
    small = 1e4*eps(typeof(Jtotx))
    anglelim = 0.995*π
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1
    tmpγ = sqrt(x21^2+z21^2)

    # Return 0 if two corners are too close
    if tmpγ < small
        return 0.0
    end

    #--------------------
    
    # check if den ≠ 0.0
    if tmpγ!=0.0
        γx = x21 / tmpγ
        γz = z21 / tmpγ
    else
        return 0.0
    end

    # if the segment is horizontal it provides no contribution!
    if z21==0.0
        return 0.0
    end
    
    #------------
    g = x21/z21

    if x1 >= g*z1  
        δ = 1.0
    elseif x1 < g*z1
        δ = -1.0
    end

    #--------------------
    # Get the angles
    α1 = atan(δ*(z1+g*x1),(x1-g*z1))
    α2 = atan(δ*(z2+g*x2),(x2-g*z2))
    

    #In the case polygon sides cross the x axis
    αdiff = α2 - α1
    
    if αdiff < -π
        αdiff = αdiff + 2.0*π
    elseif αdiff > π
        αdiff = αdiff - 2.0*π
    end 

    #----------------------------------

    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if x1 < small && z1 < small
        x1 = small
        z1 = small
        @warn  "A corner is too close to an observation point (calculation continues)"
    end
    
    if x2 < small && z2 < small
        x2 = small
        z2 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end

    r1 = sqrt(x1^2+z1^2)
    r2 = sqrt(x2^2+z2^2)
   
    lor21 = log(r2)-log(r1)
    
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(αdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
           
    #------------------------------------
    
    P = γz*γx*lor21 + δ*(γz^2)*(αdiff)
    Q = (γz^2)*lor21 - δ*γx*γz*(αdiff)
    
    ## horizonatl and vertical field components
    H = 1.0/(2.0*pi) * (Jtotz*Q + Jtotx*P)
    V = 1.0/(2.0*pi) * (Jtotx*Q - Jtotz*P)
    
    ## total field anomaly 
    totfield = V*sin(Iind)+H*cos(Iind)*cos(Cnorth-Dind)

    return totfield
end

########################################################################################

""" 
$(TYPEDSIGNATURES)

 Total magnetic field (2D) for a ribbon. Talwani & Heitzler (1964) modified by Kravchinsky et al. (2019).
"""
function tmagtalwanired(x1::Real,z1::Real,x2::Real,z2::Real,
                        Jx::Real,Jz::Real,Iind::Real,Dind::Real,C::Real)


    # Quantities for errors definitions
    small = 1e4*eps(typeof(Jx))
    anglelim = 0.995*π
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1
    s = sqrt(x21^2+z21^2)
    
    # Return 0 if two corners are too close
    if s < small
        return 0.0
    end

    #-----------------------
    
    # if the segment is horizontal it provides no contribution!
    if z21 != 0.0
        g = -x21/z21
    else
        return 0.0
    end    
        
    ϕ = acot(g)
    
    den1 = x1+z1*cot(ϕ)
    den2 = x2+z2*cot(ϕ)
    num1 = z1-x1*cot(ϕ)
    num2 = z2-x2*cot(ϕ)

    # Controls on signs of atan argument (abs in den1 and den2)
    #-----------------------
    if den1 < 0.0
        den1 = -den1
        δ = -1.0
        θ1 = atan(num1,den1)
    else
        δ = 1.0
        θ1 = atan(num1,den1)
    end 

    if den2 < 0.0
        den2 = -den2
        θ2 = atan(num2,den2)
    else
        θ2 = atan(num2,den2)        
    end
    #-----------------------

    # In the case polygon sides cross the x axis
    θdiff = θ2-θ1
    
    if θdiff < -π
        θdiff = θdiff + 2.0*π
    elseif θdiff > π
        θdiff = θdiff - 2.0*π
    end    

    #-------------------------------------------------

    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if x1 < small && z1 < small
        x1 = small
        z1 = small
        @warn "A corner is too close to an observation point (calculation continues)\n"
    end
    
    if x2 < small && z2 < small
        x2 = small
        z2 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end

    r1 = sqrt(x1^2+z1^2)
    r2 = sqrt(x2^2+z2^2)

    flog = log(r2)-log(r1)
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #----------------------------------------------------

    # vertical component    
    V = 2.0*sin(ϕ) * (Jx * (δ*(θdiff)*cos(ϕ) + sin(ϕ)*flog)-
                    Jz * (δ*(θdiff)*sin(ϕ) - cos(ϕ)*flog) )
 
    # horizontal component
    H = 2.0*sin(ϕ) * (Jx * (δ*(θdiff)*sin(ϕ) - cos(ϕ)*flog)+
                    Jz * (δ*(θdiff)*cos(ϕ) + sin(ϕ)*flog) )
    
    ## total field anomaly divided by 4π to take into account algorithm formulation in emu units
    totfield = (1.0/(4.0*π)) * (H*cos(Iind)*cos(C-Dind) + V*sin(Iind))
    
    return totfield  
end


#######################################################################################################################

""" 
$(TYPEDSIGNATURES)

 Total magnetic field (2D) for a line segment. Formulas from Won & Bevis (1987).
"""
function tmagwonbev(x1::Real,z1::Real,x2::Real,z2::Real,
                    modJind::Real,modJrem::Real,Iind::Real,Dind::Real,
                    Irem::Real,Drem::Real,C::Real)

    # β is angle among North and profle direction
    βi = Dind - C + π/2
    βr = Drem - C + π/2
    
    #-------------------

    # Quantities for errors definitions
    small = 1e4*eps(typeof(modJrem))
    anglelim = 0.995*π

    #----------------------
    x21 = x2-x1
    z21 = z2-z1
    R = sqrt(x21^2+z21^2)

    # Return 0 if two corners are too close
    if R < small
        return 0.0
    end

    #------------------------

    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if x1 < small && z1 < small
        x1 = small
        z1 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if x2 < small && z2 < small
        x2 = small
        z2 = small
        @warn "A corner is too close to an observation point (calculation continues)"
    end

    r1 = sqrt(x1^2+z1^2)
    r2 = sqrt(x2^2+z2^2)

    lor21 = log(r2) - log(r1)

    θ1 = atan(z1,x1) 
    θ2 = atan(z2,x2)

    # In the case polygon sides cross the x axis
    if sign(z1) != sign(z2)
        test = x1*z2 - x2*z1
        if test > 0.0
            if z1 >= 0.0 
                θ2 = θ2 + 2π
            end
        elseif test < 0.0
            if z2 >= 0.0
                θ1 = θ1 + 2π
            end
        else
            return 0.0 
        end
    end
    
    # Error if the side is too close to the observation point (calculation continues)
    θdiff = θ1-θ2
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #-----------------------
    
    P = (1/R^2)*(x1*z2 - x2*z1)*(((x1*x21 - z1*z21)/(r1^2))-
                                 ((x2*x21 - z2*z21)/(r2^2)))

    Q = (1/R^2)*(x1*z2 - x2*z1)*(((x1*z21 + z1*x21)/(r1^2))-
                                 ((x2*z21 + z2*x21)/(r2^2)))
    
    if x21 != 0.0
        
        g = z21/x21
        derZz = ((x21^2)/(R^2))*((θdiff) + g*lor21) - P
        derZx = -((x21*z21)/(R^2))*((θdiff) + g*lor21) + Q
        derXz = -((x21^2)/(R^2))*(g*(θdiff) - lor21) + Q
        derXx = ((x21*z21)/(R^2))*(g*(θdiff) - lor21) + P
    
    else

        derZz = -P
        derZx = -((z21^2)/(R^2))*lor21 + Q
        derXz = Q
        derXx = ((z21^2)/(R^2))*(θdiff) + P
        
    end

    # Magnetic strenght components due to induced magnetization
    ΔHzind = 2.0*modJind*(sin(Iind)*derZz + sin(βi)*cos(Iind)*derZx) 
    ΔHxind = 2.0*modJind*(sin(Iind)*derXz + sin(βi)*cos(Iind)*derXx) 

    # Magnetic strenght components due to remnant magnetization
    ΔHzrem = 2.0*modJrem*(sin(Irem)*derZz + sin(βr)*cos(Irem)*derZx) 
    ΔHxrem = 2.0*modJrem*(sin(Irem)*derXz + sin(βr)*cos(Irem)*derXx) 

    ΔHztot = ΔHzind + ΔHzrem
    ΔHxtot = ΔHxind + ΔHxrem

    ## total field anomaly divided by 4π to take into account algorithm formulation in emu units
    ΔHtot = -(1.0/(4.0*π))*(ΔHztot*sin(Iind) + ΔHxtot*sin(βi)*cos(Iind))
   
    return ΔHtot
end

###################################################################################
