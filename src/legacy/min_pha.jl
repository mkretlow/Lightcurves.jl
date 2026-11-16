"""
## Purpose
Find date of minimum solar phase angle of small bodies (Asteroids, TNOs, Centaurs, etc.).

## functions
   min_pha : calculate table with minimum solar phase angle of small bodies (e.g. for opposition effect)
"""

export min_pha

rad = ac.rad
deg = ac.deg

struct Eph
    tjd :: Real
    ra  :: Real
    de  :: Real
    r   :: Real
    Δ   :: Real
    α   :: Real
    elg :: Real
    mag :: Real
    Eph() = new()
    Eph(tjd,ra,de,r,Δ,α,elg,mag) = new(tjd,ra,de,r,Δ,α,elg,mag)
end

# Orbital elements in rad
function ephem_twobody_classic(pv_earth_t1t2,teph1,teph2,epoch,mano0,sma,ecc,peri,node,incl,absmag,slope=0.15,dteph=1)

    pvh = zeros(6)
    pve = zeros(6)

    P,Q = gauss_vector(peri, node, incl)

    eph = Eph[]

    mdm = ac.gk0/sma^1.5 # daily motion in rad

    mass = 0

    # Generate ephemeris for time span teph1 -> teph2
    i = 0
    for t = teph1:dteph:teph2

        i+=1

        # Mean anomaly for t
        mano = mano0 + mdm*(t-epoch)

        if ecc < 1.0
            gk  = ac.gk0 * sqrt(1.0 + mass)
            ea  = ADTK.exanom(mano,ecc)
            S,C = sincos(ea)
            fac = sqrt(1.0 - ecc^2)
            x   = sma * (C - ecc)
            y   = sma * fac * S
            vx  = -gk * S / sqrt(sma) / (1.0 - C*ecc)
            vy  = +gk * C / sqrt(sma) * fac / (1.0 - C*ecc)

        else
            return push!(eph,Eph(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
        end

        pvh[1:3] = P*x + Q*y
        pvh[4:6] = P*vx + Q*vy

        # Heliocentric distance
        r = norm(pvh[1:3])

        # Get cartesian ecliptical earth coordinates.
        pv_earth = pv_earth_t1t2[i]

        # Transform to geocenter
        pve = pvh - pv_earth

        Δ = norm((pve[1:3]))

        # LT correction
        tau  = Δ / ac.clight
        pve[1:3] = pve[1:3] - tau * pvh[4:6]
        pvh[1:3] = pvh[1:3] - r/ac.clight * pvh[4:6]
        r = norm(pvh[1:3])

        # PV liegt nun im System der Bahnelemente (i.a. ekliptikal) vor. Trafo in das Äquatorsystem.
        #obl = eeps(1976,equin)
        obl = rad * ac.eps00_76
        pve = so3pv(1,-obl,pve)

        # Equatorial coordinates
        ra,de,Δ = rec2sph(pve[1:3])

        # Apparent magnitude, phase angle and elongation
        vmag,phw,elg = ADTK.ephmag(0,slope,absmag,Δ,r,norm(pv_earth[1:3])^2)
        
        # Save values in Array{Eph()}
        push!(eph,Eph(t,ra/15*deg,de*deg,r,Δ,deg*phw,deg*elg,vmag))
    end

    return eph
end


"""
   `min_pha(dat1,dat2)`

Find date of minimum solar phase angle of Small Bodies (Asteroids, TNOs, Centaurs, etc.)

## Args

   dat1 : start date as number in format yyyymmdd(.d)
   dat2 : end date as number in format yyyymmdd(.d)

## Return
   Prints table to stdout

"""
function min_pha(dat1 = 2023_01_01.0, dat2 = 2024_01_01.0; phalim = 0.1, maglim = 22.0, asteledb = "/data1/astele/astele.sdb")

    if !isfile(asteledb) 
        @error "Cannot open DB at $(asteledb). Path correct?"
        return
    end

    eph1 = ydate(dat1)
    eph2 = ydate(dat2)
    dteph = 0.1

    db1 = SQLite.DB(asteledb)    

    df1 = DataFrame(DBInterface.execute(db1, "SELECT id,rdes,iepoch,ma0,sma,ecc,peri,node,incl,h10,slope,diam,arcspan FROM astorb WHERE sma > 30.0 LIMIT 10000"));

    # Earth coordinates for whole eph time span
    pv_earth = []

    for tt in eph1:dteph:eph2
        pv,ierr = vsop87(tt, 1, 3, 0.0)
        push!(pv_earth,pv)
    end

    N = 0
    for i in 1:size(df1)[1]

        if df1[i,:ecc] >= 1
            @printf("%-20s : Object skipped because e >= 1\n", df1[i,:rdes])
        end

        eph = ephem_twobody_classic(pv_earth,eph1,eph2,ydate(df1[i,:iepoch]),rad*df1[i,:ma0],df1[i,:sma],df1[i,:ecc], rad*df1[i,:peri],rad*df1[i,:node],rad*df1[i,:incl], df1[i,:h10], df1[i,:slope],dteph)

        # Find Minimum in phase angle. minimum(x->x.α, eph) # liefert aber nicht den Index
        minpha,idx = findmin(getfield.(eph, :α))
        if minpha < phalim && eph[idx].elg > 30.0 && eph[idx].mag < maglim
        
            # Elongation to Moon and Moon phase for rejecting unobservable events
            k = mphase(eph[idx].tjd)
            
            k > 0.75 && continue
            
            ra_m,de_m,ρ = minimoon(eph[idx].tjd)
            elm = acos(sin(rad*eph[idx].de)*sin(de_m)+cos(rad*eph[idx].de)*cos(de_m)*cos(15*rad*eph[idx].ra-ra_m))

            @printf("%-20s | %s | %0.1f | Pha = %0.3f° | Elg = %6.2f° | Mag = %0.1f | RA = %04.1fh  DE = %+03d° | Elg.Moon = %03d°  k =%3d%% | a=%0.1f e=%0.2f i=%04.1f°\n", 
            df1[i,:rdes], outdat(caldat(eph[idx].tjd),10), eph[idx].tjd, eph[idx].α, eph[idx].elg, eph[idx].mag, eph[idx].ra, eph[idx].de, deg*elm, 100*k, 
            df1[i,:sma], df1[i,:ecc], df1[i,:incl])
            N += 1
        end
    end

    println("\n$N objects with a > 30 au, pha < $(phalim)°, elongation > 30° and mag < $(maglim) found")
end

