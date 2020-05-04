function Helmert_X(X, Y, Z, DX, Y_rot, Z_rot, S) {
    /*Computed Helmert transformed X coordinate.
        Input: - _
        cartesian XYZ coords (X,Y,Z), X translation (DX) all in meters ; _
        Y and Z rotations in seconds of arc (Y_Rot, Z_Rot) and scale in ppm (s).

        Convert rotations to radians and ppm scale to a factor
    */
    var sfactor = S * 0.000001;
    var RadY_Rot = (Y_rot / 3600) * (Math.PI / 180);
    var RadZ_Rot = (Z_rot / 3600) * (Math.PI / 180);
    return X + (X * sfactor) - (Y * RadZ_Rot) + (Z * RadY_Rot) + DX;
}

function Helmert_Y(X, Y, Z, DY, X_rot, Z_rot, S) {
    /*    
        Computed Helmert transformed Y coordinate.
        Input: - _
        cartesian XYZ coords(X, Y, Z), Y translation(DY) all in meters; _
        X and Z rotations in seconds of arc(X_Rot, Z_Rot) and scale in ppm(s).
        
        Convert rotations to radians and ppm scale to a factor
    */
    var sfactor = S * 0.000001;
    var RadX_Rot = (X_rot / 3600) * (Math.PI / 180);
    var RadZ_Rot = (Z_rot / 3600) * (Math.PI / 180);

    //Compute transformed Y coord
    return (X * RadZ_Rot) + Y + (Y * sfactor) - (Z * RadX_Rot) + DY;
}




function Helmert_Z(X, Y, Z, DZ, X_rot, Y_rot, S) {
    /*
        Computed Helmert transformed Z coordinate.
        Input: - _
        cartesian XYZ coords (X,Y,Z), Z translation (DZ) all in meters ; _
        X and Y rotations in seconds of arc (X_Rot, Y_Rot) and scale in ppm (s).
        
        Convert rotations to radians and ppm scale to a factor
    */

    var sfactor = S * 0.000001;
    var RadX_Rot = (X_rot / 3600) * (Math.PI / 180);
    var RadY_Rot = (Y_rot / 3600) * (Math.PI / 180);
    //Compute transformed Z coord
    return  (-1 * X * RadY_Rot) + (Y * RadX_Rot) + Z + (Z * sfactor) + DZ;
}



function XYZ_to_Lat(X, Y, Z, a, b) {
    /*
        Convert XYZ to Latitude (PHI) in Dec Degrees.
        Input: - _
        XYZ cartesian coords (X,Y,Z) and ellipsoid axis dimensions (a & b), all in meters.
        
        THIS FUNCTION REQUIRES THE "Iterate_XYZ_to_Lat" FUNCTION
        THIS FUNCTION IS CALLED BY THE "XYZ_to_H" FUNCTION
    */

    var RootXYSqr = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2));
    var e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
    var PHI1 = Math.atan(Z / (RootXYSqr * (1 - e2)));

    var PHI = Iterate_XYZ_to_Lat(a, e2, PHI1, Z, RootXYSqr);

    return PHI * (180 / Math.PI);
}



function Iterate_XYZ_to_Lat(a, e2, PHI1, Z, RootXYSqr) {
    /*
        Iteratively computes Latitude (PHI).
        Input: - _
            ellipsoid semi major axis (a) in meters; _
            eta squared (e2); _
            estimated value for latitude (PHI1) in radians; _
            cartesian Z coordinate (Z) in meters; _
            RootXYSqr computed from X & Y in meters.
            
        THIS FUNCTION IS CALLED BY THE "XYZ_to_PHI" FUNCTION
        THIS FUNCTION IS ALSO USED ON IT'S OWN IN THE _
        "Projection and Transformation Calculations.xls" SPREADSHEET
            
    */

    var V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHI1)), 2))));
    var PHI2 = Math.atan((Z + (e2 * V * (Math.sin(PHI1)))) / RootXYSqr);

    while (Math.abs(PHI1 - PHI2) > 0.000000001) {
        PHI1 = PHI2;
        V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHI1)), 2))));
        PHI2 = Math.atan((Z + (e2 * V * (Math.sin(PHI1)))) / RootXYSqr);
    }

    return PHI2;
}


function XYZ_to_Long(X, Y) {
    /*
        Convert XYZ to Longitude (LAM) in Dec Degrees.
        Input: - _
        X and Y cartesian coords in meters.
    */


    //tailor the output to fit the equatorial quadrant as determined by the signs of X and Y

    if (X >= 0) //longitude is in the W90 thru 0 to E90 hemisphere
        return (Math.atan(Y / X)) * (180 / Math.PI);

    if (X < 0 && Y >= 0) //longitude is in the E90 to E180 quadrant
        return ((Math.atan(Y / X)) * (180 / Math.PI)) + 180;

    if (X < 0 && Y < 0)//longitude is in the E180 to W90 quadrant
        return ((Math.atan(Y / X)) * (180 / Math.PI)) - 180;
}



function XYZ_to_H(X, Y, Z, a, b) {
    /*
        Convert XYZ to Ellipsoidal Height.
        Input: - _
        XYZ cartesian coords (X,Y,Z) and ellipsoid axis dimensions (a & b), all in meters.

        REQUIRES THE "XYZ_to_Lat" FUNCTION
    */

    //Compute PHI (Dec Degrees) first
    var PHI = XYZ_to_Lat(X, Y, Z, a, b);

    //Convert PHI radians
    var RadPHI = PHI * (Math.PI / 180);

    //Compute H
    var RootXYSqr = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2));
    var e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
    var V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));
    var H = (RootXYSqr / Math.cos(RadPHI)) - V;

    return H;
}



function Lat_Long_H_to_X(PHI, LAM, H, a, b) {
    /*
        Convert geodetic coords lat (PHI), long (LAM) and height (H) to cartesian X coordinate.
        Input: - _
        Latitude (PHI)& Longitude (LAM) both in decimal degrees; _
        Ellipsoidal height (H) and ellipsoid axis dimensions (a & b) all in meters.
    */

    //Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180);
    var RadLAM = LAM * (Math.PI / 180);

    //Compute eccentricity squared and nu
    var e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
    var V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));

    //Compute X
    return (V + H) * (Math.cos(RadPHI)) * (Math.cos(RadLAM));
}




function Lat_Long_H_to_Y(PHI, LAM, H, a, b) {
    /*
        Convert geodetic coords lat (PHI), long (LAM) and height (H) to cartesian Y coordinate.
        Input: - _
            Latitude (PHI)& Longitude (LAM) both in decimal degrees; _
            Ellipsoidal height (H) and ellipsoid axis dimensions (a & b) all in meters.
    */

    //'Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180);
    var RadLAM = LAM * (Math.PI / 180);

    //Compute eccentricity squared and nu
    var e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
    var V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));

    //Compute Y
    return (V + H) * (Math.cos(RadPHI)) * (Math.sin(RadLAM));
}



function Lat_H_to_Z(PHI, H, a, b) {
    /*
        Convert geodetic coord components latitude (PHI) and height (H) to cartesian Z coordinate.
        Input: - _
            Latitude (PHI) decimal degrees; _
            Ellipsoidal height (H) and ellipsoid axis dimensions (a & b) all in meters.
    */

    //Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180);

    //Compute eccentricity squared and nu
    var e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
    var V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));

    //Compute X
    return ((V * (1 - e2)) + H) * (Math.sin(RadPHI));
}





function Lat_Long_to_East(PHI, LAM, a, b, E0, f0, PHI0, LAM0) {
    /*
    Project Latitude and longitude to Transverse Mercator eastings.
    Input: - _
     Latitude (PHI) and Longitude (LAM) in decimal degrees; _
     ellipsoid axis dimensions (a & b) in meters; _
     eastings of false origin (e0) in meters; _
     central meridian scale factor (f0); _
     latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
    */

    //Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180);
    var RadLAM = LAM * (Math.PI / 180);
    var RadPHI0 = PHI0 * (Math.PI / 180);
    var RadLAM0 = LAM0 * (Math.PI / 180);
    var af0 = a * f0;
    var bf0 = b * f0;
    var e2 = (Math.pow(af0, 2) - Math.pow(bf0, 2)) / Math.pow(af0, 2);
    var N = (af0 - bf0) / (af0 + bf0);
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(RadPHI)), 2));
    var eta2 = (nu / rho) - 1;
    var p = RadLAM - RadLAM0;

    var IV = nu * (Math.cos(RadPHI));
    var V = (nu / 6) * Math.pow((Math.cos(RadPHI)), 3) * ((nu / rho) - (Math.pow(Math.tan(RadPHI), 2)));
    var VI = (nu / 120) * Math.pow((Math.cos(RadPHI)), 5) * (5 - (18 * Math.pow((Math.tan(RadPHI)), 2)) + Math.pow((Math.tan(RadPHI)), 4) + (14 * eta2) - (58 * Math.pow((Math.tan(RadPHI)), 2) * eta2));

    return E0 + (p * IV) + (Math.pow(p, 3) * V) + (Math.pow(p, 5) * VI);

}




function Lat_Long_to_North(PHI, LAM, a, b, E0, N0, f0, PHI0, LAM0) {
    /*
        Project Latitude and longitude to Transverse Mercator northings
        'Input: - _
        Latitude (PHI) and Longitude (LAM) in decimal degrees; _
        ellipsoid axis dimensions (a & b) in meters; _
        eastings (e0) and northings (n0) of false origin in meters; _
        central meridian scale factor (f0); _
        latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
        
        'REQUIRES THE "Marc" FUNCTION
    */

    //Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180);
    var RadLAM = LAM * (Math.PI / 180);
    var RadPHI0 = PHI0 * (Math.PI / 180);
    var RadLAM0 = LAM0 * (Math.PI / 180);

    var af0 = a * f0;
    var bf0 = b * f0;
    var e2 = (Math.pow(af0, 2) - Math.pow(bf0, 2)) / Math.pow(af0, 2);
    var N = (af0 - bf0) / (af0 + bf0);
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(RadPHI)), 2));
    var eta2 = (nu / rho) - 1;
    var p = RadLAM - RadLAM0;
    var m = Marc(bf0, N, RadPHI0, RadPHI);

    var i = m + N0;
    var II = (nu / 2) * (Math.sin(RadPHI)) * (Math.cos(RadPHI));
    var III = ((nu / 24) * (Math.sin(RadPHI)) * Math.pow((Math.cos(RadPHI)), 3)) * (5 - Math.pow((Math.tan(RadPHI)), 2) + (9 * eta2));
    var IIIA = ((nu / 720) * (Math.sin(RadPHI)) * ((Math.cos(RadPHI)) /**/, 5)) * (61 - (58 * Math.pow((Math.tan(RadPHI)), 2)) + Math.pow((Math.tan(RadPHI)), 4));

    return i + (Math.pow(p /**/, 2) * II) + (Math.pow(p /**/, 4) * III) + (Math.pow(p /**/, 6) * IIIA);

}



function E_N_to_Lat(East, North, a, b, E0, N0, f0, PHI0, LAM0) {
    /*
        'Un-project Transverse Mercator eastings and northings back to latitude.
        'Input: - _
            eastings (East) and northings (North) in meters; _
            ellipsoid axis dimensions (a & b) in meters; _
            eastings (e0) and northings (n0) of false origin in meters; _
            central meridian scale factor (f0) and _
            latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
        
        'REQUIRES THE "Marc" AND "InitialLat" FUNCTIONS
    */

    //Convert angle measures to radians
    var RadPHI0 = PHI0 * (Math.PI / 180);
    var RadLAM0 = LAM0 * (Math.PI / 180);

    //Compute af0, bf0, e squared (e2), n and Et
    var af0 = a * f0;
    var bf0 = b * f0;
    var e2 = (Math.pow(af0 /**/, 2) - Math.pow(bf0 /**/, 2)) / Math.pow(af0 /**/, 2);
    var N = (af0 - bf0) / (af0 + bf0);
    var Et = East - E0;

    //Compute initial value for latitude (PHI) in radians
    var PHId = InitialLat(North, N0, af0, RadPHI0, N, bf0);

    //Compute nu, rho and eta2 using value for PHId
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)) /**/, 2))));
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) /**/, 2));
    var eta2 = (nu / rho) - 1;

    //Compute Latitude
    var VII = (Math.tan(PHId)) / (2 * rho * nu);
    var VIII = ((Math.tan(PHId)) / (24 * rho * Math.pow(nu /**/, 3))) * (5 + (3 * Math.pow((Math.tan(PHId)) /**/, 2)) + eta2 - (9 * eta2 * Math.pow((Math.tan(PHId)) /**/, 2)));
    var IX = ((Math.tan(PHId)) / (720 * rho * Math.pow(nu /**/, 5))) * (61 + (90 * Math.pow((Math.tan(PHId)) /**/, 2)) + (45 * Math.pow((Math.tan(PHId)) /**/, 4)));

    return (180 / Math.PI) * (PHId - (Math.pow(Et /**/, 2) * VII) + (Math.pow(Et /**/, 4) * VIII) - (Math.pow(Et /**/, 6) * IX));
}




function E_N_to_Long(East, North, a, b, E0, N0, f0, PHI0, LAM0) {
    /*
        Un-project Transverse Mercator eastings and northings back to longitude.
        Input: - _
            eastings (East) and northings (North) in meters; _
            ellipsoid axis dimensions (a & b) in meters; _
            eastings (e0) and northings (n0) of false origin in meters; _
            central meridian scale factor (f0) and _
            latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
    */

    //REQUIRES THE "Marc" AND "InitialLat" FUNCTIONS

    //Convert angle measures to radians
    var RadPHI0 = PHI0 * (Math.PI / 180);
    var RadLAM0 = LAM0 * (Math.PI / 180);

    //Compute af0, bf0, e squared (e2), n and Et
    var af0 = a * f0;
    var bf0 = b * f0;
    var e2 = (Math.pow(af0 /**/, 2) - Math.pow(bf0 /**/, 2)) / Math.pow(af0 /**/, 2);
    var N = (af0 - bf0) / (af0 + bf0)
    var Et = East - E0

    //Compute initial value for latitude (PHI) in radians
    var PHId = InitialLat(North, N0, af0, RadPHI0, N, bf0)

    //Compute nu, rho and eta2 using value for PHId
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) /**/, 2))
    var eta2 = (nu / rho) - 1

    //Compute Longitude
    var X = Math.pow((Math.cos(PHId)) /**/, -1) / nu
    var XI = (Math.pow((Math.cos(PHId)) /**/, -1) / (6 * Math.pow(nu /**/, 3))) * ((nu / rho) + (2 * Math.pow((Math.tan(PHId)) /**/, 2)))
    var XII = (Math.pow((Math.cos(PHId)) /**/, -1) / (120 * Math.pow(nu /**/, 5))) * (5 + (28 * Math.pow((Math.tan(PHId)) /**/, 2)) + (24 * Math.pow((Math.tan(PHId)) /**/, 4)))
    var XIIA = (Math.pow((Math.cos(PHId)) /**/, -1) / (5040 * Math.pow(nu /**/, 7))) * (61 + (662 * Math.pow((Math.tan(PHId)) /**/, 2)) + (1320 * Math.pow((Math.tan(PHId)) /**/, 4)) + (720 * Math.pow((Math.tan(PHId)) /**/, 6)))

    return (180 / Math.PI) * (RadLAM0 + (Et * X) - (Math.pow(Et /**/, 3) * XI) + (Math.pow(Et /**/, 5) * XII) - (Math.pow(Et /**/, 7) * XIIA))
}



function InitialLat(North, N0, afo, PHI0, N, bfo) {
    /*
    Compute initial value for Latitude (PHI) IN RADIANS.
    Input: - _
        northing of point (North) and northing of false origin (n0) in meters; _
        semi major axis multiplied by central meridian scale factor (af0) in meters; _
        latitude of false origin (PHI0) IN RADIANS; _
        n (computed from a, b and f0) and _
        ellipsoid semi major axis multiplied by central meridian scale factor (bf0) in meters.
    
    REQUIRES THE "Marc" FUNCTION
    THIS FUNCTION IS CALLED BY THE "E_N_to_Lat", "E_N_to_Long" and "E_N_to_C" FUNCTIONS
    THIS FUNCTION IS ALSO USED ON IT'S OWN IN THE  "Projection and Transformation Calculations.xls" SPREADSHEET
    */

    //First PHI value (PHI1)
    var PHI1 = ((North - N0) / afo) + PHI0

    //Calculate M
    var m = Marc(bfo, N, PHI0, PHI1)

    //Calculate new PHI value (PHI2)
    var PHI2 = ((North - N0 - m) / afo) + PHI1

    //Iterate to get final value for InitialLat
    while (Math.abs(North - N0 - m) > 0.00001) {
        PHI2 = ((North - N0 - m) / afo) + PHI1
        m = Marc(bfo, N, PHI0, PHI2)
        PHI1 = PHI2
    }

    return PHI2
}






function Marc(bf0, N, PHI0, PHI) {
    /*
        Compute meridional arc.
        Input: - _
        ellipsoid semi major axis multiplied by central meridian scale factor (bf0) in meters; _
        n (computed from a, b and f0); _
        lat of false origin (PHI0) and initial or final latitude of point (PHI) IN RADIANS.

        THIS FUNCTION IS CALLED BY THE - _
        "Lat_Long_to_North" and "InitialLat" FUNCTIONS
        THIS FUNCTION IS ALSO USED ON IT'S OWN IN THE "Projection and Transformation Calculations.xls" SPREADSHEET
    */

    return bf0 * (((1 + N + ((5 / 4) * Math.pow(N /**/, 2)) + ((5 / 4) * Math.pow(N /**/, 3))) * (PHI - PHI0))
        - (((3 * N) + (3 * Math.pow(N /**/, 2)) + ((21 / 8) * Math.pow(N /**/, 3))) * (Math.sin(PHI - PHI0)) * (Math.cos(PHI + PHI0)))
        + ((((15 / 8) * Math.pow(N /**/, 2)) + ((15 / 8) * Math.pow(N /**/, 3))) * (Math.sin(2 * (PHI - PHI0))) * (Math.cos(2 * (PHI + PHI0))))
        - (((35 / 24) * Math.pow(N /**/, 3)) * (Math.sin(3 * (PHI - PHI0))) * (Math.cos(3 * (PHI + PHI0)))));
}




function Lat_Long_to_C(PHI, LAM, LAM0, a, b, f0) {
    /*
    'Compute convergence (in decimal degrees) from latitude and longitude
    'Input: - _
     latitude (PHI), longitude (LAM) and longitude (LAM0) of false origin in decimal degrees; _
     ellipsoid axis dimensions (a & b) in meters; _
     central meridian scale factor (f0).
     */

    //'Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180)
    var RadLAM = LAM * (Math.PI / 180)
    var RadLAM0 = LAM0 * (Math.PI / 180)

    //'Compute af0, bf0 and e squared (e2)
    var af0 = a * f0
    var bf0 = b * f0
    var e2 = (Math.pow(af0 /**/, 2) - Math.pow(bf0 /**/, 2)) / Math.pow(af0 /**/, 2)

    //'Compute nu, rho, eta2 and p
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(RadPHI)) /**/, 2))
    eta2 = (nu / rho) - 1
    var p = RadLAM - RadLAM0

    //'Compute Convergence
    var XIII = Math.sin(RadPHI)
    var XIV = ((Math.sin(RadPHI) * Math.pow((Math.cos(RadPHI)) /**/, 2)) / 3) * (1 + (3 * eta2) + (2 * Math.pow(eta2 /**/, 2)))
    var XV = ((Math.sin(RadPHI) * Math.pow((Math.cos(RadPHI)) /**/, 4)) / 15) * (2 - Math.pow((Math.tan(RadPHI)) /**/, 2))

    return (180 / Math.PI) * ((p * XIII) + (Math.pow(p /**/, 3) * XIV) + (Math.pow(p /**/, 5) * XV))

}



function E_N_to_C(East, North, a, b, E0, N0, f0, PHI0) {
    /*
    'Compute convergence (in decimal degrees) from easting and northing
    'Input: - _
     Eastings (East) and Northings (North) in meters; _
     ellipsoid axis dimensions (a & b) in meters; _
     easting (e0) and northing (n0) of true origin in meters; _
     central meridian scale factor (f0); _
     latitude of central meridian (PHI0) in decimal degrees.
    
    'REQUIRES THE "Marc" AND "InitialLat" FUNCTIONS
    */

    //'Convert angle measures to radians
    var RadPHI0 = PHI0 * (Math.PI / 180)

    // 'Compute af0, bf0, e squared (e2), n and Et
    var af0 = a * f0
    var bf0 = b * f0
    var e2 = (Math.pow(af0 /**/, 2) - Math.pow(bf0 /**/, 2)) / Math.pow(af0 /**/, 2)
    var N = (af0 - bf0) / (af0 + bf0)
    var Et = East - E0

    // 'Compute initial value for latitude (PHI) in radians
    var PHId = InitialLat(North, N0, af0, RadPHI0, N, bf0)

    // 'Compute nu, rho and eta2 using value for PHId
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) /**/, 2))
    var eta2 = (nu / rho) - 1

    // 'Compute Convergence
    var XVI = (Math.tan(PHId)) / nu
    var XVII = ((Math.tan(PHId)) / (3 * Math.pow(nu /**/, 3))) * (1 + Math.pow((Math.tan(PHId)) /**/, 2) - eta2 - Math.pow(2 * (eta2 /**/, 2)))
    var XVIII = ((Math.tan(PHId)) / (15 * Math.pow(nu /**/, 5))) * (2 + (5 * Math.pow((Math.tan(PHId)) /**/, 2)) + (3 * Math.pow((Math.tan(PHId)) /**/, 4)))

    return (180 / Math.PI) * ((Et * XVI) - (Math.pow(Et /**/, 3) * XVII) + (Math.pow(Et /**/, 5) * XVIII))
}




function Lat_Long_to_LSF(PHI, LAM, LAM0, a, b, f0) {
    /*
    'Compute local scale factor from latitude and longitude
    'Input: - _
     latitude (PHI), longitude (LAM) and longitude (LAM0) of false origin in decimal degrees; _
     ellipsoid axis dimensions (a & b) in meters; _
     central meridian scale factor (f0).
     */

    //  'Convert angle measures to radians
    var RadPHI = PHI * (Math.PI / 180)
    var RadLAM = LAM * (Math.PI / 180)
    var RadLAM0 = LAM0 * (Math.PI / 180)

    // 'Compute af0, bf0 and e squared (e2)
    var af0 = a * f0
    var bf0 = b * f0
    var e2 = (Math.pow(af0 /**/, 2) - (bf0 /**/, 2)) / Math.pow(af0 /**/, 2)

    // 'Compute nu, rho, eta2 and p
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(RadPHI)) /**/, 2))
    var eta2 = (nu / rho) - 1
    var p = RadLAM - RadLAM0

    // 'Compute local scale factor
    var XIX = (Math.pow(Math.cos(RadPHI) /**/, 2) / 2) * (1 + eta2)
    var XX = (Math.pow(Math.cos(RadPHI) /**/, 4) / 24) * (5 - (4 * Math.pow((Math.tan(RadPHI)) /**/, 2)) + (14 * eta2) - (28 * Math.pow((Math.tan(RadPHI * eta2)) /**/, 2)))

    return f0 * (1 + (Math.pow(p /**/, 2) * XIX) + (Math.pow(p /**/, 4) * XX))

}



function E_N_to_LSF(East, North, a, b, E0, N0, f0, PHI0) {
    /*
        'Compute local scale factor from from easting and northing
        'Input: - _
        Eastings (East) and Northings (North) in meters; _
        ellipsoid axis dimensions (a & b) in meters; _
        easting (e0) and northing (n0) of true origin in meters; _
        central meridian scale factor (f0); _
        latitude of central meridian (PHI0) in decimal degrees.
        'REQUIRES THE "Marc" AND "InitialLat" FUNCTIONS
     */

    // 'Convert angle measures to radians
    var RadPHI0 = PHI0 * (Math.PI / 180)

    // 'Compute af0, bf0, e squared (e2), n and Et
    var af0 = a * f0
    var bf0 = b * f0
    var e2 = (Math.pow(af0 /**/, 2) - (bf0 /**/, 2)) / Math.pow(af0 /**/, 2)
    var N = (af0 - bf0) / (af0 + bf0)
    var Et = East - E0

    // 'Compute initial value for latitude (PHI) in radians
    var PHId = InitialLat(North, N0, af0, RadPHI0, N, bf0)

    // 'Compute nu, rho and eta2 using value for PHId
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) /**/, 2))
    var eta2 = (nu / rho) - 1

    // 'Compute local scale factor
    var XXI = 1 / (2 * rho * nu)
    var XXII = (1 + (4 * eta2)) / (24 * Math.pow(rho /**/, 2) * (nu /**/, 2))

    return f0 * (1 + (Math.pow(Et /**/, 2) * XXI) + (Math.pow(Et /**/, 4) * XXII))
}




function E_N_to_t_minus_T(AtEast, AtNorth, ToEast, ToNorth, a, b, E0, N0, f0, PHI0) {
    /*
    'Compute (t-T) correction in decimal degrees at point (AtEast, AtNorth) to point (ToEast,ToNorth)
    'Input: - _
     Eastings (AtEast) and Northings (AtNorth) in meters, of point where (t-T) is being computed; _
     Eastings (ToEast) and Northings (ToNorth) in meters, of point at other end of line to which (t-T) is being computed; _
     ellipsoid axis dimensions (a & b) and easting & northing (e0 & n0) of true origin in meters; _
     central meridian scale factor (f0); _
     latitude of central meridian (PHI0) in decimal degrees.
    
    'REQUIRES THE "Marc" AND "InitialLat" FUNCTIONS
    */

    // 'Convert angle measures to radians
    var RadPHI0 = PHI0 * (Math.PI / 180)

    // 'Compute af0, bf0, e squared (e2), n and Nm (Northing of mid point)
    var af0 = a * f0
    var bf0 = b * f0
    var e2 = (Math.pow(af0 /**/, 2) - Math.pow(bf0 /**/, 2)) / Math.pow(af0 /**/, 2)
    var N = (af0 - bf0) / (af0 + bf0)
    var Nm = (AtNorth + ToNorth) / 2

    // 'Compute initial value for latitude (PHI) in radians
    var PHId = InitialLat(Nm, N0, af0, RadPHI0, N, bf0)

    // 'Compute nu, rho and eta2 using value for PHId
    var nu = af0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)) /**/, 2))))
    var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) /**/, 2))

    // 'Compute (t-T)
    var XXIII = 1 / (6 * nu * rho)

    return (180 / Math.PI) * ((2 * (AtEast - E0)) + (ToEast - E0)) * (AtNorth - ToNorth) * XXIII

}




function TrueAzimuth(AtEast, AtNorth, ToEast, ToNorth, a, b, E0, N0, f0, PHI0) {
    /*
    'Compute true azimuth in decimal degrees at point (AtEast, AtNorth) to point (ToEast,ToNorth)
    'Input: - _
     Eastings (AtEast) and Northings (AtNorth) in meters, of point where true azimuth is being computed; _
     Eastings (ToEast) and Northings (ToNorth) in meters, of point at other end of line to which true azimuth is being computed; _
     ellipsoid axis dimensions (a & b) and easting & northing (e0 & n0) of true origin in meters; _
     central meridian scale factor (f0); _
     latitude of central meridian (PHI0) in decimal degrees.
    
    'REQUIRES THE "Marc", "InitialLat", "E_N_to_t_minus_T" and "E_N_to_C" FUNCTIONS
    */

    // 'Compute eastings and northings differences
    var Diffe = ToEast - AtEast
    var Diffn = ToNorth - AtNorth
    var GridBearing;

    // 'Compute grid bearing
    if (Diffe == 0) {
        if (Diffn < 0)
            GridBearing = 180
        else
            GridBearing = 0
    }
    else {

        var Ratio = Diffn / Diffe
        var GridAngle = (180 / Math.PI) * Math.atan(Ratio)

        if (Diffe > 0)
            GridBearing = 90 - GridAngle
        if (Diffe < 0)
            GridBearing = 270 - GridAngle
    }

    EndOfComputeBearing:

    // 'Compute convergence
    var Convergence = E_N_to_C(AtEast, AtNorth, a, b, E0, N0, f0, PHI0)

    // 'Compute (t-T) correction
    var t_minus_T = E_N_to_t_minus_T(AtEast, AtNorth, ToEast, ToNorth, a, b, E0, N0, f0, PHI0)

    // 'Compute initial azimuth
    var InitAzimuth = GridBearing + Convergence - t_minus_T

    // 'Set TrueAzimuth >=0 and <=360
    if (InitAzimuth < 0)
        return InitAzimuth + 360
    else if (InitAzimuth > 360)
        return InitAzimuth - 360
    else
        return InitAzimuth
}


var slieu = {
    a: 6378137.000,
    b: 6356752.314,
    DX: 191.9044143,
    DY: 39.30318279,
    DZ: 111.450,
    E0: 500000,
    f: 0.003352811,
    F0: 0.999900,
    LAM0: 107.75,
    N0: 0,
    PHI0: 0,
    Scale: 0.252906278,
    X_Rotation: -0.00928836,
    Y_Rotation: 0.0198,
    Z_Rotation: -0.00427372,
    Range: function (name) {
        return this[name]
    }
};


function NBT_to_VN2000_X(Latitude, longitude, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var VN2000X, VN2000Y, VN2000H
    // ' Lay cac thong so co ban cua VN2000
    var a = slieu.Range("a")
    var b = slieu.Range("b")
    var F = slieu.Range("f")
    if (b == 0) b = a * (1 - F);

    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // ' Lat long H to XYZ WGS84
    X = Lat_Long_H_to_X(Latitude, longitude, Height, a, b)
    Y = Lat_Long_H_to_Y(Latitude, longitude, Height, a, b)
    Z = Lat_H_to_Z(Latitude, Height, a, b)
    // 'XYZ to X1Y1Z1 VN2000
    X1 = Helmert_X(X, Y, Z, DX, Y_rot, Z_rot, S)
    Y1 = Helmert_Y(X, Y, Z, DY, X_rot, Z_rot, S)
    Z1 = Helmert_Z(X, Y, Z, DZ, X_rot, Y_rot, S)
    // 'X1Y1Z1 to Lat Long H VN2000
    B1 = XYZ_to_Lat(X1, Y1, Z1, a, b)
    L1 = XYZ_to_Long(X1, Y1)
    H1 = XYZ_to_H(X1, Y1, Z1, a, b)
    // 'Lat Long H to E N
    VN2000Y = Lat_Long_to_East(B1, L1, a, b, E0, muichieu, VTruc, KTtruc)
    VN2000X = Lat_Long_to_North(B1, L1, a, b, E0, N0, muichieu, VTruc, KTtruc)
    VN2000H = H1

    return VN2000X
}


function NBT_to_VN2000_Y(Latitude, longitude, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var VN2000X, VN2000Y, VN2000H
    // Lay cac thong so co ban cua VN2000
    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // Lat long H to XYZ WGS84
    X = Lat_Long_H_to_X(Latitude, longitude, Height, a, b)
    Y = Lat_Long_H_to_Y(Latitude, longitude, Height, a, b)
    Z = Lat_H_to_Z(Latitude, Height, a, b)
    //XYZ to X1Y1Z1 VN2000
    X1 = Helmert_X(X, Y, Z, DX, Y_rot, Z_rot, S)
    Y1 = Helmert_Y(X, Y, Z, DY, X_rot, Z_rot, S)
    Z1 = Helmert_Z(X, Y, Z, DZ, X_rot, Y_rot, S)
    //X1Y1Z1 to Lat Long H VN2000
    B1 = XYZ_to_Lat(X1, Y1, Z1, a, b)
    L1 = XYZ_to_Long(X1, Y1)
    H1 = XYZ_to_H(X1, Y1, Z1, a, b)
    //Lat Long H to E N
    VN2000Y = Lat_Long_to_East(B1, L1, a, b, E0, muichieu, VTruc, KTtruc)
    VN2000X = Lat_Long_to_North(B1, L1, a, b, E0, N0, muichieu, VTruc, KTtruc)
    VN2000H = H1

    return VN2000Y
}





function NBT_to_VN2000_Z(Latitude, longitude, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var VN2000X, VN2000Y, VN2000H
    // Lay cac thong so co ban cua VN2000

    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // Lat long H to XYZ WGS84
    X = Lat_Long_H_to_X(Latitude, longitude, Height, a, b)
    Y = Lat_Long_H_to_Y(Latitude, longitude, Height, a, b)
    Z = Lat_H_to_Z(Latitude, Height, a, b)
    //XYZ to X1Y1Z1 VN2000
    X1 = Helmert_X(X, Y, Z, DX, Y_rot, Z_rot, S)
    Y1 = Helmert_Y(X, Y, Z, DY, X_rot, Z_rot, S)
    Z1 = Helmert_Z(X, Y, Z, DZ, X_rot, Y_rot, S)
    //X1Y1Z1 to Lat Long H VN2000
    B1 = XYZ_to_Lat(X1, Y1, Z1, a, b)
    L1 = XYZ_to_Long(X1, Y1)
    H1 = XYZ_to_H(X1, Y1, Z1, a, b)
    //Lat Long H to E N
    VN2000Y = Lat_Long_to_East(B1, L1, a, b, E0, muichieu, VTruc, KTtruc)
    VN2000X = Lat_Long_to_North(B1, L1, a, b, E0, N0, muichieu, VTruc, KTtruc)
    VN2000H = H1

    return VN2000H
}



function NBT_to_WGS84_Lat(North, East, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var WGS84Lat, WGS84Long, WGS84H
    // Lay cac thong so co ban cua VN2000

    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // E N to Lat Long VN2000
    B1 = E_N_to_Lat(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    L1 = E_N_to_Long(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    H1 = Height
    // Lat long H to XYZ VN2000
    X1 = Lat_Long_H_to_X(B1, L1, Height, a, b)
    Y1 = Lat_Long_H_to_Y(B1, L1, Height, a, b)
    Z1 = Lat_H_to_Z(B1, L1, a, b)
    //X1Y1Z1 to XYZ WGS84
    X = Helmert_X(X1, Y1, Z1, DX, Y_rot, Z_rot, S)
    Y = Helmert_Y(X1, Y1, Z1, DY, X_rot, Z_rot, S)
    Z = Helmert_Z(X1, Y1, Z1, DX, X_rot, Y_rot, S)
    //XYZ to Lat Long H WGS84
    WGS84Lat = XYZ_to_Lat(X, Y, Z, a, b)
    WGS84Long = XYZ_to_Long(X, Y)
    WGS84H = XYZ_to_H(X, Y, Z, a, b)

    return WGS84Lat - 2.90417781181418E-03 - 0.00006/5

}

function NBT_to_WGS84_Long(North, East, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var WGS84Lat, WGS84Long, WGS84H
    // Lay cac thong so co ban cua VN2000

    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // E N to Lat Long VN2000
    B1 = E_N_to_Lat(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    L1 = E_N_to_Long(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    H1 = Height
    // Lat long H to XYZ VN2000
    X1 = Lat_Long_H_to_X(B1, L1, Height, a, b)
    Y1 = Lat_Long_H_to_Y(B1, L1, Height, a, b)
    Z1 = Lat_H_to_Z(B1, L1, a, b)
    //X1Y1Z1 to XYZ WGS84
    X = Helmert_X(X1, Y1, Z1, DX, Y_rot, Z_rot, S)
    Y = Helmert_Y(X1, Y1, Z1, DY, X_rot, Z_rot, S)
    Z = Helmert_Z(X1, Y1, Z1, DX, X_rot, Y_rot, S)
    //XYZ to Lat Long H WGS84
    WGS84Lat = XYZ_to_Lat(X, Y, Z, a, b)
    WGS84Long = XYZ_to_Long(X, Y)
    WGS84H = XYZ_to_H(X, Y, Z, a, b)

    return WGS84Long + 3.76088254250817E-03 - 0.00013*1.8
}

function NBT_to_WGS84_Long_Lat(North, East, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var WGS84Lat, WGS84Long, WGS84H
    // Lay cac thong so co ban cua VN2000

    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // E N to Lat Long VN2000
    B1 = E_N_to_Lat(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    L1 = E_N_to_Long(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    H1 = Height
    // Lat long H to XYZ VN2000
    X1 = Lat_Long_H_to_X(B1, L1, Height, a, b)
    Y1 = Lat_Long_H_to_Y(B1, L1, Height, a, b)
    Z1 = Lat_H_to_Z(B1, L1, a, b)
    //X1Y1Z1 to XYZ WGS84
    X = Helmert_X(X1, Y1, Z1, DX, Y_rot, Z_rot, S)
    Y = Helmert_Y(X1, Y1, Z1, DY, X_rot, Z_rot, S)
    Z = Helmert_Z(X1, Y1, Z1, DX, X_rot, Y_rot, S)
    //XYZ to Lat Long H WGS84
    WGS84Lat = XYZ_to_Lat(X, Y, Z, a, b)
    WGS84Long = XYZ_to_Long(X, Y)
    WGS84H = XYZ_to_H(X, Y, Z, a, b)
    return [WGS84Long + 3.76088254250817E-03 - 0.00013*1.8, WGS84Lat - 2.90417781181418E-03 - 0.00006/2]

}

function NBT_to_WGS84_H(North, East, Height) {
    var a, b, F
    var KTtruc, VTruc
    var muichieu
    var E0, N0
    var DX, DY, DZ
    var X_rot, Y_rot, Z_rot
    var S
    var X, Y, Z
    var X1, Y1, Z1
    var B1, L1, H1
    var WGS84Lat, WGS84Long, WGS84H
    // Lay cac thong so co ban cua VN2000
    a = slieu.Range("a")
    b = slieu.Range("b")
    F = slieu.Range("f")
    if (b == 0) b = a * (1 - F)
    KTtruc = slieu.Range("LAM0")
    VTruc = slieu.Range("PHI0")
    E0 = slieu.Range("E0")
    N0 = slieu.Range("N0")
    muichieu = slieu.Range("F0")
    DX = slieu.Range("DX")
    DY = slieu.Range("DY")
    DZ = slieu.Range("DZ")
    X_rot = slieu.Range("X_Rotation")
    Y_rot = slieu.Range("Y_Rotation")
    Z_rot = slieu.Range("Z_Rotation")
    S = slieu.Range("Scale")
    // E N to Lat Long VN2000
    B1 = E_N_to_Lat(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    L1 = E_N_to_Long(East, North, a, b, E0, N0, muichieu, VTruc, KTtruc)
    H1 = Height
    // Lat long H to XYZ VN2000
    X1 = Lat_Long_H_to_X(B1, L1, Height, a, b)
    Y1 = Lat_Long_H_to_Y(B1, L1, Height, a, b)
    Z1 = Lat_H_to_Z(B1, L1, a, b)
    //X1Y1Z1 to XYZ WGS84
    X = Helmert_X(X1, Y1, Z1, DX, Y_rot, Z_rot, S)
    Y = Helmert_Y(X1, Y1, Z1, DY, X_rot, Z_rot, S)
    Z = Helmert_Z(X1, Y1, Z1, DX, X_rot, Y_rot, S)
    //XYZ to Lat Long H WGS84
    WGS84Lat = XYZ_to_Lat(X, Y, Z, a, b)
    WGS84Long = XYZ_to_Long(X, Y)
    WGS84H = XYZ_to_H(X, Y, Z, a, b)


    return WGS84H - 91.4520812714472 - 4.2632564145606E-14 + 47

}

/*********************************************************************************************
 Nguyen Bach Thao
 Phong Ky thuat-Cong ty Dien luc Quang Nam
 Tel: 0983130081
 Email: thaoeqn1983@gmail.com
 *********************************************************************************************/
