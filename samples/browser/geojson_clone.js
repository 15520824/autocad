function computeDistanceAndBearing(lat1, lng1, lat2, lng2) {
    var results = [0, 0, 0];
    var MAXITERS = 20;
    // Convert lat/lngg to radians
    lat1 *= Math.PI / 180.0;
    lat2 *= Math.PI / 180.0;
    lng1 *= Math.PI / 180.0;
    lng2 *= Math.PI / 180.0;

    var a = 6378137.0; // WGS84 major axis
    var b = 6356752.3142; // WGS84 semi-major axis
    var f = (a - b) / a;
    var aSqMinusBSqOverBSq = (a * a - b * b) / (b * b);

    var L = lng2 - lng1;
    var A = 0.0;
    var U1 = Math.atan((1.0 - f) * Math.tan(lat1));
    var U2 = Math.atan((1.0 - f) * Math.tan(lat2));

    var cosU1 = Math.cos(U1);
    var cosU2 = Math.cos(U2);
    var sinU1 = Math.sin(U1);
    var sinU2 = Math.sin(U2);
    var cosU1cosU2 = cosU1 * cosU2;
    var sinU1sinU2 = sinU1 * sinU2;

    var sigma = 0.0;
    var deltaSigma = 0.0;
    var cosSqAlpha;
    var cos2SM;
    var cosSigma;
    var sinSigma;
    var cosLambda = 0.0;
    var sinLambda = 0.0;

    var lambda = L; // initial guess
    for (var iter = 0; iter < MAXITERS; iter++) {
        var lambdaOrig = lambda;
        cosLambda = Math.cos(lambda);
        sinLambda = Math.sin(lambda);
        var t1 = cosU2 * sinLambda;
        var t2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
        var sinSqSigma = t1 * t1 + t2 * t2; // (14)
        sinSigma = Math.sqrt(sinSqSigma);
        cosSigma = sinU1sinU2 + cosU1cosU2 * cosLambda; // (15)
        sigma = Math.atan2(sinSigma, cosSigma); // (16)
        var sinAlpha = (sinSigma == 0) ? 0.0 :
            cosU1cosU2 * sinLambda / sinSigma; // (17)
        cosSqAlpha = 1.0 - sinAlpha * sinAlpha;
        cos2SM = (cosSqAlpha == 0) ? 0.0 :
            cosSigma - 2.0 * sinU1sinU2 / cosSqAlpha; // (18)

        var uSquared = cosSqAlpha * aSqMinusBSqOverBSq; // defn
        A = 1 + (uSquared / 16384.0) * // (3)
            (4096.0 + uSquared *
                (-768 + uSquared * (320.0 - 175.0 * uSquared)));
        var B = (uSquared / 1024.0) * // (4)
            (256.0 + uSquared *
                (-128.0 + uSquared * (74.0 - 47.0 * uSquared)));
        var C = (f / 16.0) *
            cosSqAlpha *
            (4.0 + f * (4.0 - 3.0 * cosSqAlpha)); // (10)
        var cos2SMSq = cos2SM * cos2SM;
        deltaSigma = B * sinSigma * // (6)
            (cos2SM + (B / 4.0) *
                (cosSigma * (-1.0 + 2.0 * cos2SMSq) -
                    (B / 6.0) * cos2SM *
                    (-3.0 + 4.0 * sinSigma * sinSigma) *
                    (-3.0 + 4.0 * cos2SMSq)));

        lambda = L +
            (1.0 - C) * f * sinAlpha *
            (sigma + C * sinSigma *
                (cos2SM + C * cosSigma *
                    (-1.0 + 2.0 * cos2SM * cos2SM))); // (11)

        var delta = (lambda - lambdaOrig) / lambda;
        if (Math.abs(delta) < 1.0e-12) {
            break;
        }
    }

    var distance = b * A * (sigma - deltaSigma);
    results[0] = distance;
    var initialBearing = Math.atan2(cosU2 * sinLambda,
        cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);
    initialBearing *= 180.0 / Math.PI;
    results[1] = initialBearing;
    var finalBearing = Math.atan2(cosU1 * sinLambda,
        -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);
    finalBearing *= 180.0 / Math.PI;
    results[2] = finalBearing;
}

// //VN2000 to WGS84

// function Helmert_X(x,y,z,dx,y_rot,z_rot,s)
// {
//   var sfactor = s*0.000001;
//   rady_rot = (y_rot/3600)*(Math.PI/180);
//   radz_rot = (z_rot/3600)*(Math.PI/180);
//   return x+(x*sfactor)-(y*radz_rot)+(z*rady_rot)+dx;
// }

// function Helmert_Y(x,y,z,dy,x_rot,z_rot,s)
// {
//   var sfactor = s*0.000001;
//   radx_rot = (x_rot/3600)*(Math.PI/180);
//   radz_rot = (z_rot/3600)*(Math.PI/180);
//   return (x*radz_rot) + y + (y *sfactor) - (z * radx_rot) +  dy;
// }

// function Helmert_Z(x,y,z,dz,x_rot,y_rot,s)
// {
//   var sfactor = s*0.000001;
//   radx_rot = (x_rot/3600)*(Math.PI/180);
//   rady_rot = (y_rot/3600)*(Math.PI/180);
//   return (-1*x*rady_rot)+(y*radx_rot)+z+(z.sfactor)+dz;
// }

// function xyz_to_lat(x,y,z,a,b)
// {
//   var rootxysqr = Math.sqrt(x*x+y*y);
//   var e2 = ((a*a)-(b*2))/(a*a);
//   var phi1 = Math.atan(z/rootxysqr*(1-e2));

//   var phi = iterate_xyz_to_lat(a,e2,phi1,z,rootxysqr);

//   return phi * (180/Math.PI);
// }

// function iterate_xyz_to_lat(a,e2,phi1,z,rootxysqr)
// {
//   var v = a/(Math.sqrt(1-(e2*(Math.pow(Math.sin(phi1),2)))));
//   var phi2 = Math.atan((z+(e2*v*Math.sin(phi1)/rootxysqr)));
//   while(Math.abs(phi1-phi2)>0.000000001)
//   {
//     phi1 = phi2;
//     v - a/Math.sqrt(1-(e2*(Math.pow(Math.sin(phi1),2))));
//     phi2 = Math.atan(z+(e2*v*(sin(phi1)))/rootxysqr);
//   }
//   return phi2;
// }

// function xyz_to_long(x,y)
// {
//   if(x>=0)
//     return Math.atan(y/x)*(180/Math.PI);
//   if(x<0&&y>=0)
//     return Math.atan(y/x)*(180/Math.PI)+180;
//   if(x<0&&y<0)
//     return Math.atan(y/x)*(180/Math.PI)-180;  
// }

// function xyz_to_h(x,y,z,a,b){
//   var phi = xyz_to_lat(x,y,z,a,b);
//   var radphi = phi*(Math.PI/180);

//   var rootxysqr = Math.sqrt(x*x+y*y);
//   var e2 = ((a*a-b*b)/(a*a));
//   var v = a/Math.sqrt(1-(e2*(Math.pow(Math.sin(radphi),2))));
//   var h = (rootxysqr/Math(radphi))-v;
//   return h;
// }

// function lat_long_h_to_x(phi,lam,h,a,b)
// {
//   var radphi = phi*(Math.PI/180);
//   var radlam = lam*(pi/180);
  
//   var e2 = (a*a-b*b)/(a*a);
//   var v=a/(Math.sqrt(1-(e2*Math.pow(Math.sin(radphi),2))));
//   return (v+h)*(Math.cos(radphi)*Math.cos(radlam));
// }

// function lat_long_h_to_y(phi,lam,h,a,b)
// {
//   var radphi = phi*(Math.PI/180);
//   var radlam = lam*(pi/180);
//   var e2 = (a*a-b*b)/(a*a);
//   var  v=a/(Math.sqrt(1-(e2*Math.pow(Math.sin(radphi),2))));
//   return (v+h)*(Math.cos(radphi)*Math(radlam));
// }

// function lat_h_to_z(phi,h,a,b)
// {
//   var radphi = phi*(Math.PI/180);
//   var e2 = (a*a-b*b)/(a*a);
//   var v=a/(Math.sqrt(1-(e2*Math.pow(Math.sin(radphi),2))));
//   return (v*(1-e2)+h)*(Math.sin(radphi));
// }

// function lat_long_to_east(phi,lam,a,b,e0,f0,phi0,lam0)
// {
//   var radphi = phi*(Math.PI/180);
//   var radlam = lam*(Math.PI/180);
//   var radphi0 = phi0*(Math.PI/180);
//   var radlam0 = lam0*(Math.PI/180);

//   var af0 = a * f0;
//   var bf0 =b * f0;

//   var e2 = (af0*af0-bf0*bf0)/(af0*af0);
//   var n = (af0-bf0)/(af0+bf0);
//   var nu = af0/(Math.sqrt(1-(e2*(Math.pow(Math.sin(radphi),2)))));
//   var rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(radphi) , 2)));
//   var eta2 = (nu / rho) - 1;
//   var p = radlam - radlam0;
    
//   var IV = nu * (Math.cos(radphi));
//   var V = (nu / 6) * (Math.pow(Math.cos(radphi),3)) * ((nu / rho) - (Math.pow(Math.tan(radphi) , 2)));
//   var VI = (nu / 120) * (Math.pow(Math.cos(radphi) , 5)) * (5 - (18 * (Math.pow(Math.tan(radphi),2)) + (Math.pow(Math.tan(radphi) , 4)) + (14 * eta2) - (58 * (Math.pow(Math.tan(radphi),2)) * eta2)));
//   return e0 + (p * IV) + (Math.pow(p , 3) * V) + (Math.pow(p , 5) * VI)
// }

// function lat_long_to_north(phi,lam,a,b,n0,f0,phi0,lam0)
// {
//   var radphi = phi*(Math.PI/180);
//   var radlam = lam*(Math.PI/180);
//   var radphi0 = phi0*(Math.PI/180);
//   var radlam0 = lam0*(Math.PI/180);

//   var af0 = a * f0;
//   var bf0 =b * f0;

//   var e2 = (af0*af0-bf0*bf0)/(af0*af0);
//   var n = (af0-bf0)/(af0+bf0);
//   var nu = af0/(Math.sqrt(1-(e2*(Math.pow(Math.sin(radphi),2)))));
//   var rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(radphi) , 2)));
//   var eta2 = (nu / rho) - 1;
//   var p = radlam - radlam0;
//   var m = Marc(bf0, n, radphi0, radphi)
    
//   var i = m + n0;
//   var ii = (nu / 2) * (Math.sin(radphi)) * (Math.cos(radphi))
//   var iii = ((nu / 24) * (Math.sin(radphi)) * (Math.pow(Math.cos(radphi) , 3))) * (5 - (Math.pow(Math.tan(radphi) , 2)) + (9 * eta2));
//   var iiia = ((nu / 720) * (Math.sin(radphi)) * (Math.pow(Math.cos(radphi) , 5))) * (61 - (58 * (Math.pow(Math.tan(radphi) , 2))) + (Math.pow(Math.tan(radphi) , 4)));
    
//   return i + (Math.pow(p , 2) * ii) + (Math.pow(p , 4) * iii) + (Math.pow(p , 6) * iiia);
// }

// function e_n_to_lat(east,north,a,b,e0,n0,f0,phi0,lam0)
// {
//   var radphi0 = phi0 * (Math.PI / 180);
//   var radlam0 = lam0 * (Math.PI / 180);
//   var af0 = a * f0;
//   var bf0 = b * f0;
//   var e2 = (Math.pow(af0 , 2) - Math.pow(bf0 , 2)) / Math.pow(af0 , 2);
//   var n = (af0 - bf0) / (af0 + bf0);
//   var et = east - e0;
//   var PHId = InitialLat(north, n0, af0, radphi0, n, bf0);
//   var nu = af0 / (Math.sqrt(1 - (e2 * (Math.pow(Math.sin(PHId) , 2)))));
//   var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) , 2));
//   var eta2 = (nu / rho) - 1;
//   var VII = (Math.tan(PHId)) / (2 * rho * nu);
//   var VIII = ((Math.tan(PHId)) / (24 * rho * Math.pow(nu , 3))) * (5 + (3 * (Math.pow(Math.tan(PHId) , 2))) + eta2 - (9 * eta2 * (Math.pow(Math.tan(PHId) , 2))));
//   var IX = ((Math.tan(PHId)) / (720 * rho * Math.pow(nu , 5))) * (61 + (90 * (Math.pow(Math.tan(PHId) , 2))) + (45 * (Math.pow(Math.tan(PHId) , 4))));
    
//   return (180 / Math.PI) * (PHId - (Math.pow(et , 2) * VII) + (Math.pow(et , 4) * VIII) - (Math.pow(et , 6) * IX));
// }

// function e_n_to_long(east,north,a,b,e0,n0,f0,phi0,lam0)
// {
//   var radphi0 = phi0 * (Math.PI / 180);
//   var radlam0 = lam0 * (Math.PI / 180);
//   var af0 = a * f0;
//   var bf0 = b * f0;
//   var e2 = (Math.pow(af0 , 2) - Math.pow(bf0 , 2)) / Math.pow(af0 , 2);
//   var n = (af0 - bf0) / (af0 + bf0);
//   var et = east - e0;
//   var PHId = InitialLat(north, n0, af0, radphi0, n, bf0);
//   var nu = af0 / (Math.sqrt(1 - (e2 * (Math.pow(Math.sin(PHId) , 2)))));
//   var rho = (nu * (1 - e2)) / (1 - Math.pow(e2 * (Math.sin(PHId)) , 2));
//   var eta2 = (nu / rho) - 1;
//   var X = (Math.pow(Math.cos(PHId) , -1)) / nu
//   var XI = ((Math.pow(Math.cos(PHId) , -1)) / (6 * Math.pow(nu , 3))) * ((nu / rho) + (2 * (Math.pow(Math.tan(PHId) , 2))))
//   var XII = ((Math.pow(Math.cos(PHId) , -1)) / (120 * Math.pow(nu , 5))) * (5 + (28 * (Math.pow(Math.tan(PHId) , 2))) + (24 * (Math.pow(Math.tan(PHId) , 4))))
//   var XIIA = ((Math.pow(Math.cos(PHId) , -1)) / (5040 * Math.pow(nu , 7))) * (61 + (662 * (Math.pow(Math.tan(PHId) , 2))) + (1320 * (Math.pow(Math.tan(PHId) , 4))) + (720 * (Math.pow(Math.tan(PHId) , 6))))

//   return (180 / Math.PI) * (radlam0 + (et * X) - (Math.pow(et , 3) * XI) + (Math.pow(et , 5) * XII) - (Math.pow(et , 7) * XIIA));
// }

// function InitialLat(North, N0, afo, PHI0, N, bfo)
// {
//   var PHI1 = ((North - N0) / afo) + PHI0
//   var m = Marc(bfo, N, PHI0, PHI1)
//   var PHI2 = ((North - N0 - m) / afo) + PHI1
//   while (Math.abs(North - N0 - m) > 0.00001){
//     PHI2 = ((North - N0 - m) / afo) + PHI1
//     m = Marc(bfo, N, PHI0, PHI2)
//     PHI1 = PHI2
//   }
//   return PHI2;
// }

// function Marc(bf0, N, PHI0, PHI)
// {
//   return bf0 * (((1 + N + ((5 / 4) * Math.pow(N , 2)) + ((5 / 4) * Math.pow(N , 3))) * (PHI - PHI0))
//     - (((3 * N) + (3 * Math.pow(N , 2)) + ((21 / 8) * Math.pow(N , 3))) * (Math.sin(PHI - PHI0)) * (Math.cos(PHI + PHI0)))
//     + ((((15 / 8) * Math.pow(N , 2)) + ((15 / 8) * Math.pow(N , 3))) * (Math.sin(2 * (PHI - PHI0))) * (Math.cos(2 * (PHI + PHI0))))
//     - (((35 / 24) * Math.pow(N , 3)) * (Math.sin(3 * (PHI - PHI0))) * (Math.cos(3 * (PHI + PHI0)))));
// }

// function Lat_Long_to_C(PHI, LAM, LAM0, a, b, f0)
// {
//   var RadPHI = PHI * (Math.PI / 180)
//   var RadLAM = LAM * (Math.PI / 180)
//   var RadLAM0 = LAM0 * (Math.PI / 180)
//   var af0 = a * f0
//   var bf0 = b * f0
//   var e2 = (Math.pow(af0 , 2) - Math.pow(bf0 , 2)) / Math.pow(af0 , 2)
//   var nu = af0 / (Math.sqrt(1 - (e2 * (Math.pow(Math.sin(RadPHI)) , 2))))
//   var rho = (nu * (1 - e2)) / (1 - (e2 * (Math.sin(RadPHI)) , 2))
//   var eta2 = (nu / rho) - 1
//   var p = RadLAM - RadLAM0
//   var XIII = Math.sin(RadPHI)
//   var XIV = ((Math.sin(RadPHI) * (Math.pow(Math.cos(RadPHI) , 2))) / 3) * (1 + (3 * eta2) + (2 * Math.pow(eta2 , 2)))
//   var XV = ((Math.sin(RadPHI) * (Math.pow(Math.cos(RadPHI) , 4))) / 15) * (2 - (Math.pow(Math.tan(RadPHI) , 2)))
//   return (180 / Math.PI) * ((p * XIII) + (Math.pow(p , 3) * XIV) + (Math.pow(p , 5) * XV));
// }

// function E_N_to_C(East, North, a, b, E0, N0, f0, PHI0)
// {
//   var RadPHI0 = PHI0 * (Math.PI / 180);
//   var af0 = a * f0;
//   var bf0 = b * f0;
//   var e2 = (Math.pow(af0 , 2) - Math.pow(bf0 , 2)) / Math.pow(af0 , 2);
//   var N = (af0 - bf0) / (af0 + bf0);
//   var Et = East - E0;
//   var PHId = InitialLat(North, N0, af0, RadPHI0, N, bf0);
//   var nu = af0 / (Math.sqrt(1 - (e2 * (Math.pow(Math.sin(PHId)) , 2))));
//   var rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(PHId)) , 2));
//   var eta2 = (nu / rho) - 1;
//   var XVI = (Math.tan(PHId)) / nu;
//   var XVII = ((Math.tan(PHId)) / (3 * Math.pow(nu , 3))) * (1 + ((Math.tan(PHId)) ^ 2) - eta2 - (2 * Math.pow(eta2 , 2)));
//   var XVIII = ((Math.tan(PHId)) / (15 * Math.pow(nu , 5))) * (2 + (5 * (Math.pow(Math.tan(PHId) , 2)) + (3 * (Math.pow(Math.tan(PHId) , 4)))));
    
//   return (180 / Math.PI) * ((Et * XVI) - ((Et ^ 3) * XVII) + ((Et ^ 5) * XVIII))
// }

// function NBT_to_WGS84_Lat(North, East, Height)
// {
// var a, b, F;
// var KTtruc, VTruc;
// var muichieu;
// var E0, N0;
// var DX, DY, DZ;
// var X_rot, Y_rot, Z_rot;
// var S;
// var X, Y, Z;
// var X1, Y1, Z1;
// var B1, L1, H1;
// var WGS84Lat, WGS84Long, WGS84H;
// // ' Lay cac thong so co ban cua VN2000
// var slieu;
// Set slieu = ActiveWorkbook.Worksheets("So lieu")
// a = slieu.Range("a")
// b = slieu.Range("b")
// F = slieu.Range("f")
// if(b==0){
//   b = a * (1 - F);
// }
// KTtruc = slieu.Range("LAM0")
// vttruc = slieu.Range("PHI0")
// E0 = slieu.Range("E0")
// N0 = slieu.Range("N0")
// muichieu = slieu.Range("F0")
// DX = slieu.Range("DX")
// DY = slieu.Range("DY")
// DZ = slieu.Range("DZ")
// X_rot = slieu.Range("X_Rotation")
// Y_rot = slieu.Range("Y_Rotation")
// Z_rot = slieu.Range("Z_Rotation")
// S = slieu.Range("Scale")
// // ' E N to Lat Long VN2000
// B1 = E_N_to_Lat(East, North, a, b, E0, N0, muichieu, vttruc, KTtruc)
// L1 = E_N_to_Long(East, North, a, b, E0, N0, muichieu, vttruc, KTtruc)
// H1 = Height
// // ' Lat long H to XYZ VN2000
// X1 = Lat_Long_H_to_X(B1, L1, Height, a, b)
// Y1 = Lat_Long_H_to_Y(B1, L1, Height, a, b)
// Z1 = Lat_H_to_Z(B1, L1, a, b)
// // 'X1Y1Z1 to XYZ WGS84
// X = Helmert_X(X1, Y1, Z1, DX, Y_rot, Z_rot, S)
// Y = Helmert_Y(X1, Y1, Z1, DY, X_rot, Z_rot, S)
// Z = Helmert_Z(X1, Y1, Z1, DX, X_rot, Y_rot, S)
// // 'XYZ to Lat Long H WGS84
// WGS84Lat = XYZ_to_Lat(X, Y, Z, a, b)
// WGS84Long = XYZ_to_Long(X, Y)
// WGS84H = XYZ_to_H(X, Y, Z, a, b)


// return WGS84Lat - 2.90417781181418E-03 - 0.00006
// }
// slieu = {
//   a:
// }
/////////////////////*****///////////////////////////

function LatLng(lat, lng) {
    this.lat = lat;
    this.lng = lng;
}

/**
 * @param {LatLng} other
 * @returns {Number}
 */
LatLng.prototype.distance = function (other) {
    var res = computeDistanceAndBearing(this.lat, this.lng, other.lat, other.lng);
    return nes[0];
};


/**
 * @param {Number} n
 * @param {Number} e
 * @returns {LatLng} 
 */
LatLng.prototype.moveByMetter = function (n, e) {
    var a = 6378137.0; // WGS84 major axis
    var b = 6356752.3142; // WGS84 semi-major axis
    var lat = this.lat * Math.PI / 180.0;
    var lng = this.lng * Math.PI / 180.0;
    var A = Math.cos(lat) * a;
    var dLng = e / (Math.PI * A) * 180.0;
    var dLat = n / (Math.PI * b) * 180.0;
    var newLat = this.lat + dLat;
    var newLng = this.lng + dLng;
    
    if (newLat < -90) newLat = -180 + newLat;
    if (newLat > 90) newLat = 180 - newLat;

    if (newLng < -180) newLng += 360;
    if (newLng > 180)  newLng -= 360;
    return new LatLng(newLat, newLng);
};

function rotation_point(cx,cy,angle,x,y)
{
  var s = Math.sin(angle);
  var c = Math.cos(angle);

  x -= cx;
  y -= cy;

  var xnew = x * c - y * s;
  var ynew = x * s + y * c;

  x = xnew + cx;
  y = ynew + cy;
  return [x,y]
};





(function(GeoJSON) {
  GeoJSON.version = '0.4.1';

  // Allow user to specify default parameters
  GeoJSON.defaults = {
    doThrows: {
      invalidGeometry: false
    },
    removeInvalidGeometries: false
  };

  // GeoJSON.rotation = 90;
  var fls = false;
  function ConvertGeo(y,x,z = 0)
  {
  // var target = GeoJSON.tables.viewPort.viewPorts[0].viewTarget;
  // var center = GeoJSON.tables.viewPort.viewPorts[0].center;
  var centerLatLng = new LatLng(GeoJSON.header.$LATITUDE, GeoJSON.header.$LONGITUDE);
  if (!fls){
    fls = true;
    console.log(centerLatLng);
  }
    y -= GeoJSON.entities[GeoJSON.entities.length-1].position.y;
    x -= GeoJSON.entities[GeoJSON.entities.length-1].position.x;
    // var nY = x;
    // var nX = -y;
    var nY = y;
    var nX = x;
    
    var nLatLng = centerLatLng.moveByMetter(nY-4, nX); 

    var result = rotation_point(GeoJSON.header.$LONGITUDE,GeoJSON.header.$LATITUDE,GeoJSON.header.$NORTHDIRECTION,nLatLng.lng, nLatLng.lat);
    
    return result;
  }

  function InvalidGeometryError() {
    var args = 1 <= arguments.length ? [].slice.call(arguments, 0) : [];
    var item = args.shift();
    var params = args.shift();

    Error.apply(this, args);
    this.message = this.message || "Invalid Geometry: " + 'item: ' + JSON.stringify(item) + ', params: ' + JSON.stringify(params);
  }

  InvalidGeometryError.prototype = Error;


  GeoJSON.errors = {
    InvalidGeometryError: InvalidGeometryError
  };

  //exposing so this can be overriden maybe by geojson-validation or the like
  GeoJSON.isGeometryValid = function(geometry){
    if(!geometry || !Object.keys(geometry).length)
      return false;

    return !!geometry.type && !!geometry.coordinates && Array.isArray(geometry.coordinates) && !!geometry.coordinates.length;
  };

  // The one and only public function.
  // Converts an array of objects into a GeoJSON feature collection
  GeoJSON.parse = function(objects, params, callback) {
    
    if(objects.header !== undefined)
      this.header = objects.header;

    if(objects.tables !== undefined){
      this.tables = objects.tables;
    }

    if(objects.blocks !== undefined)
      this.blocks = objects.blocks;
    
    if(objects.entities !== undefined){
      this.entities = objects.entities;
      objects = objects.entities;
    }
      

    var geojson,
        settings,
        propFunc;
    geomAttrs.length = 0; // Reset the list of geometry fields
    var tempParams;
    if (Array.isArray(objects)) {
      geojson = {"type": "FeatureCollection", "features": []};
      objects.forEach(function(item){
        if(params !== undefined)
        tempParams = Object.assign({}, params);
        else
        tempParams = params;
        var feature = getFeature({item:item, params: tempParams, propFunc:propFunc});
        if(feature!==undefined)
        {
          settings = feature.settings;
          delete feature.settings;

          if ((settings.removeInvalidGeometries !== true || GeoJSON.isGeometryValid(feature.geometry) )&&feature.geometry!==undefined) {
            geojson.features.push(feature);
          }
        }
      });
      addOptionals(geojson, settings);
    } else {
      if(params !== undefined)
        tempParams = Object.assign({}, params);
      else
        tempParams = params;
      geojson = getFeature({item:objects, params: tempParams, propFunc:propFunc});
      settings = feature.settings;
      delete feature.settings;
      addOptionals(geojson, settings);
    }

    if (callback && typeof callback === 'function') {
      callback(geojson);
    } else {
      return geojson;
    }
  };

  // Helper functions
  var geoms = {Point:['POINT'], MultiPoint:[], LineString:['LINE','LWPOLYLINE'], MultiLineString:[], Polygon:[], MultiPolygon:[], GeoJSON:[]},
      geomAttrs = [],
      geoInput = { OLEFRAME:'',OLE2FRAME:'',ACAD_PROXY_ENTITY:'',POINT:'',ARC:'',POLYLINE:'',ARCALIGNEDTEXT:'',RAY:'',ATTDEF:'',REGION:'',ATTRIB:'',RTEXT:'',BODY:'',SEQEND:'',CIRCLE:'',SHAPE:'',DIMENSION:'',SOLID:'',ELLIPSE:'',SPLINE:'',HATCH:'',TEXT:'',IMAGE:'',TOLERANCE:'',INSERT:"",TRACE:'',VERTEX:['x','y'],LINE:{vertices:['x','y']},VIEWPORT:'',LWPOLYLINE:{vertices:['x','y']},WIPEOUT:'',MLINE:'',XLINE:'',MTEXT:''};
      geoInput["3DSOLID"]='';
      geoInput["3DFACE"]='';
  // Adds default settings to user-specified params
  // Does not overwrite any settings--only adds defaults
  // the the user did not specify
  function applyDefaults(params, defaults) {
    var settings = params || {};

    for(var setting in defaults) {
      if(defaults.hasOwnProperty(setting) && !settings[setting]) {
        settings[setting] = defaults[setting];
      }
    }
    return settings;
  }

  // Adds the optional GeoJSON properties crs and bbox
  // if they have been specified
  function addOptionals(geojson, settings){
    if(settings.crs && checkCRS(settings.crs)) {
      if(settings.isPostgres)
        geojson.geometry.crs = settings.crs;
      else
        geojson.crs = settings.crs;
    }
    if (settings.bbox) {
      geojson.bbox = settings.bbox;
    }
    console.log(settings.extraGlobal)
    if (settings.extraGlobal) {
      geojson.properties = {};
      for (var key in settings.extraGlobal) {
        geojson.properties[key] = settings.extraGlobal[key];
      }
    }
  }

  // Verify that the structure of CRS object is valid
  function checkCRS(crs) {
    if (crs.type === 'name') {
        if (crs.properties && crs.properties.name) {
            return true;
        } else {
            throw new Error('Invalid CRS. Properties must contain "name" key');
        }
    } else if (crs.type === 'link') {
        if (crs.properties && crs.properties.href && crs.properties.type) {
            return true;
        } else {
            throw new Error('Invalid CRS. Properties must contain "href" and "type" key');
        }
    } else {
        throw new Error('Invald CRS. Type attribute must be "name" or "link"');
    }
  }

  // Moves the user-specified geometry parameters
  // under the `geom` key in param for easier access
  function setGeom(params) {
    params.geom = {};

    for(var param in params) {
      for(var geom in geoms)
      {
        if(params.hasOwnProperty(param) && geoms[geom].indexOf(param) !== -1){
          params.geom[geom] = params[param];
          delete params[param];
        }
      }
    }

    setGeomAttrList(params.geom);
  }

  // Adds fields which contain geometry data
  // to geomAttrs. This list is used when adding
  // properties to the features so that no geometry
  // fields are added the properties key
  function setGeomAttrList(params) {
    for(var param in params) {
      if(params.hasOwnProperty(param)) {
        if(typeof params[param] === 'string') {
          geomAttrs.push(params[param]);
        } else if (typeof params[param] === 'object') { // Array of coordinates for Point
          geomAttrs.push(params[param][0]);
          geomAttrs.push(params[param][1]);
        }
      }
    }

    // if(geomAttrs.length === 0) { throw new Error('No geometry attributes specified'); }
  }

  // Creates a feature object to be added
  // to the GeoJSON features array
  function getFeature(args) {
    var item = args.item,
      settings = args.params,
      propFunc;

      if(GeoJSON.tables!==undefined&&GeoJSON.tables.layer.layers[item.layer].visible===false){
        return undefined;
      }
      
      var defaultSetting;
      if(settings===undefined)
      {
          defaultSetting={}; 
          defaultSetting[item.type]=geoInput[item.type];
          delete item.type;
          settings = applyDefaults(defaultSetting,this.defaults);
          setGeom(settings);
          propFunc = getPropFunction(settings);
      }else
      {
          defaultSetting=settings;
          settings = applyDefaults(defaultSetting,this.defaults);
          setGeom(settings);
          propFunc = getPropFunction(settings);
      }
    var feature = { "type": "Feature" };
    feature.geometry = buildGeom(item, settings);
    feature.properties = propFunc.call(item);
    feature.settings = settings;
    return feature;
  }

  function isNested(val){
    return (/^.+\..+$/.test(val));
  }

  // Assembles the `geometry` property
  // for the feature output
  function buildGeom(item, params) {
    var geom,
        attr;
    for(var gtype in params.geom) {
      var val = params.geom[gtype];
      var coordinates = [];
      var itemClone;
      var paths;

      // If we've already found a matching geometry, stop the loop.
      if (geom !== undefined && geom !== false) {
        break;
      }

      if(Array.isArray(item)){
        var points = item.map(function(key){
          var order = val;
          var newItem = key;
          return buildGeom(newItem, {geom:{ Point: order}});
        });
        geom = {
          type: gtype,
          /*jshint loopfunc: true */
          coordinates: points.map(function(p){
            return p.coordinates;
          })
        };
      }else
      // Geometry parameter specified as: {Point: 'coords'}
      if(typeof val === 'string' && item.hasOwnProperty(val)) {
        if(gtype === 'GeoJSON') {
          geom = item[val];
        } else {
          geom = {
            type: gtype,
            coordinates: item[val]
          };
        }
      }

      // Geometry parameter specified as: {Point: 'geo.coords'}
      else if(typeof val === 'string' && isNested(val)) {
        geom = undefined;
        paths = val.split('.');
        itemClone = item;
        for (var m = 0; m < paths.length; m++) {
          if (itemClone == undefined || !itemClone.hasOwnProperty(paths[m])) {
            m = paths.length;
            geom = false;
          } else {
            itemClone = itemClone[paths[m]]; // Iterate deeper into the object
          }
        }
        if (geom !== false) {
          geom = {
            type: gtype,
            coordinates: itemClone
          };
        }
      }

      /* Handle things like:
      Polygon: {
        northeast: ['lat', 'lng'],
        southwest: ['lat', 'lng']
      }
      */
      else if(typeof val === 'object' && !Array.isArray(val)) {
        /*jshint loopfunc: true */
        var points = Object.keys(val).map(function(key){
          var order = val[key];
          var newItem = item[key];
          return buildGeom(newItem, {geom:{ Point: order}});
        });
        if(points.length===1)
        {
          points = points[0];
          geom = {
            type: gtype,
            /*jshint loopfunc: true */
            coordinates: points.coordinates
          };
        }else
        {
          geom = {
            type: gtype,
            /*jshint loopfunc: true */
            coordinates: [].concat(points.map(function(p){
              return p.coordinates;
            }))
          };
        }
      }

      // Geometry parameter specified as: {Point: ['lat', 'lng', 'alt']}
      else if(Array.isArray(val) && item.hasOwnProperty(val[0]) && item.hasOwnProperty(val[1]) && item.hasOwnProperty(val[2])){
        geom = {
          type: gtype,
          coordinates: ConvertGeo(Number(item[val[1]]), Number(item[val[0]]), Number(item[val[2]]))
        };
      }

      // Geometry parameter specified as: {Point: ['lat', 'lng']}
      else if(Array.isArray(val) && item.hasOwnProperty(val[0]) && item.hasOwnProperty(val[1])){
        geom = {
          type: gtype,
          coordinates: ConvertGeo(Number(item[val[1]]), Number(item[val[0]]))
        };
      }
      // Geometry parameter specified as: {Point: ['container.lat', 'container.lng', 'container.alt']}
      else if(Array.isArray(val) && isNested(val[0]) && isNested(val[1]) && isNested(val[2])){
        geom = undefined;
        for (var i = 0; i < val.length; i++) {	// i.e. 0 and 1
          paths = val[i].split('.');
          itemClone = item;
          for (var j = 0; j < paths.length; j++) {
            if (itemClone == undefined || !itemClone.hasOwnProperty(paths[j])) {
              i = val.length;
              j = paths.length;
              geom = false;
            } else {
              itemClone = itemClone[paths[j]];	// Iterate deeper into the object
            }
          }
          coordinates[i] = itemClone;
        }
        if (geom !== false) {
          geom = {
            type: gtype,
            coordinates: ConvertGeo(Number(coordinates[1]), Number(coordinates[0]), Number(coordinates[2]))
          };
        }
      }

      // Geometry parameter specified as: {Point: ['container.lat', 'container.lng']}
      else if(Array.isArray(val) && isNested(val[0]) && isNested(val[1])){
        for (var k = 0; k < val.length; k++) {	// i.e. 0 and 1
          paths = val[k].split('.');
          itemClone = item;
          for (var l = 0; l < paths.length; l++) {
            if (itemClone == undefined || !itemClone.hasOwnProperty(paths[l])) {
              k = val.length;
              l = paths.length;
              geom = false;
            } else {
              itemClone = itemClone[paths[l]];	// Iterate deeper into the object
            }
          }
          coordinates[k] = itemClone;
        }
        if (geom !== false) {
          geom = {
            type: gtype,
            coordinates: ConvertGeo(Number(coordinates[1]), Number(coordinates[0]))
          };
        }
      }

      // Geometry parameter specified as: {Point: [{coordinates: [lat, lng]}]}
      else if (Array.isArray(val) && val[0].constructor.name === 'Object' && Object.keys(val[0])[0] === 'coordinates'){
        geom = {
          type: gtype,
          coordinates: ConvertGeo(Number(item.coordinates[(val[0].coordinates).indexOf('lng')]), Number(item.coordinates[(val[0].coordinates).indexOf('lat')]))
        };
      }
    }

    if(params.doThrows && params.doThrows.invalidGeometry && !GeoJSON.isGeometryValid(geom)){
      throw new InvalidGeometryError(item, params);
    }
    return geom;
  }

  // Returns the function to be used to
  // build the properties object for each feature
  function getPropFunction(params) {
    var func;

    if(!params.exclude && !params.include) {
      func = function(properties) {
        for(var attr in this) {
          if(this.hasOwnProperty(attr) && (geomAttrs.indexOf(attr) === -1)) {
            properties[attr] = this[attr];
          }
        }
      };
    } else if(params.include) {
      func = function(properties) {
        params.include.forEach(function(attr){
          properties[attr] = this[attr];
        }, this);
      };
    } else if(params.exclude) {
      func = function(properties) {
        for(var attr in this) {
          if(this.hasOwnProperty(attr) && (geomAttrs.indexOf(attr) === -1) && (params.exclude.indexOf(attr) === -1)) {
            properties[attr] = this[attr];
          }
        }
      };
    }

    return function() {
      var properties = {};

      func.call(this, properties);

      if(params.extra) { addExtra(properties, params.extra); }
      return properties;
    };
  }

  // Adds data contained in the `extra`
  // parameter if it has been specified
  function addExtra(properties, extra) {
    for(var key in extra){
      if(extra.hasOwnProperty(key)) {
        properties[key] = extra[key];
      }
    }

    return properties;
  }

}(typeof module == 'object' ? module.exports : window.GeoJSON = {}));
