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
    x -=  GeoJSON.entities[GeoJSON.entities.length-1].position.x;
    // var newPos - 
    // var nY = x;
    // var nX = -y;
    var nX = x;
    var nY = -y;
    var nX1 = -nY;
    var nY1 = nX;

    nX1 += GeoJSON.header.VN2000_X;
    nY1 += GeoJSON.header.VN2000_Y;
    // var nLatLng = centerLatLng.moveByMetter(nY-4, nX); 
    var result = rotation_point(GeoJSON.header.VN2000_X,GeoJSON.header.VN2000_Y,-GeoJSON.header.$NORTHDIRECTION,nX1, nY1);
    nX1 = result[0];
    nY1 = result[1];
    // result = [NBT_to_WGS84_Long(nX1, nY1, 0),NBT_to_WGS84_Lat(nX1, nY1, 0)];
    return NBT_to_WGS84_Long_Lat(nX1, nY1, 0);
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
    if(objects.entities !== undefined)
      this.entities = objects.entities;
    if(objects.header !== undefined){
      this.header = objects.header;
      this.header.VN2000_X = NBT_to_VN2000_X(this.header.$LATITUDE,this.header.$LONGITUDE,0);
      this.header.VN2000_Y = NBT_to_VN2000_Y(this.header.$LATITUDE,this.header.$LONGITUDE,0);
      this.header.distance_X = this.header.VN2000_X - objects.entities[objects.entities.length-1].position.x;
      this.header.distance_Y = this.header.VN2000_Y - objects.entities[objects.entities.length-1].position.y;
      this.EXTMIN = ConvertGeo(this.header.$EXTMIN.y,this.header.$EXTMIN.x);
      this.EXTMAX = ConvertGeo(this.header.$EXTMAX.y,this.header.$EXTMAX.x);
    }

    if(objects.tables !== undefined){
      this.tables = objects.tables;
    }

    if(objects.blocks !== undefined)
      this.blocks = objects.blocks;
    
    
    objects = objects.entities;
    
    window.dcel = new DCEL();

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
  var geoms = {Point:['POINT'], MultiPoint:[], LineString:['LINE','LWPOLYLINE','SPLINE'], MultiLineString:[], Polygon:[], MultiPolygon:[], GeoJSON:[]},
      geomAttrs = [],
      geoInput = { OLEFRAME:'',OLE2FRAME:'',ACAD_PROXY_ENTITY:'',POINT:'',ARC:'',POLYLINE:'',ARCALIGNEDTEXT:'',RAY:'',ATTDEF:'',REGION:'',ATTRIB:'',RTEXT:'',BODY:'',SEQEND:'',CIRCLE:'',SHAPE:'',DIMENSION:'',SOLID:'',ELLIPSE:'',SPLINE:{controlPoints:['x','y']},HATCH:'',TEXT:'',IMAGE:'',TOLERANCE:'',INSERT:"",TRACE:'',VERTEX:['x','y'],LINE:{vertices:['x','y']},VIEWPORT:'',LWPOLYLINE:{vertices:['x','y']},WIPEOUT:'',MLINE:'',XLINE:'',MTEXT:''};
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
    if(feature.geometry!==undefined)
    if(feature.geometry.type==="LineString"){
      for(var k=1;k<feature.geometry.coordinates.length;k++)
      {
        window.dcel.stackLine([feature.geometry.coordinates[k],feature.geometry.coordinates[k-1]]);
      }
    }
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
