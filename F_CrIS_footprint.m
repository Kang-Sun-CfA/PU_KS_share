function outp = F_CrIS_footprint(inp)

% calculate CrIS footprint accurately using geometry. 
% Repackaged by Kang Sun on 2018/01/23

outp = [];

geodeticLat = inp.latc; 
longitude = inp.lonc;
range = inp.range;
azimuth = inp.azi;
zenith = inp.zen;
fovDia = inp.fovDia;

if ~isfield(inp,'nanglestep')
    nanglestep = 37;
else
    nanglestep = inp.nanglestep;
end
%
% Description:
%
%   Given a Earth radius, a reference location  (latRef, lonRef) in
%   geographical coordinate, and a location of interest( lat, lon), 
%   this function return the distance between the 2 locations.
%
% Usage:
%
%   [distance] = spherical_distance( lat1, lon1, lat2, lon2, radius)
%
% Input:
%
%     lat1    : latitude  of location 1 ( geographical coordinates).
%     lon1    : longitude of location 1 ( geographical coordinates).
%     lat2    : latitude  of location 2 ( geographical coordinates).
%     lon2    : longitude of location 2 ( geographical coordinates).
%     radius  : Earth radius in Km
%
% Output:
%
%     distance :  Distance between location 1 and 2, same units as
%                 the given radius.
%
% Prerequisite
%
%   
% Called routines
%
%
% Authors: 
%   Denis Tremblay, Science Data Processing Inc.
%
% Change Record:
%   Date         By   Description
%   23-Sep-2010  DAT  Initial RCS version.
%
% COPYRIGHT (C) SDPI
%
%   Science Data Processing Inc.
%


   earthRadius = 6378137.0e0;
   flatFact = 1.0e0 / 298.257223563e0;
   deg2rad = pi / 180.0e0;

   geocentricPoint = zeros( 1, 3);
   geodeticEquatorialVector = zeros( 1, 3);
   lonRad = longitude * deg2rad;
   azimuthRad = azimuth * deg2rad;
   zenithRad = zenith * deg2rad;
   geodeticLatRad = geodeticLat * deg2rad;
   fovRadius = fovDia / 2.0e0;  % here fov diameter is in radians.

   %  Find the geocentric latitude.
   %  1 = geodetic to geocentric, 2 = geocentric to geodetic.
   %
   geocentricLat = latitude_transform( geodeticLat, flatFact, 1); 
   geocentricLatRad = geocentricLat * deg2rad ; 

   geocentricRange  = geocentric_range(geocentricLat, earthRadius, flatFact);  

   geocentricTempX = geocentricRange * cos(geocentricLatRad);
   geocentricPoint(3) = geocentricRange * sin(geocentricLatRad);
   geocentricPoint(1) = geocentricRange * cos(geocentricLatRad) * cos( lonRad);
   geocentricPoint(2) = geocentricRange * cos(geocentricLatRad) * sin( lonRad);


   %  Find the satellite line-of-sight ( look vector) assuming
   %  that the geolocation point is (1, 0, 0). The next step will be
   %  to perform 2 rotation to find the LOS in the geocentric frame
   %  or Earth Centered Reference (ECR).
   %  
   %  Here, the X axis is outward at longitude zero, Z is toward the North pole
   %  and Y is orthogonal to X and Z.
   %
   losPointX = zeros(1, 3);
   losPointX = [ -1.0e0, 0.0e0, 0.0e0];

   rotationVector = [ 0.0e0, -1.0e0, 0.0e0];
   losPointX1 = rotate_vector( rotationVector, zenithRad, losPointX); 
   rotationVector = [ -1.0e0, 0.0e0, 0.0e0];
   losPointX2 = rotate_vector( rotationVector, azimuthRad, losPointX1); 

   %  Now find the LOS in the geocentric frame. This is done by
   %  performing 2 rotations: 
   %   1) geodetic latitude around [0,-1,0],
   %   2) longitude around [ 0, 0, 1].
   %

%keyboard
   rotationVector = [ 0.0e0, -1.0e0, 0.0e0];
   losPointX3 = rotate_vector( rotationVector, geodeticLatRad, losPointX2); 
   rotationVector = [ 0.0e0, 0.0e0, 1.0e0];

   if( lonRad < 0.0)
      useLonRad = lonRad + 2.0*pi;
   else
      useLonRad = lonRad;
   end
   losPointX4 = rotate_vector( rotationVector, useLonRad, losPointX3); 

%keyboard

   %  Find the satellite position in the ECR.
   %
   satellitePosition = geocentricPoint - (losPointX4 * range);



   %  Find an orthogonal vector to the LOS.
   %  Here, it does not matter the direction.
   %
   orthoVectorLOS = cross( losPointX4, [ 0.0e0, 0.0e0, 1.0e0]);
   orthoVectorLOS = orthoVectorLOS / ( sqrt( dot(orthoVectorLOS, orthoVectorLOS))); 

   %  Rotatate the orthoVector to LOS by the FOV radius ( in radians).
   %  
   fovVector = rotate_vector( orthoVectorLOS, fovRadius, losPointX4);


   %  Rotate the fovVector in step of 10 degrees, and find the ellipsoid
   %  geolocation ( in geodetic latitude/longitude coordiante).
   %  The final result is a vector of longitude ( 37 data) and 
   %  latitude ( 37 entries). The last entry (index 37) is identical
   %  to the first index, hence allowing the plotting using plotm routine.
   % 
   angleStep = (0:36) * (10.0e0 * deg2rad);

   ellLat = zeros( 1, 37);
   ellLon = zeros( 1, 37);

   for iAngleStep = 1: nanglestep

      curAngle = angleStep(iAngleStep);
      curFovVector = rotate_vector( losPointX4, curAngle, fovVector);

      %  Find the geodetic latitude and longitude ( geolocation)
      %  for the given current fov Vector and satellite position.
      %
      [ geolocationLat, geolocationLon ] = compute_geolocation( satellitePosition, curFovVector, earthRadius, flatFact);

%format long eng
%satellitePosition
%losPointX4
%range
%lonRad
%geocentricLat

%      [ geolocationLat, geolocationLon ] = compute_geolocation( satellitePosition, losPointX4, earthRadius, flatFact);
      ellLat(iAngleStep) = geolocationLat;
      ellLon(iAngleStep) = geolocationLon;

   end

   outp.latr = ellLat;
   outp.lonr = ellLon;

return;

function [geoLat, geoLon] = compute_geolocation(satellitePosition, lineOfSight, earthRadius, flatFact  )
%
% Description:
%
%   Given the satellite position and the line-of-sight (or lookup vector), this function returns
%   the geodetic (or geographic) latitude and longitude of the geolocation point on the Earth
%   ellipsoid characterized by the Earth radius and flattening factor. 
%
% Usage:
%
%   [geoLat, geoLon] = compute_geolocation(satellitePosition, lineOfSight, earthRadius, flatFact  )
%
% Input:
%
%     satellitePosition
%             : Three elements vector of the satellite position in the Earth Centered
%               Reference (ECR) frame (float or double)
%     lineOfSight
%             : Three elements unit vector of lookup vector ( viewing vector toward the Earth surface).
%               in the Earth Centered Reference (ECR) frame (float or double)
%     earthRadius
%             : Earth equatorial radius in meters (float/double). For WGS84, the value
%               is 6378137.0 meters.
%     flatFact: Earth ellipsoid flattening factor. For WGS84, the value is 1/298.257223563 ( unitless).
%
% Output:
%
%     geoLat  :  Geodetic (geographic) latitude of the Earth geolocation point. 
%     geoLon  :  Longitude of the Earth geolocation point. 
%
% Prerequisite
%
%   
% Called routines
%
%
% Authors: 
%   Denis Tremblay, Science Data Processing Inc.
%
% Change Record:
%   Date         By   Description
%   12-Oct-2010  DAT  Initial RCS version.
%
% COPYRIGHT (C) SDPI
%
%   Science Data Processing Inc.
%
 
   %  The basic equations are:
   %
   %  P + lambda LOS = G  ( P = satellite position, LOS = line of sight, G = geolocation point).
   %  and lambda is the slant range.
   %
   %  and
   %
   %  Gx^2 / a^2 + Gy^2 / a^2 + Gz^2 / c^2 = 1    Earth ellipsoid equation, a = equatorial radius.
   %  c is polar radius ( c = ( 1-f) * a) where f is flattening factor.
   %

   %  Initialize working variable.
   %
   geolocationPoint = zeros( 1, 3);
   deg2rad = pi / 180.0;

   %  The geolocation vector position is the solution of
   %  a quadratic equation where  A lambda^2 + B lambda + C = 0, here x is the slant range.
   %
   polarRadius = earthRadius * ( 1.0e0 - flatFact);
   termA = ((lineOfSight(1) / earthRadius)^2) +...
           ((lineOfSight(2) / earthRadius)^2) +...
           ((lineOfSight(3) / polarRadius)^2);

   termB =  ( satellitePosition(1) * lineOfSight(1) / (earthRadius^2)) +...
            ( satellitePosition(2) * lineOfSight(2) / (earthRadius^2)) +...
            ( satellitePosition(3) * lineOfSight(3) / (polarRadius^2))  ;
   termB = termB * 2.0e0;

   termC =  (satellitePosition(1)/earthRadius)^2 + ...
            (satellitePosition(2)/earthRadius)^2 + ...
            (satellitePosition(3)/polarRadius)^2 - 1.0e0;

   radical = termB^2 - (4.0e0 * termA * termC);

   if( radical < 0.0e0)
      %  The line of sight does not intercept the Earth ellipsoid.
      %
      ellLat = -999.0;
      ellLon = -999.0;
   end
   if( radical == 0.0e0)
      %  The line of sight does not intercept the Earth ellipsoid tangentially.
      %
      slantRange = -termB / (2.0 * termA);
      geolocationPoint = satellitePosition + slantRange * lineOfSight;
   end
   if( radical > 0.0e0)
      %  The line of sight intercepts the Earth ellipsoid at 2 point, the solution
      %  is the shorter slant range.
      %
      slantRange1 = (-termB - sqrt(radical))  / (2.0 * termA);
      slantRange2 = (-termB + sqrt(radical))  / (2.0 * termA);
      slantRange = min( [ slantRange1, slantRange2]);
      geolocationPoint = satellitePosition + slantRange * lineOfSight;
   end

   geoLon = atan2( geolocationPoint(2) , geolocationPoint(1)) / deg2rad;
   magXY = sqrt( geolocationPoint(1)^2 +  geolocationPoint(2)^2); 
   geocentricLat = atan2( geolocationPoint(3) , magXY) / deg2rad;
   geoLat = latitude_transform( geocentricLat, flatFact, 2);
            

return;

function [range] = geocentric_range(geocentricLatitude, equatorialRadius, flatteningFactor)
%
% Description:
%
%   Given a Earth equatorial radius, the geocentric latitude, and
%   the flattening factor, find the geocentric range.
%
% Usage:
%
%   [range] = geocentric_range( latitude, equatorialRadius, flatteningFactor )
%
% Input:
%
%     latitude    : Geocentric latitude  of location 1, in degrees ( geographical coordinates).
%     equatorialRadius
%                 : Earth equaorial radius in meters ( WGS84 is 6378137.0 meters).
%     flatteningFactor
%                 : Flattening factor of the Earth, WGS84 value is 1/298.257223560
%
% Output:
%
%     range       :  Geocetric distance at the given geocentric latitude.
%
% Prerequisite
%
%   
% Called routines
%
%
% Authors: 
%   Denis Tremblay, Science Data Processing Inc.
%
% Change Record:
%   Date         By   Description
%   23-Sep-2010  DAT  Initial RCS version.
%
% COPYRIGHT (C) SDPI
%
%   Science Data Processing Inc.
%

   deg2rad = pi / 180.0e0;
   polarRadius = (1.0e0 - flatteningFactor) * equatorialRadius;
   lat = geocentricLatitude * deg2rad;
   range = sqrt( (equatorialRadius * cos(lat))^2 + (polarRadius * sin( lat))^2);

return;

function [latitudeOutput] = latitude_transform( latRef, flatteningFactor, option)
%
% Description:
%
%   Given a latitude, transform into geographic( geodetic) or into
%   geocentric.
%
% Usage:
%
%   [latitudeOutput] = lattitude_transform( latitudeInput, flatteningFactor, option)
%
% Input:
%
%     latitudeInput
%             : Latitude input
%     flatteningFactor
%             : flattening factor defined as (a-c)/a where a is the equatorial
%               radius, and c is the polar radius. The WGS84 values is
%               1/298.2573223560
%     option  : Transformation : 1 = geodetic to geocentric
%                                2 = geocentric to geodetic
%
% Output:
%
%     latitudeOutput
%             : Latitude output
%
% Prerequisite
%
%   
% Called routines
%
%
% Authors: 
%   Denis Tremblay, Science Data Processing Inc.
%
% Change Record:
%   Date         By   Description
%   23-Sep-2010  DAT  Initial RCS version.
%
% COPYRIGHT (C) SDPI
%
%   Science Data Processing Inc.
%
   deg2rad = pi / 180.0;
   latIn = latRef * deg2rad;
   ff = 1.0 - flatteningFactor;

   if( option == 1)

      latOut = atan( ((ff)^2) * tan(latIn));
      latitudeOutput = latOut / deg2rad;
   else

      latOut = atan( (1.0/ (ff^2)) * tan(latIn));
      latitudeOutput = latOut / deg2rad;

   end

return;

function [newVector] = rotate_vector( rotationAxis, angle, oldVector )
%
% Description:
%
%   Given an axis or rotation, the angle (right hand screw motion), and
%   the original ( old) vector, this function rotates the old vector 
%   around the rotation axis by the amount of specified angle. 
%   geocentric.
%
% Usage:
%
%   [newVector] = rotate_vector( rotationAxis, angle, oldVector )
%
% Input:
%
%     rotationAxis
%             : Three elements unit vector of the rotation axis ( cartesian coordinate)
%               (float or double)
%               Example, rotate around the Y axis, then the rotation axis is [ 0, 1, ].
%     oldVector
%             : Three elements vector to be rotated  ( cartesian coordinate)
%               (float or double)
%     angle   : Angle of ratation in radian.
%               (float or double)
%
% Output:
%
%     newVector
%             : Three elements vector after rotation.
%               (float or double)
%
% Prerequisite
%
%   
% Called routines
%
%
% Authors: 
%   Denis Tremblay, Science Data Processing Inc.
%
% Change Record:
%   Date         By   Description
%   12-Oct-2010  DAT  Initial RCS version.
%
% COPYRIGHT (C) SDPI
%
%   Science Data Processing Inc.
%

   cosAngle = cos( angle);
   sinAngle = sin( angle);

   firstTerm = oldVector * cosAngle;
   secondTerm =  cross( oldVector, rotationAxis) * sinAngle;
   thirdTerm = (dot( oldVector, rotationAxis) * rotationAxis) * ( 1.0e0 - cosAngle);
   newVector = firstTerm - secondTerm + thirdTerm;
return;


