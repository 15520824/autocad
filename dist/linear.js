(function () {
    window.linear = {
        checkLineIntersection : function(line1StartX, line1StartY, line1EndX, line1EndY, line2StartX, line2StartY, line2EndX, line2EndY) {
            // console.log(line1StartX, line1StartY,",,,,,,,,,,,,,", line1EndX, line1EndY,"|||||||||||||", line2StartX, line2StartY,",,,,,,,,,,,,,", line2EndX, line2EndY)
            // if the lines intersect, the result contains the x and y of the intersection (treating the lines as infinite) and booleans for whether line segment 1 or line segment 2 contain the point
            var denominator, a, b, numerator1, numerator2, result = [null,null];
            denominator = ((line2EndY - line2StartY) * (line1EndX - line1StartX)) - ((line2EndX - line2StartX) * (line1EndY - line1StartY));
            if (denominator == 0) {
                if(this.lineIntersection(line1StartX, line1StartY, line1EndX, line1EndY, line2StartX, line2StartY, line2EndX, line2EndY))
                        return result;
                return false;
            }
            a = line1StartY - line2StartY;
            b = line1StartX - line2StartX;
            numerator1 = ((line2EndX - line2StartX) * a) - ((line2EndY - line2StartY) * b);
            numerator2 = ((line1EndX - line1StartX) * a) - ((line1EndY - line1StartY) * b);
            a = numerator1 / denominator;
            b = numerator2 / denominator;
            // if line1 is a segment and line2 is infinite, they intersect if:
            if (!(a >= 0 && a <= 1)) {
                return false;
            }
                
            // if line2 is a segment and line1 is infinite, they intersect if:
            if (!(b >= 0 && b <= 1)) {
                return false;
            }
        
            // if we cast these lines infinitely in both directions, they intersect here:
            result[0] = line1StartX + (a * (line1EndX - line1StartX));
            result[1] = line1StartY + (a * (line1EndY - line1StartY));
            
        /*
                // it is worth noting that this should be the same as:
                x = line2StartX + (b * (line2EndX - line2StartX));
                y = line2StartX + (b * (line2EndY - line2StartY));
                */
            
            // if line1 and line2 are segments, they intersect if both of the above are true
            return result;
        },
        lineIntersection: function(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y) {

            var s1_x, s1_y, s2_x, s2_y;
            s1_x = p1_x - p0_x;
            s1_y = p1_y - p0_y;
            s2_x = p3_x - p2_x;
            s2_y = p3_y - p2_y;
            
            var s, t;
            s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
            t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
            
            m1 = s1_y/s1_x; 
            m2 = s2_y/s2_x;
            
            //If they have the same slope check for the points to intersect
            if(m1 == m2)
                return ((p0_x - p2_x)*(p0_y - p3_y) - (p0_x - p3_x)*(p0_y - p2_y) == 0 || (p1_x - p2_x)*(p1_y - p3_y) - (p1_x - p3_x)*(p1_y - p2_y) == 0)|| 
                (p2_x - p0_x)*(p2_y - p1_y) - (p2_x - p1_x)*(p2_y - p0_y) == 0||(p3_x - p0_x)*(p3_y - p1_y) - (p3_x - p1_x)*(p3_y - p0_y) == 0;
            return (s >= 0 && s <= 1 && t >= 0 && t <= 1);
            }
    }
}());