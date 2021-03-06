
function distance(p1, p2) {
    var dx = p2[0] - p1[0], dy = p2[1] - p1[1];
    return Math.sqrt(dx * dx + dy * dy);
}

function subtractPoint(p1, p2) {
    return [p1[0] - p2[0], p1[1] - p2[1]];
}

function dotProduct(p1, p2) {
    return p1[0] * p2[0] + p1[1] * p2[1];
}

function crossProduct(p1, p2) {
    return p1[0] * p2[1] - p1[1] * p2[0];
}

function vecLength(v) {
    return distance(v, v);
}

function nextVertexInTri(i) {
    if (i % 3 === 2)
        return i - 2;
    else
        return i + 1;
}

function pointFromTri(del, pointId) {
    return [del.coords[2*pointId], del.coords[2*pointId+1]];
}

function getDelaunayTri(del, triId) {
    var p1 = pointFromTri(del, del.triangles[3*triId]);
    var p2 = pointFromTri(del, del.triangles[3*triId+1]);
    var p3 = pointFromTri(del, del.triangles[3*triId+2]);
    return [p1, p2, p3];
}

function barycenterOfTri(p1, p2, p3) {
    return [(p1[0] + p2[0] + p3[0]) / 3, (p1[1] + p2[1] + p3[1]) / 3];
}

function signedArea(p1, p2, p3) {
    var s1 = subtractPoint(p2, p1); // p1->p2
    var s2 = subtractPoint(p1, p3);
    return crossProduct(s1, s2) / 2;
}

// pointInsideTri assumes p1, p2, p3 in CCW order.
function pointInsideTri(p1, p2, p3, p) {
    if (signedArea(p1, p2, p) < 0) return false;
    if (signedArea(p2, p3, p) < 0) return false;
    if (signedArea(p3, p1, p) < 0) return false;
    return true;
}

// https://en.wikipedia.org/wiki/Circumscribed_circle#Circumcenter_coordinates
function circumcenter(a, b, c) {
    var ad = a[0] * a[0] + a[1] * a[1];
    var bd = b[0] * b[0] + b[1] * b[1];
    var cd = c[0] * c[0] + c[1] * c[1];
    var d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    var cx = 1 / d * (ad * (b[1] - c[1]) + bd * (c[1] - a[1]) + cd * (a[1] - b[1]));
    var cy = 1 / d * (ad * (c[0] - b[0]) + bd * (a[0] - c[0]) + cd * (b[0] - a[0]));
    return [
        cx, cy, distance([cx, cy], a)
    ];
}

// barycenterInsideShape walks the triangle mesh monotonically until
// finding a shape edge, then checks a local condition.

// barycenterInsideShape assumes that the shape is implicitly given by
// the edges 0-1-2-3-4-5-6-7-8-...-n-0, and the shape is given in CCW
// order.

function barycenterInsideShape(del, pts, i)
{
    // first: is this an easy triangle?
    var nPoints = pts.length / 2;
    var di1 = del.triangles[3*i+1] - del.triangles[3*i];
    var di2 = del.triangles[3*i+2] - del.triangles[3*i+1];
    var di3 = del.triangles[3*i  ] - del.triangles[3*i+2];
    if (di1 === 1 || di2 === 1 || di3 === 1 ||
        di1 === pts.length-1 || di2 === pts.length-1 || di3 === pts.length-1) {
        // we found a CCW edge: we're inside the shape
        return true;
    }
    if (di1 === -1 || di2 === -1 || di3 === -1 ||
        di1 === 1-pts.length || di2 === 1-pts.length || di3 === 1-pts.length) {
        // we found a CW edge: we're outside the shape
        return false;
    }

    // else, walk to the triangle that increases the x coordinate
    // the barycenter, repeat.

    var thisTriangle = getDelaunayTri(del, i);
    var thisBarycenter = barycenterOfTri(
        thisTriangle[0], thisTriangle[1], thisTriangle[2]);
    // find a half-edge that's available
    for (var h=0; h<3; ++h) {
        var he = del.halfedges[3*i+h];
        if (he === -1)
            continue;
        var triangleOfHalfEdge = ~~(he / 3);
        var neighborTri = getDelaunayTri(del, triangleOfHalfEdge);
        var neighborBarycenter = barycenterOfTri(
            neighborTri[0], neighborTri[1], neighborTri[2]);
        if (neighborBarycenter[0] > thisBarycenter[0]) {
            return barycenterInsideShape(del, pts, triangleOfHalfEdge);
        }
    }

    // couldn't find an available half-edge, that means we ran into
    // the convex hull of the DT without crossing the shape, so we
    // must have been outside to begin with

    return false;
}

function determineTriangleSides(del, pts)
{
    var triangleSides = d3.range(del.triangles.length / 3).map((i) => "unknown");
    var nPts = pts.length / 2;
    var queue = d3.range(del.triangles.length / 3);
    while (queue.length > 0) {
        var triId = queue.pop();
        if (triangleSides[triId] !== "unknown")
            continue;
        // we don't know about this one, so let's run the slow test.
        var inside = barycenterInsideShape(del, pts, triId);

        // now we flood-fill the graph with the result
        var stack = [triId];
        while (stack.length > 0) {
            var triIndex = stack.pop();
            if (triangleSides[triIndex] !== "unknown") {
                if (triangleSides[triIndex] !== inside) {
                    console.error("Doesn't look like we've sampled this shape enough.");
                    debugger;
                }
                continue;
            }
            triangleSides[triIndex] = inside;
            // find neighbor triangles
            for (var h=0; h<3; ++h) {
                var dIndex = Math.abs(del.triangles[3*triIndex+h%3] -
                                      del.triangles[3*triIndex+(h+1)%3]);
                if (dIndex == 1 || dIndex == nPts - 1)
                    continue; // we hit a shape edge
                var he = del.halfedges[3*triIndex+h];
                if (he === -1)
                    continue; // we hit a boundary of the DT
                var triangleOfHalfEdge = ~~(he / 3);
                if (triangleSides[triangleOfHalfEdge] === "unknown") {
                    stack.push(triangleOfHalfEdge);
                } else if (triangleSides[triangleOfHalfEdge] !== inside) {
                    console.error("Doesn't look like we've sampled this shape enough.");
                    debugger;
                }
            }
        }
    }

    return triangleSides;
}
