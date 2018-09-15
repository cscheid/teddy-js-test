(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global.Delaunator = factory());
}(this, (function () { 'use strict';

    var EPSILON = Math.pow(2, -52);

    var Delaunator = function Delaunator(coords) {
        var this$1 = this;

        var minX = Infinity;
        var minY = Infinity;
        var maxX = -Infinity;
        var maxY = -Infinity;

        var n = coords.length >> 1;
        var ids = this.ids = new Uint32Array(n);

        if (n > 0 && typeof coords[0] !== 'number') { throw new Error('Expected coords to contain numbers.'); }

        this.coords = coords;

        for (var i = 0; i < n; i++) {
            var x = coords[2 * i];
            var y = coords[2 * i + 1];
            if (x < minX) { minX = x; }
            if (y < minY) { minY = y; }
            if (x > maxX) { maxX = x; }
            if (y > maxY) { maxY = y; }
            ids[i] = i;
        }

        var cx = (minX + maxX) / 2;
        var cy = (minY + maxY) / 2;

        var minDist = Infinity;
        var i0, i1, i2;

        // pick a seed point close to the centroid
        for (var i$1 = 0; i$1 < n; i$1++) {
            var d = dist(cx, cy, coords[2 * i$1], coords[2 * i$1 + 1]);
            if (d < minDist) {
                i0 = i$1;
                minDist = d;
            }
        }
        var i0x = coords[2 * i0];
        var i0y = coords[2 * i0 + 1];

        minDist = Infinity;

        // find the point closest to the seed
        for (var i$2 = 0; i$2 < n; i$2++) {
            if (i$2 === i0) { continue; }
            var d$1 = dist(i0x, i0y, coords[2 * i$2], coords[2 * i$2 + 1]);
            if (d$1 < minDist && d$1 > 0) {
                i1 = i$2;
                minDist = d$1;
            }
        }
        var i1x = coords[2 * i1];
        var i1y = coords[2 * i1 + 1];

        var minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (var i$3 = 0; i$3 < n; i$3++) {
            if (i$3 === i0 || i$3 === i1) { continue; }
            var r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i$3], coords[2 * i$3 + 1]);
            if (r < minRadius) {
                i2 = i$3;
                minRadius = r;
            }
        }
        var i2x = coords[2 * i2];
        var i2y = coords[2 * i2 + 1];

        if (minRadius === Infinity) {
            throw new Error('No Delaunay triangulation exists for this input.');
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
            var i$4 = i1;
            var x$1 = i1x;
            var y$1 = i1y;
            i1 = i2;
            i1x = i2x;
            i1y = i2y;
            i2 = i$4;
            i2x = x$1;
            i2y = y$1;
        }

        var center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
        this._cx = center.x;
        this._cy = center.y;

        // sort the points by distance from the seed triangle circumcenter
        quicksort(ids, coords, 0, ids.length - 1, center.x, center.y);

        // initialize a hash table for storing edges of the advancing convex hull
        this._hashSize = Math.ceil(Math.sqrt(n));
        this._hash = new Int32Array(this._hashSize).fill(-1);

        // initialize arrays for tracking the edges of the advancing convex hull
        var hullPrev = this.hullPrev = new Uint32Array(n); // edge to prev edge
        var hullNext = this.hullNext = new Uint32Array(n); // edge to next edge
        var hullTri = this.hullTri = new Uint32Array(n); // edge to adjacent triangle

        // set up the seed triangle as the starting hull
        this.hullStart = i0;
        var hullSize = 3;

        hullNext[i0] = hullPrev[i2] = i1;
        hullNext[i1] = hullPrev[i0] = i2;
        hullNext[i2] = hullPrev[i1] = i0;

        hullTri[i0] = 0;
        hullTri[i1] = 1;
        hullTri[i2] = 2;

        this._hash[this._hashKey(i0x, i0y)] = i0;
        this._hash[this._hashKey(i1x, i1y)] = i1;
        this._hash[this._hashKey(i2x, i2y)] = i2;

        var maxTriangles = 2 * n - 5;
        var triangles = this.triangles = new Uint32Array(maxTriangles * 3);
        var halfedges = this.halfedges = new Int32Array(maxTriangles * 3);

        this.trianglesLen = 0;

        this._addTriangle(i0, i1, i2, -1, -1, -1);

        for (var k = 0, xp = (void 0), yp = (void 0); k < ids.length; k++) {
            var i$5 = ids[k];
            var x$2 = coords[2 * i$5];
            var y$2 = coords[2 * i$5 + 1];

            // skip near-duplicate points
            if (k > 0 && Math.abs(x$2 - xp) <= EPSILON && Math.abs(y$2 - yp) <= EPSILON) { continue; }
            xp = x$2;
            yp = y$2;

            // skip seed triangle points
            if (i$5 === i0 || i$5 === i1 || i$5 === i2) { continue; }

            // find a visible edge on the convex hull using edge hash
            var start = 0;
            for (var j = 0, key = this._hashKey(x$2, y$2); j < this._hashSize; j++) {
                start = this$1._hash[(key + j) % this$1._hashSize];
                if (start !== -1 && start !== hullNext[start]) { break; }
            }

            start = hullPrev[start];
            var e = start, q = (void 0);
            while (q = hullNext[e], !orient(x$2, y$2, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
                e = q;
                if (e === start) {
                    e = -1;
                    break;
                }
            }
            if (e === -1) { continue; } // likely a near-duplicate point; skip it

            // add the first triangle from the point
            var t = this$1._addTriangle(e, i$5, hullNext[e], -1, -1, hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i$5] = this$1._legalize(t + 2);
            hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize++;

            // walk forward through the hull, adding more triangles and flipping recursively
            var n$1 = hullNext[e];
            while (q = hullNext[n$1], orient(x$2, y$2, coords[2 * n$1], coords[2 * n$1 + 1], coords[2 * q], coords[2 * q + 1])) {
                t = this$1._addTriangle(n$1, i$5, q, hullTri[i$5], -1, hullTri[n$1]);
                hullTri[i$5] = this$1._legalize(t + 2);
                hullNext[n$1] = n$1; // mark as removed
                hullSize--;
                n$1 = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e === start) {
                while (q = hullPrev[e], orient(x$2, y$2, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
                    t = this$1._addTriangle(q, i$5, e, -1, hullTri[e], hullTri[q]);
                    this$1._legalize(t + 2);
                    hullTri[q] = t;
                    hullNext[e] = e; // mark as removed
                    hullSize--;
                    e = q;
                }
            }

            // update the hull indices
            this$1.hullStart = hullPrev[i$5] = e;
            hullNext[e] = hullPrev[n$1] = i$5;
            hullNext[i$5] = n$1;

            // save the two new edges in the hash table
            this$1._hash[this$1._hashKey(x$2, y$2)] = i$5;
            this$1._hash[this$1._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
        }

        this.hull = new Uint32Array(hullSize);
        for (var i$6 = 0, e$1 = this.hullStart; i$6 < hullSize; i$6++) {
            this$1.hull[i$6] = e$1;
            e$1 = hullNext[e$1];
        }
        this.hullPrev = this.hullNext = this.hullTri = null; // get rid of temporary arrays

        // trim typed triangle mesh arrays
        this.triangles = triangles.subarray(0, this.trianglesLen);
        this.halfedges = halfedges.subarray(0, this.trianglesLen);
    };

    Delaunator.from = function from (points, getX, getY) {
        if (!getX) { getX = defaultGetX; }
        if (!getY) { getY = defaultGetY; }

        var n = points.length;
        var coords = new Float64Array(n * 2);

        for (var i = 0; i < n; i++) {
            var p = points[i];
            coords[2 * i] = getX(p);
            coords[2 * i + 1] = getY(p);
        }

        return new Delaunator(coords);
    };

    Delaunator.prototype._hashKey = function _hashKey (x, y) {
        return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
    };

    Delaunator.prototype._legalize = function _legalize (a) {
            var this$1 = this;

        var ref = this;
            var triangles = ref.triangles;
            var coords = ref.coords;
            var halfedges = ref.halfedges;

        var b = halfedges[a];

        /* if the pair of triangles doesn't satisfy the Delaunay condition
         * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
         * then do the same check/flip recursively for the new pair of triangles
         *
         *       pl                pl
         *      /||\              /  \
         *   al/ || \bl        al/\a
         *    /  ||  \          /  \
         *   /  a||b  \flip/___ar___\
         * p0\   ||   /p1   =>   p0\---bl---/p1
         *    \  ||  /          \  /
         *   ar\ || /br         b\/br
         *      \||/              \  /
         *       pr                pr
         */
        var a0 = a - a % 3;
        var b0 = b - b % 3;

        var al = a0 + (a + 1) % 3;
        var ar = a0 + (a + 2) % 3;
        var bl = b0 + (b + 2) % 3;

        if (b === -1) { return ar; }

        var p0 = triangles[ar];
        var pr = triangles[a];
        var pl = triangles[al];
        var p1 = triangles[bl];

        var illegal = inCircle(
            coords[2 * p0], coords[2 * p0 + 1],
            coords[2 * pr], coords[2 * pr + 1],
            coords[2 * pl], coords[2 * pl + 1],
            coords[2 * p1], coords[2 * p1 + 1]);

        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            var hbl = halfedges[bl];

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if (hbl === -1) {
                var e = this.hullStart;
                do {
                    if (this$1.hullTri[e] === bl) {
                        this$1.hullTri[e] = a;
                        break;
                    }
                    e = this$1.hullNext[e];
                } while (e !== this.hullStart);
            }
            this._link(a, hbl);
            this._link(b, halfedges[ar]);
            this._link(ar, bl);

            var br = b0 + (b + 1) % 3;

            this._legalize(a);
            return this._legalize(br);
        }

        return ar;
    };

    Delaunator.prototype._link = function _link (a, b) {
        this.halfedges[a] = b;
        if (b !== -1) { this.halfedges[b] = a; }
    };

    // add a new triangle given vertex indices and adjacent half-edge ids
    Delaunator.prototype._addTriangle = function _addTriangle (i0, i1, i2, a, b, c) {
        var t = this.trianglesLen;

        this.triangles[t] = i0;
        this.triangles[t + 1] = i1;
        this.triangles[t + 2] = i2;

        this._link(t, a);
        this._link(t + 1, b);
        this._link(t + 2, c);

        this.trianglesLen += 3;

        return t;
    };

    // monotonically increases with real angle, but doesn't need expensive trigonometry
    function pseudoAngle(dx, dy) {
        var p = dx / (Math.abs(dx) + Math.abs(dy));
        return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
    }

    function dist(ax, ay, bx, by) {
        var dx = ax - bx;
        var dy = ay - by;
        return dx * dx + dy * dy;
    }

    function orient(px, py, qx, qy, rx, ry) {
        return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0;
    }

    function inCircle(ax, ay, bx, by, cx, cy, px, py) {
        var dx = ax - px;
        var dy = ay - py;
        var ex = bx - px;
        var ey = by - py;
        var fx = cx - px;
        var fy = cy - py;

        var ap = dx * dx + dy * dy;
        var bp = ex * ex + ey * ey;
        var cp = fx * fx + fy * fy;

        return dx * (ey * cp - bp * fy) -
               dy * (ex * cp - bp * fx) +
               ap * (ex * fy - ey * fx) < 0;
    }

    function circumradius(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = dx * ey - dy * ex;

        var x = (ey * bl - dy * cl) * 0.5 / d;
        var y = (dx * cl - ex * bl) * 0.5 / d;

        return bl && cl && d && (x * x + y * y) || Infinity;
    }

    function circumcenter(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = dx * ey - dy * ex;

        var x = ax + (ey * bl - dy * cl) * 0.5 / d;
        var y = ay + (dx * cl - ex * bl) * 0.5 / d;

        return {x: x, y: y};
    }

    function quicksort(ids, coords, left, right, cx, cy) {
        var i, j, temp;

        if (right - left <= 20) {
            for (i = left + 1; i <= right; i++) {
                temp = ids[i];
                j = i - 1;
                while (j >= left && compare(coords, ids[j], temp, cx, cy) > 0) { ids[j + 1] = ids[j--]; }
                ids[j + 1] = temp;
            }
        } else {
            var median = (left + right) >> 1;
            i = left + 1;
            j = right;
            swap(ids, median, i);
            if (compare(coords, ids[left], ids[right], cx, cy) > 0) { swap(ids, left, right); }
            if (compare(coords, ids[i], ids[right], cx, cy) > 0) { swap(ids, i, right); }
            if (compare(coords, ids[left], ids[i], cx, cy) > 0) { swap(ids, left, i); }

            temp = ids[i];
            while (true) {
                do { i++; } while (compare(coords, ids[i], temp, cx, cy) < 0);
                do { j--; } while (compare(coords, ids[j], temp, cx, cy) > 0);
                if (j < i) { break; }
                swap(ids, i, j);
            }
            ids[left + 1] = ids[j];
            ids[j] = temp;

            if (right - i + 1 >= j - left) {
                quicksort(ids, coords, i, right, cx, cy);
                quicksort(ids, coords, left, j - 1, cx, cy);
            } else {
                quicksort(ids, coords, left, j - 1, cx, cy);
                quicksort(ids, coords, i, right, cx, cy);
            }
        }
    }

    function compare(coords, i, j, cx, cy) {
        var d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
        var d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
        return (d1 - d2) || (coords[2 * i] - coords[2 * j]) || (coords[2 * i + 1] - coords[2 * j + 1]);
    }

    function swap(arr, i, j) {
        var tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    function defaultGetX(p) {
        return p[0];
    }
    function defaultGetY(p) {
        return p[1];
    }

    return Delaunator;

})));
