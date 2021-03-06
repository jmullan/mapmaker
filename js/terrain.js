"use strict";

var INF = 99999999,
    MAX_HEIGHT = 9000,
    MIN_HEIGHT = -11000;

function randBetween(lo, hi) {
    // Get a random number between two values
    return lo + (Math.random() * (hi - lo));
}

function randomVector(scale) {
    // get two values that are tied together
    var x1 = 0;
    var x2 = 0;
    var w = 2.0;
    while (w >= 1) {
        x1 = randBetween(-1, 1);
        x2 = randBetween(-1, 1);
        w = x1 * x1 + x2 * x2;
    }
    w = Math.sqrt(-2 * Math.log(w) / w) * scale;
    // return [scale * x1 * w, scale * x2 * w];
    return scaleByVector([x1, x2], [w, w]);
}

function scaleByVector(vector1, vector2) {
    return [vector1[0] * vector2[0], vector1[1] * vector2[1]];
}

function generatePoints(config) {
    var points = [],
        scaleVector = [
            config.extent.width,
            config.extent.height
        ],
        i;
    points.config = config;
    for (i = 0; i < config.pointCount; i++) {
        points.push(scaleByVector(randomVector(1), scaleVector));
    }
    validatePoints(points);
    return points;
}

function centroid(points) {
    // finds the center of a set of points
    var x = 0,
        y = 0;
    points.forEach(function (point) {
        x += point[0];
        y += point[1];
    });
    return [x / points.length, y / points.length];
}

function validatePoints(points) {
    console.assert(points.length);
    console.assert(points[0], points);
    console.assert(points[0][0]);
}

function validateHeights(h) {
    for (var i = 0; i < h.length; i++) {
        if (h[i] === undefined || h[i] === null || isNaN(h[i])) {
            console.log([i, h[i]]);
            console.assert(h[i]);
            break;
        }
    }
}

function getImprovedPoints(points) {
    validatePoints(points);
    var polygons = voronoi(points, points.extent)
                    .polygons(points);
    var pointsOut = [];
    polygons.forEach(function (points) {
        pointsOut.push(centroid(points));
    });
    pointsOut.config = points.config;
    validatePoints(pointsOut);
    return pointsOut;
}

function generateGoodPoints(config) {
    var points = generatePoints(config);
    validatePoints(points);
    points = points.sort(function (a, b) {
        if (a[0] === b[0]) {
            return a[1] - b[1];
        }
        return a[0] - b[0];
    });
    validatePoints(points);
    return getImprovedPoints(points);
}

function generateGoodMesh(config) {
    var points = generateGoodPoints(config);
    validatePoints(points);
    return makeMesh(points);
}

function makeMesh(points) {
    var vor = voronoi(points);
    var vxs = [];
    var vxids = {};
    var adj = [];
    var edges = [];
    var triangles = [];
    var gunks = [];
    var i;
    for (i = 0; i < vor.edges.length; i++) {
        var edge = vor.edges[i];
        if (edge == undefined) {
            continue;
        }
        edge[0].border = !edge.right || !edge.left;
        edge[1].border = !edge.right || !edge.left;
        edge.border = !edge.right || !edge.left;
        var e0 = vxids[edge[0]];
        var e1 = vxids[edge[1]];
        if (e0 == undefined) {
            e0 = vxs.length;
            vxids[edge[0]] = e0;
            vxs.push(edge[0]);
        }
        if (e1 == undefined) {
            e1 = vxs.length;
            vxids[edge[1]] = e1;
            vxs.push(edge[1]);
        }
        adj[e0] = adj[e0] || [];
        adj[e0].push(e1);
        adj[e1] = adj[e1] || [];
        adj[e1].push(e0);
        edges.push([e0, e1, edge.left, edge.right]);
        triangles[e0] = triangles[e0] || [];
        if (!triangles[e0].includes(edge.left)) {
            triangles[e0].push(edge.left);
        }
        if (edge.right && !triangles[e0].includes(edge.right)) {
            triangles[e0].push(edge.right);
        }
        triangles[e1] = triangles[e1] || [];
        if (!triangles[e1].includes(edge.left)) {
            triangles[e1].push(edge.left);
        }
        if (edge.right && !triangles[e1].includes(edge.right)) {
            triangles[e1].push(edge.right);
        }
    }
    var mesh = {
        points: points,
        vor: vor,
        vxs: vxs,
        adj: adj,
        triangles: triangles,
        edges: edges,
        config: points.config,
        gunks: gunks
    }
    mesh.edgew = mesh.config.extent.width * 0.45;
    mesh.edgeh = mesh.config.extent.height * 0.45;
    mesh.map = function (f) {
        var mapped = vxs.map(f);
        mapped.mesh = mesh;
        return mapped;
    }
    return mesh;
}

function isBorderEdge(mesh, i) {
    return (mesh.adj[i].length < 3);
}

function isNearEdge(mesh, i) {
    return mesh.vxs[i].border;
}

function getNeighbors(mesh, i) {
    var onbs = mesh.adj[i];
    var nbs = [];
    for (var i = 0; i < onbs.length; i++) {
        nbs.push(onbs[i]);
    }
    return nbs;
}

function mesh_distance(mesh, i, j) {
    var p = mesh.vxs[i];
    var q = mesh.vxs[j];
    return distance(p, q);
}

function distance(p, q) {
    return Math.sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
}

function getSlopes(h) {
    // map input points to the amount of slope per point
    var dh = getDownhills(h);
    var slopes = flatten(h.mesh);
    var slopeVectors;
    for (var i = 0; i < h.length; i++) {
        slopeVectors = getSlopeVectors(h, i);
        slopes[i] = Math.sqrt(
            slopeVectors[0] * slopeVectors[0] +
                slopeVectors[1] * slopeVectors[1]);
    }
    return slopes;
}

function getSlopeVectors(heights, i) {
    var nbs = getNeighbors(heights.mesh, i);
    if (nbs.length != 3) {
        return [0, 0];
    }
    var p0 = heights.mesh.vxs[nbs[0]];
    var p1 = heights.mesh.vxs[nbs[1]];
    var p2 = heights.mesh.vxs[nbs[2]];

    var x1 = p1[0] - p0[0];
    var x2 = p2[0] - p0[0];
    var y1 = p1[1] - p0[1];
    var y2 = p2[1] - p0[1];

    var det = x1 * y2 - x2 * y1;
    var d1 = heights[nbs[1]] - heights[nbs[0]];
    var d2 = heights[nbs[2]] - heights[nbs[0]];

    return [
        (y2 * d1 - y1 * d2) / det,
        (-x2 * d1 + x1 * d2) / det
    ];
}

function quantile(points, q) {
    var sortedPoints = [];
    for (var i = 0; i < points.length; i++) {
        sortedPoints[i] = points[i];
    }
    sortedPoints.sort(d3.ascending);
    sortedPoints.config = points.config;
    return d3.quantile(sortedPoints, q);
}

function flatten(mesh) {
    var z = [];
    z.downhill = null;
    for (var i = 0; i < mesh.vxs.length; i++) {
        z[i] = 0;
    }
    z.mesh = mesh;
    z.config = mesh.config;
    return z;
}

function tiltMap(mesh, direction) {
    // Add
    return mesh.map(function (point) {
        return point[0] * direction[0] + point[1] * direction[1];
    });
}

function cone(mesh, slope) {
    return mesh.map(function (point) {
        return Math.pow(point[0] * point[0] + point[1] * point[1], 0.5) * slope;
    });
}

function map(h, f) {
    var newh = h.map(f);
    newh.mesh = h.mesh;
    return newh;
}

function normalize(h) {
    var lo = d3.min(h),
        hi = d3.max(h),
        newh = flatten(h.mesh),
        range = MAX_HEIGHT - MIN_HEIGHT,
        i;
    for (i = 0; i < h.length; i++) {
        newh[i] = ((h[i] - lo) / (hi - lo) * range) + MIN_HEIGHT;
    };
    return newh;
}

function peaky(h) {
    var lo = d3.min(h),
        hi = d3.max(h),
        newh = flatten(h.mesh),
        range = MAX_HEIGHT - MIN_HEIGHT,
        i;
    for (i = 0; i < h.length; i++) {
        newh[i] = Math.sqrt((h[i] - MIN_HEIGHT) / range);
        newh[i] = ((h[i] - lo) / (hi - lo) * range) + MIN_HEIGHT;
    };
    for (i = 0; i < h.length; i++) {

    }
    return newh;
}

function add() {
    var n = arguments[0].length;
    var newvals = flatten(arguments[0].mesh);
    for (var i = 0; i < n; i++) {
        for (var j = 0; j < arguments.length; j++) {
            newvals[i] += arguments[j][i];
        }
    }
    return newvals;
}

function relax(h) {
    var newh = flatten(h.mesh);
    for (var i = 0; i < h.length; i++) {
        var nbs = getNeighbors(h.mesh, i);
        if (nbs.length < 3) {
            newh[i] = 0;
            continue;
        }
        newh[i] = d3.mean(nbs.map(function (j) {return h[j]}));
    }
    return newh;
}


function findSinks(h) {
    var dh = getDownhills(h);
    var sinks = [];
    for (var i = 0; i < dh.length; i++) {
        var node = i;
        while (true) {
            if (isBorderEdge(h.mesh, node)) {
                sinks[i] = -2;
                break;
            }
            if (dh[node] == -1) {
                sinks[i] = node;
                break;
            }
            node = dh[node];
        }
    }
}

function getDownhills(h) {
    if (h.downhill !== null && h.downhill !== undefined) {
        return h.downhill;
    }
    function downfrom(i) {
        var best = i;
        var besth = h[i];
        var nbs = getNeighbors(h.mesh, i);
        for (var j = 0; j < nbs.length; j++) {
            if (h[nbs[j]] < besth) {
                besth = h[nbs[j]];
                best = nbs[j];
            }
        }
        return best;
    }
    var downs = flatten(h.mesh);
    for (var i = 0; i < h.length; i++) {
        downs[i] = downfrom(i);
    }
    h.downhill = downs;
    return downs;
}

function by(elevations, ascending) {
    var inverted = ascending ? 1 : -1;
    return (function (a, b) {
        return inverted * (elevations[a] - elevations[b]);
    });
}

function getWaterDepth(elevations) {
    var startTime = new Date();

    var flux = getFlux(elevations);
    var downhills = getDownhills(elevations);

    var water = flatten(elevations.mesh);
    water.minima = {};
    water.catchments = {};

    var catchments = [];

    var i;
    var indexes = [];
    var index;

    var toCheck = [];

    for (i = 0; i < elevations.length; i++) {
        // find obvious local minima and catchments
        water.minima[i] = null;
        if (elevations[i] < 0) {
            // elevations below 0 are water by default, and count as a minimum
            water[i] = -elevations[i];
            water.minima[i] = i;
            water.catchments[i] = {
                sill: i,
                nodes: [i]
            };
            catchments.push(i);
        } else if (isNearEdge(elevations.mesh, i)) {
            // border nodes drain off the border
            water[i] = 0.1;
            water.minima[i] = i;
            water.catchments[i] = {
                sill: i,
                nodes: [i]
            };
            catchments.push(i);
        } else if (downhills[i] == i) {
            // a local minimum
            water[i] = 0.1;
            water.minima[i] = i;
            water.catchments[i] = {
                sill: null,
                nodes: [i]
            };
            catchments.push(i);
        } else {
            toCheck.push(i);
        }
    }

    var checking;
    var receiver;
    while (toCheck.length) {
        // take the top item from the stack
        checking = toCheck.pop();
        if (water.minima[checking] !== null) {
            // this already has a minimum
            continue;
        }
        receiver = downhills[checking];
        if (water.minima[receiver] !== null) {
            // the receiver for this already has a minimum
            water.minima[checking] = water.minima[receiver];
            continue;
        }
        // we have not found the minimum for this item, so put
        // it back on the stack and check its receiver first.
        toCheck.push(checking);
        toCheck.push(receiver);
    }

    // everything should either be or have a minimum now

    var catchment;
    var neighbors;
    var neighbor;

    for (i = 0; i < elevations.length; i++) {
        indexes[i] = i;
    }
    indexes.sort(by(elevations, true));

    for (i = 0; i < elevations.length; i++) {
        if (elevations[i] > 0 && water[i] <= 0.1) {
            water[i] = 0;
        }
        if (elevations[i] < 0 && water[i] <= 0) {
            water[i] = -elevations[i];
        }
    }

    var endTime = new Date();
    var timeDiff = endTime - startTime;
    timeDiff /= 1000;
    return water;
}


function getFlux(points) {
    var downhills = getDownhills(points);
    var idxs = [];
    var flux = flatten(points.mesh);
    var surfacePointCount = countGround(points);
    var i;
    var j;
    for (i = 0; i < points.length; i++) {
        idxs[i] = i;
        if (points[i] > 0) {
            flux[i] = 1 / surfacePointCount;
        }
    }
    idxs.sort(by(points, false));
    for (i = 0; i < idxs.length; i++) {
        j = idxs[i];
        if (downhills[j] != j) {
            if (points[j] >= 0) {
                flux[downhills[j]] += flux[j];
            }
        }
    }
    return flux;
}

function fillSinks(h) {
    // everything must be epsilon amount higher than at least one of its
    // neighbors, unless it is at the edge
    var newh = flatten(h.mesh),
        i,
        j,
        epsilon = 10,
        indexes = [],
        index = 0,
        epsilon_height,
        changed = true,
        max_height = -INF;;
    for (i = 0; i < h.length; i++) {
        indexes[i] = i;
        if (isNearEdge(h.mesh, i)) {
            newh[i] = h[i];
        } else {
            newh[i] = INF;
        }
    }
    indexes.sort(by(h, false));
    while (changed) {
        changed = false;
        for (index = 0; index < indexes.length; index++) {
            i = indexes[index];
            if (newh[i] == h[i]) {
                continue;
            }
            var nbs = getNeighbors(h.mesh, i);
            for (j = 0; j < nbs.length; j++) {
                epsilon_height = newh[nbs[j]] + epsilon;
                if (h[i] >= epsilon_height) {
                    // if the existing height is greater than any neighbor plus
                    // an epsilon, use that
                    newh[i] = h[i];
                    changed = true;
                    break;
                }
                if ((newh[i] > epsilon_height) && (epsilon_height > h[i])) {
                    // if the new height is greater than the neighbor plus an
                    // epsilon, but the neighbor plus epsilon is already
                    // higher than the original site, use the neighbor plus
                    // epsilon
                    newh[i] = epsilon_height;
                    changed = true;
                }
            }
        }
    }
    for (i = 0; i < newh.length; i++) {
        if (newh[i] > max_height && newh[i] < INF) {
            max_height = newh[i];
        }
    }
    for (i = 0; i < newh.length; i++) {
        if (newh[i] >= INF) {
            // newh[i] = max_height;
        }
    }
    return newh;
}

function getErosionRate(h, slopes) {
    var flux = getFlux(h),
        erosion = flatten(h.mesh),
        max_slope = -INF,
        i;
    for (i = 0; i < h.length; i++) {
        if (h[i] > 0 && slopes[i] > max_slope) {
            max_slope = slopes[i];
        }
    }
    for (i = 0; max_slope > 0 && i < h.length; i++) {
        if (h[i] > 0 && slopes[i]) {
            erosion[i] = slopes[i] / max_slope;
        }
    }
    validateHeights(erosion);
    return erosion;
}

function erode(h, amount) {
    var slopes = getSlopes(h),
        erosionRate = getErosionRate(h, slopes),
        max_r = d3.max(erosionRate) || 1,
        downhills = getDownhills(h),
        max_slope = 0,
        newh = flatten(h.mesh),
        i,
        j,
        k,
        gunk,
        gunkOut,
        lowerNeighbors,
        neighbor,
        neighborSlope,
        idxs = [];

    for (i = 0; i < h.length; i++) {
        idxs[i] = i;
        newh[i] = h[i] - amount * erosionRate[i];
        h.mesh.gunks[i] = h[i] - newh[i];

        if (slopes[i] > max_slope) {
            max_slope = slopes[i];
        }
    }
    validateHeights(h.mesh.gunks);
    idxs.sort(by(h, false));
    for (k = 0; k < idxs.length; k++) {
        i = idxs[k];
        newh[i] += h.mesh.gunks[i];
        h.mesh.gunks[i] = 0;
        lowerNeighbors = getLowerNeighbors(h, i);
        neighborSlope = 0;
        for (j = 0; j < lowerNeighbors.length; j++) {
            neighbor = lowerNeighbors[j];
            neighborSlope += slopes[neighbor];
        }
        gunkOut = h.mesh.gunks[i];
        for (j = 0; j < lowerNeighbors.length; j++) {
            neighbor = lowerNeighbors[j];
            if (neighborSlope === 0 || slopes[neighbor] === 0) {
                gunk = gunkOut / lowerNeighbors.length;
            } else {
                gunk = gunkOut / slopes[neighbor] / neighborSlope;
            }
            h.mesh.gunks[neighbor] += gunk;
            h.mesh.gunks[i] -= gunk;
        }
        if (gunkOut) {
            newh[i] += gunkOut;
        }
    }
    validateHeights(h.mesh.gunks);
    validateHeights(newh);
    return newh;
}

function previewErode(h) {
    var deltas = flatten(h.mesh),
        eroded = erode(h, 0.1),
        i;
    for (i = 0; i < h.length; i++) {
        deltas[i] = eroded[i] - h[i];
    }
    validateHeights(deltas);
    return deltas;
}

function getLowerNeighbors(h, i) {
    var j,
        nbrs = getNeighbors(h.mesh, i),
        lowerNeighbors = [];
    for (j = 0; j < nbrs.length; j++) {
        if (h[nbrs[j]] < h[i]) {
            lowerNeighbors.push(nbrs[j]);
        }
    }
    return lowerNeighbors;
}

function doErosion(h, amount, n) {
    n = n || 1;
    // h = fillSinks(h);
    for (var i = 0; i < n; i++) {
        h = erode(h, amount);
        // h = fillSinks(h);
    }
    return h;
}

function voronoi(points) {
    var w = points.config.extent.width / 2;
    var h = points.config.extent.height / 2;
    validatePoints(points);
    return d3.voronoi().extent([[-w, -h], [w, h]])(points);
}

function bell(center, width) {
    var updown = Math.random() > 0.5 ? 1 : -1;

    return ((Math.random() * Math.random() * width * updown) + center);
}

function roughen(heights, sections, density) {
    var i,
        selected_points = [],
        new_heights = [],
        newheight = 0,
        plates_diagram,
        p,
        cell,
        neighbor,
        neighbor_closeness,
        closeness,
        diagonal = 0,
        scaleVector = [
            heights.mesh.points.config.extent.width,
            heights.mesh.points.config.extent.height
        ];
    diagonal = Math.sqrt((scaleVector[0] * scaleVector[0]) +
                         (scaleVector[1] * scaleVector[1]))

    selected_points.config = heights.mesh.points.config;

    for (i = 0; i < sections; i++) {
        selected_points.push(scaleByVector(randomVector(1), scaleVector));
        if (Math.random() < density) {
            new_heights.push(bell(0, 1000 / density));
        } else {
            new_heights.push(0);
        }
    }
    plates_diagram = voronoi(selected_points);
    for (i = 0; i < heights.length; i++) {
        p = heights.mesh.vxs[i];
        cell = bestCell(plates_diagram, p[0], p[1], 100);
        if (cell && cell.site) {
            closeness = distance(p, cell.site);
            neighbor = neighboring_cell(plates_diagram, p, cell);
            if (neighbor && neighbor.site) {
                neighbor_closeness = distance(p, neighbor.site);
                newheight = ((
                    (new_heights[cell.site.index] * closeness) +
                        (new_heights[neighbor.site.index] * neighbor_closeness)
                ) / (closeness + neighbor_closeness));

            } else {
                newheight = new_heights[cell.site.index] / diagonal * closeness;
            }
        } else {
            newheight = 0;
        }
        heights[i] = heights[i] + newheight;
    }
    return heights;
}

function neighboring_cell(diagram, p, cell) {
    var closest = null,
        d2 = null;
    if (cell) {
        cell.halfedges.forEach(function(e) {
            var edge = diagram.edges[e],
                v = edge.left;
            if ((v === cell.site || !v) && !(v = edge.right)) {
                return;
            }
            var v2 = distance(p, v);
            if (d2 === null || v2 < d2) {
                closest = cell;
                d2 = v2;
            }
        });
    }
    return closest;
}

function bestCell(diagram, x, y, radius) {
    var cell,
        closest = null,
        d2 = (radius * radius) || null,
        dx,
        dy,
        i = diagram._found || 0;

    for (i; i < diagram.cells.length; i++) {
        cell = diagram.cells[i];
        if (cell) {
            var v2 = distance([x, y], cell.site);
            if (d2 === null || v2 < d2) {
                closest = cell;
                d2 = v2;
            }
        }
    }
    return radius == null || d2 <= radius * radius ? closest : null;
}

function setSeaLevel(points, q) {
    // var newh = flatten(points.mesh);
    var delta = quantile(points, q);
    return flood(points, -delta);
}

function flood(points, depth) {
    for (var i = 0; i < points.length; i++) {
        points[i] = points[i] + depth;
    }
    return points;
}

function cleanCoast(h, iters) {
    for (var iter = 0; iter < iters; iter++) {
        var changed = 0;
        var newh = flatten(h.mesh);
        for (var i = 0; i < h.length; i++) {
            newh[i] = h[i];
            var nbs = getNeighbors(h.mesh, i);
            if (h[i] <= 0 || nbs.length != 3) continue;
            var count = 0;
            var best = -INF;
            for (var j = 0; j < nbs.length; j++) {
                if (h[nbs[j]] > 0) {
                    count++;
                } else if (h[nbs[j]] > best) {
                    best = h[nbs[j]];
                }
            }
            if (count > 1) continue;
            newh[i] = best / 2;
            changed++;
        }
        h = newh;
        newh = flatten(h.mesh);
        for (var i = 0; i < h.length; i++) {
            newh[i] = h[i];
            var nbs = getNeighbors(h.mesh, i);
            if (h[i] > 0 || nbs.length != 3) continue;
            var count = 0;
            var best = INF;
            for (var j = 0; j < nbs.length; j++) {
                if (h[nbs[j]] <= 0) {
                    count++;
                } else if (h[nbs[j]] < best) {
                    best = h[nbs[j]];
                }
            }
            if (count > 1) continue;
            newh[i] = best / 2;
            changed++;
        }
        h = newh;
    }
    return h;
}

function cityScore(cities, zero) {
    var score = map(getFlux(zero), Math.sqrt);
    for (var i = 0; i < zero.length; i++) {
        if (zero[i] <= 0 || isNearEdge(zero.mesh, i)) {
            score[i] = -INF;
            continue;
        }
        score[i] += 0.01 / (1e-9 + Math.abs(zero.mesh.vxs[i][0]) - zero.mesh.config.extent.width/2)
        score[i] += 0.01 / (1e-9 + Math.abs(zero.mesh.vxs[i][1]) - zero.mesh.config.extent.height/2)
        for (var j = 0; j < cities.length; j++) {
            score[i] -= 0.02 / (mesh_distance(zero.mesh, cities[j], i) + 1e-9);
        }
    }
    return score;
}
function contour(zero, level) {
    level = level || 0;
    var edges = [];
    for (var i = 0; i < zero.mesh.edges.length; i++) {
        var e = zero.mesh.edges[i];
        if (e[3] == undefined) continue;
        if (isNearEdge(zero.mesh, e[0]) || isNearEdge(zero.mesh, e[1])) {
            continue;
        }
        if ((zero[e[0]] > level && zero[e[1]] <= level) ||
            (zero[e[1]] > level && zero[e[0]] <= level)) {
            edges.push([e[2], e[3]]);
        }
    }
    return mergeSegments(edges);
}

function countGround(points) {
    var above = 0;
    points.forEach(function (elevation) {
        if (elevation >= 0) {
            above++;
        }
    });
    return above;
}

function getRivers(elevations, limit) {
    var downhills = getDownhills(elevations);
    var flux = getFlux(elevations);
    var waters = getWaterDepth(elevations);
    var links = [];
    var above = countGround(elevations);
    limit *= above / elevations.length;
    for (var i = 0; i < downhills.length; i++) {
        if (isNearEdge(elevations.mesh, i)) {
            continue;
        }
        if (flux[i] > limit && elevations[i] > 0 && downhills[i] >= 0) {
            var up = elevations.mesh.vxs[i];
            var down = elevations.mesh.vxs[downhills[i]];
            if (elevations[downhills[i]] > 0) {
                links.push([up, down]);
            } else {
                links.push(
                    [
                        up,
                        [(up[0] + down[0]) / 2, (up[1] + down[1]) / 2]
                    ]
                );
            }
        }
    }
    return mergeSegments(links).map(relaxPath);
}

function placeCity(render, zero) {
    render.cities = render.cities || [];
    var score = cityScore(render.cities, zero);
    var newcity = d3.scan(score, d3.descending);
    render.cities.push(newcity);
}

function placeCities(render, zero) {
    var params = render.params;
    var n = params.ncities;
    for (var i = 0; i < n; i++) {
        placeCity(render, zero);
    }
}

function getTerritories(render, zero) {
    var cities = render.cities;
    var n = render.params.nterrs;
    if (n > render.cities.length) {
        n = render.cities.length;
    }
    var flux = getFlux(zero);
    var terr = [];
    var queue = new PriorityQueue({comparator: function (a, b) {return a.score - b.score}});
    function weight(u, v) {
        var horiz = mesh_distance(zero.mesh, u, v);
        var vert = zero[v] - zero[u];
        if (vert > 0) vert /= 10;
        var diff = 1 + 0.25 * Math.pow(vert/horiz, 2);
        diff += 100 * Math.sqrt(flux[u]);
        if (zero[u] <= 0) {
            diff = 100;
        }
        if ((zero[u] > 0) != (zero[v] > 0)) {
            return 1000;
        }
        return horiz * diff;
    }
    for (var i = 0; i < n; i++) {
        terr[cities[i]] = cities[i];
        var nbs = getNeighbors(zero.mesh, cities[i]);
        for (var j = 0; j < nbs.length; j++) {
            queue.queue({
                score: weight(cities[i], nbs[j]),
                city: cities[i],
                vx: nbs[j]
            });
        }
    }
    while (queue.length) {
        var u = queue.dequeue();
        if (terr[u.vx] != undefined) {
            continue;
        }
        terr[u.vx] = u.city;
        var nbs = getNeighbors(zero.mesh, u.vx);
        for (var i = 0; i < nbs.length; i++) {
            var v = nbs[i];
            if (terr[v] != undefined) {
                continue;
            }
            var newdist = weight(u.vx, v);
            queue.queue({
                score: u.score + newdist,
                city: u.city,
                vx: v
            });
        }
    }
    terr.mesh = zero.mesh;
    return terr;
}

function getBorders(render, zero) {
    var terr = render.terr;
    var edges = [];
    for (var i = 0; i < terr.mesh.edges.length; i++) {
        var e = terr.mesh.edges[i];
        if (e[3] == undefined) {
            continue;
        }
        if (isNearEdge(terr.mesh, e[0]) || isNearEdge(terr.mesh, e[1])) {
            continue;
        }
        if (zero[e[0]] < 0 || zero[e[1]] < 0) {
            continue;
        }
        if (terr[e[0]] != terr[e[1]]) {
            edges.push([e[2], e[3]]);
        }
    }
    return mergeSegments(edges).map(relaxPath);
}

function mergeSegments(segs) {
    var adj = {};
    for (var i = 0; i < segs.length; i++) {
        var seg = segs[i];
        var a0 = adj[seg[0]] || [];
        var a1 = adj[seg[1]] || [];
        a0.push(seg[1]);
        a1.push(seg[0]);
        adj[seg[0]] = a0;
        adj[seg[1]] = a1;
    }
    var done = [];
    var paths = [];
    var path = null;
    while (true) {
        if (path == null) {
            for (var i = 0; i < segs.length; i++) {
                if (done[i]) continue;
                done[i] = true;
                path = [segs[i][0], segs[i][1]];
                break;
            }
            if (path == null) break;
        }
        var changed = false;
        for (var i = 0; i < segs.length; i++) {
            if (done[i]) continue;
            if (adj[path[0]].length == 2 && segs[i][0] == path[0]) {
                path.unshift(segs[i][1]);
            } else if (adj[path[0]].length == 2 && segs[i][1] == path[0]) {
                path.unshift(segs[i][0]);
            } else if (adj[path[path.length - 1]].length == 2 && segs[i][0] == path[path.length - 1]) {
                path.push(segs[i][1]);
            } else if (adj[path[path.length - 1]].length == 2 && segs[i][1] == path[path.length - 1]) {
                path.push(segs[i][0]);
            } else {
                continue;
            }
            done[i] = true;
            changed = true;
            break;
        }
        if (!changed) {
            paths.push(path);
            path = null;
        }
    }
    return paths;
}

function relaxPath(path) {
    var newpath = [path[0]];
    for (var i = 1; i < path.length - 1; i++) {
        var newpt = [0.25 * path[i-1][0] + 0.5 * path[i][0] + 0.25 * path[i+1][0],
                     0.25 * path[i-1][1] + 0.5 * path[i][1] + 0.25 * path[i+1][1]];
        newpath.push(newpt);
    }
    newpath.push(path[path.length - 1]);
    return newpath;
}

function makeD3Path(path) {
    var p = d3.path();
    p.moveTo(1000 * path[0][0], 1000 * path[0][1]);
    for (var i = 1; i < path.length; i++) {
        p.lineTo(1000 * path[i][0], 1000 * path[i][1]);
    }
    return p.toString();
}

function dropEdge(h, p) {
    p = p || 4
    var newh = flatten(h.mesh);
    for (var i = 0; i < h.length; i++) {
        var v = h.mesh.vxs[i];
        var x = 2.4 * v[0] / h.mesh.config.extent.width;
        var y = 2.4 * v[1] / h.mesh.config.extent.height;
        newh[i] = h[i] - Math.exp(10*(Math.pow(Math.pow(x, p) + Math.pow(y, p), 1/p) - 1));
    }
    return newh;
}

function mountains(mesh, number, radius, height) {
    radius = radius || 0.05;
    height = height || 1;
    var mounts = [];
    var i;
    for (i = 0; i < number; i++) {
        mounts.push([
            mesh.config.extent.width * (Math.random() - 0.5),
            mesh.config.extent.height * (Math.random() - 0.5)
        ]);
    }
    var newvals = flatten(mesh);
    for (i = 0; i < mesh.vxs.length; i++) {
        var p = mesh.vxs[i];
        for (var j = 0; j < number; j++) {
            var m = mounts[j];
            newvals[i] += (Math.pow(
                Math.exp(
                    -(
                        (p[0] - m[0]) *
                        (p[0] - m[0]) +
                        (p[1] - m[1]) *
                        (p[1] - m[1])
                    ) / (2 * radius * radius)
                ),
                2
            ) * height);
        }
    }
    return newvals;
}

function generateRandomHeightmap(zero) {
    var h = add(
        tiltMap(zero.mesh, randomVector(4)),
        cone(zero.mesh, randBetween(-1, -1)),
        mountains(zero.mesh, 50)
    );
    for (var i = 0; i < 3; i++) {
        h = relax(h);
    }
    h = roughen(h, 100, 0.25);
    validateHeights(h);
    h = peaky(h);
    validateHeights(h);
    h = add(h, mountains(zero.mesh, 4));
    validateHeights(h);
    h = relax(h);
    validateHeights(h);
    h = doErosion(h, randBetween(0, 0.1), 5);
    validateHeights(h);
    h = add(h, mountains(zero.mesh, 3));
    validateHeights(h);
    h = relax(h);
    validateHeights(h);
    h = roughen(h, 100, 0.25);
    validateHeights(h);
    h = doErosion(h, randBetween(0, 0.1), 5);
    h = setSeaLevel(h, randBetween(0.2, 0.6));
    h = add(h, mountains(zero.mesh, 2));
    h = add(h, mountains(zero.mesh, 2, 0.02, -0.5));
    h = relax(h);
    // h = fillSinks(h);
    h = add(h, mountains(zero.mesh, 1));
    h = cleanCoast(h, 3);
    return h;
}

function terrCenter(h, terr, city, landOnly) {
    var x = 0;
    var y = 0;
    var n = 0;
    for (var i = 0; i < terr.length; i++) {
        if (terr[i] != city) {
            continue;
        }
        if (landOnly && h[i] <= 0) {
            continue;
        }
        x += terr.mesh.vxs[i][0];
        y += terr.mesh.vxs[i][1];
        n++;
    }
    return [x/n, y/n];
}
