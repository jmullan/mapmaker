function visualizePoints(svg, points) {
    var circle = svg.selectAll('circle').data(points);
    circle.enter()
        .append('circle');
    circle.exit().remove();
    d3.selectAll('circle')
        .attr('cx', function (d) {return 1000 * d[0]})
        .attr('cy', function (d) {return 1000 * d[1]})
        .attr('r', 100 / Math.sqrt(points.length));
}


function visualizeField(className, svg, field, lo, hi) {
    if (hi == undefined) {
        hi = d3.max(field) + 1e-9;
    }
    if (lo == undefined) {
        lo = d3.min(field) - 1e-9;
    }
    var mappedvals = field.map(
        function (x) {
            if (x > hi) {
                return 1;
            } else if (x < lo) {
                return 0;
            } else {
                return (x - lo) / (hi - lo)
            }
        });
    var tris = svg.selectAll('path.' + className).data(field.mesh.triangles)
    tris.enter()
        .append('path')
        .classed(className, true);

    tris.exit()
        .remove();

    svg.selectAll('path.' + className)
        .attr('d', makeD3Path)
        .style('fill', function (d, i) {
            if (mappedvals[i] == 0) {
                return 'transparent';
            }
            return d3.interpolateViridis(mappedvals[i]);
        })
        .on("click",
            function (datum, index) {
                console.log([index, field[index]]);
            });
}

function visualizeWater(svg, field, lo, hi) {
    visualizeField('water', svg, field, lo, hi);
}

function visualizeDownhills(h) {
    var links = getRivers(h, 0.01);
    drawPaths('river', links);
}

function drawPaths(svg, cls, paths) {
    var paths = svg.selectAll('path.' + cls).data(paths)
    paths.enter()
         .append('path')
         .classed(cls, true)
    paths.exit()
         .remove();
    svg.selectAll('path.' + cls)
         .attr('d', makeD3Path);
}

function visualizeSlopes(svg, zero) {
    var strokes = [];
    var r = 0.25 / Math.sqrt(zero.length);
    for (var i = 0; i < zero.length; i++) {
        if (zero[i] <= 0 || isNearEdge(zero.mesh, i)) {
            continue;
        }
        var nbs = getNeighbors(zero.mesh, i);
        nbs.push(i);
        var s = 0;
        var s2 = 0;
        for (var j = 0; j < nbs.length; j++) {
            var slopes = getSlopeVectors(zero, nbs[j]);
            s += slopes[0] / 10;
            s2 += slopes[1];
        }
        s /= nbs.length;
        s2 /= nbs.length;
        if (Math.abs(s) < randBetween(0.1, 0.4)) {
            continue;
        }
        var l = (
            r *
                randBetween(1, 2) *
                (1 - 0.2 * Math.pow(Math.atan(s), 2)) *
                Math.exp(s2 / 100)
        );
        var x = zero.mesh.vxs[i][0];
        var y = zero.mesh.vxs[i][1];
        if (Math.abs(l*s) > 2 * r) {
            var n = Math.floor(Math.abs(l * s / r));
            l /= n;
            if (n > 4) n = 4;
            for (var j = 0; j < n; j++) {
                var vector = randomVector(1);
                var u = vector[0] * r;
                var v = vector[1] * r;
                strokes.push(
                    [
                        [x + u - l, y + v + l * s],
                        [x + u + l, y + v - l * s]
                    ]
                );
            }
        } else {
            strokes.push([[x - l, y + l * s], [x + l, y - l * s]]);
        }
    }
    var lines = svg.selectAll('line.slope').data(strokes)
    lines.enter()
            .append('line')
            .classed('slope', true);
    lines.exit()
            .remove();
    svg.selectAll('line.slope')
        .attr('x1', function (d) {return 1000 * d[0][0]})
        .attr('y1', function (d) {return 1000 * d[0][1]})
        .attr('x2', function (d) {return 1000 * d[1][0]})
        .attr('y2', function (d) {return 1000 * d[1][1]})
}


function visualizeContour(h, level) {
    level = level || 0;
    var links = contour(h, level);
    drawPaths('coast', links);
}

function visualizeBorders(zero, cities) {
    var links = getBorders(getTerritories(zero, cities), zero);
    drawPaths('border', links);
}


function visualizeCities(svg, render, zero) {
    var cities = render.cities;
    var n = render.params.nterrs;
    var circs = svg.selectAll('circle.city').data(cities);
    circs.enter()
            .append('circle')
            .classed('city', true);
    circs.exit()
            .remove();
    svg.selectAll('circle.city')
        .attr('cx', function (d) {return 1000 * zero.mesh.vxs[d][0]})
        .attr('cy', function (d) {return 1000 * zero.mesh.vxs[d][1]})
        .attr('r', function (d, i) {return i >= n ? 4 : 10})
        .style('fill', 'white')
        .style('stroke-width', 5)
        .style('stroke-linecap', 'round')
        .style('stroke', 'black')
        .raise();
}

function visualizeVoronoi(svg, field, lo, hi) {
    visualizeField('field', svg, field, lo, hi);
}


function drawMap(svg, render, zero) {
    render.rivers = getRivers(zero, 0.01);
    render.coasts = contour(zero, 0);
    render.terr = getTerritories(render, zero);
    render.borders = getBorders(render, zero);
    drawPaths(svg, 'river', render.rivers);
    drawPaths(svg, 'coast', render.coasts);
    drawPaths(svg, 'border', render.borders);
    visualizeSlopes(svg, render);
    visualizeCities(svg, render);
    drawLabels(svg, render);
}

function drawLabels(svg, render, zero) {
    render.rivers = getRivers(zero, 0.01);
    render.coasts = contour(zero, 0);
    render.terr = getTerritories(render, zero);
    render.borders = getBorders(render, zero);

    var params = render.params;
    var terr = render.terr;
    var cities = render.cities;
    var nterrs = render.params.nterrs;
    var avoids = [render.rivers, render.coasts, render.borders];
    var lang = makeRandomLanguage();
    var citylabels = [];
    function penalty(label) {
        var pen = 0;
        if (label.x0 < -0.45 * zero.mesh.config.extent.width) {
            pen += 100;
        }
        if (label.x1 > 0.45 * zero.mesh.config.extent.width) {
            pen += 100;
        }
        if (label.y0 < -0.45 * zero.mesh.config.extent.height) {
            pen += 100;
        }
        if (label.y1 > 0.45 * zero.mesh.config.extent.height) {
            pen += 100;
        }
        for (var i = 0; i < citylabels.length; i++) {
            var olabel = citylabels[i];
            if (label.x0 < olabel.x1 && label.x1 > olabel.x0 &&
                label.y0 < olabel.y1 && label.y1 > olabel.y0) {
                pen += 100;
            }
        }

        for (var i = 0; i < cities.length; i++) {
            var c = zero.mesh.vxs[cities[i]];
            if (label.x0 < c[0] && label.x1 > c[0] && label.y0 < c[1] && label.y1 > c[1]) {
                pen += 100;
            }
        }
        for (var i = 0; i < avoids.length; i++) {
            var avoid = avoids[i];
            for (var j = 0; j < avoid.length; j++) {
                var pathToAvoid = avoid[j];
                for (var k = 0; k < pathToAvoid.length; k++) {
                    var pt = pathToAvoid[k];
                    if (pt[0] > label.x0 && pt[0] < label.x1 && pt[1] > label.y0 && pt[1] < label.y1) {
                        pen++;
                    }
                }
            }
        }
        return pen;
    }
    for (var i = 0; i < cities.length; i++) {
        var x = zero.mesh.vxs[cities[i]][0];
        var y = zero.mesh.vxs[cities[i]][1];
        var text = makeName(lang, 'city');
        var size = i < nterrs ? params.fontsizes.city : params.fontsizes.town;
        var sx = 0.65 * size/1000 * text.length;
        var sy = size/1000;
        var posslabels = [
        {
            x: x + 0.8 * sy,
            y: y + 0.3 * sy,
            align: 'start',
            x0: x + 0.7 * sy,
            y0: y - 0.6 * sy,
            x1: x + 0.7 * sy + sx,
            y1: y + 0.6 * sy
        },
        {
            x: x - 0.8 * sy,
            y: y + 0.3 * sy,
            align: 'end',
            x0: x - 0.9 * sy - sx,
            y0: y - 0.7 * sy,
            x1: x - 0.9 * sy,
            y1: y + 0.7 * sy
        },
        {
            x: x,
            y: y - 0.8 * sy,
            align: 'middle',
            x0: x - sx/2,
            y0: y - 1.9*sy,
            x1: x + sx/2,
            y1: y - 0.7 * sy
        },
        {
            x: x,
            y: y + 1.2 * sy,
            align: 'middle',
            x0: x - sx/2,
            y0: y + 0.1*sy,
            x1: x + sx/2,
            y1: y + 1.3*sy
        }
        ];
        var label = posslabels[d3.scan(posslabels, function (a, b) {return penalty(a) - penalty(b)})];
        label.text = text;
        label.size = size;
        citylabels.push(label);
    }
    var texts = svg.selectAll('text.city').data(citylabels);
    texts.enter()
        .append('text')
        .classed('city', true);
    texts.exit()
        .remove();
    svg.selectAll('text.city')
        .attr('x', function (d) {return 1000*d.x})
        .attr('y', function (d) {return 1000*d.y})
        .style('font-size', function (d) {return d.size})
        .style('text-anchor', function (d) {return d.align})
        .text(function (d) {return d.text})
        .raise();

    var reglabels = [];
    for (var i = 0; i < nterrs; i++) {
        var city = cities[i];
        var text = makeName(lang, 'region');
        var sy = params.fontsizes.region / 1000;
        var sx = 0.6 * text.length * sy;
        var lc = terrCenter(zero, terr, city, true);
        var oc = terrCenter(zero, terr, city, false);
        var best = 0;
        var bestscore = -999999;
        for (var j = 0; j < zero.length; j++) {
            var score = 0;
            var v = zero.mesh.vxs[j];
            score -= 3000 * Math.sqrt((v[0] - lc[0]) * (v[0] - lc[0]) + (v[1] - lc[1]) * (v[1] - lc[1]));
            score -= 1000 * Math.sqrt((v[0] - oc[0]) * (v[0] - oc[0]) + (v[1] - oc[1]) * (v[1] - oc[1]));
            if (terr[j] != city) score -= 3000;
            for (var k = 0; k < cities.length; k++) {
                var u = zero.mesh.vxs[cities[k]];
                if (Math.abs(v[0] - u[0]) < sx &&
                    Math.abs(v[1] - sy/2 - u[1]) < sy) {
                    score -= k < nterrs ? 4000 : 500;
                }
                if (v[0] - sx/2 < citylabels[k].x1 &&
                    v[0] + sx/2 > citylabels[k].x0 &&
                    v[1] - sy < citylabels[k].y1 &&
                    v[1] > citylabels[k].y0) {
                    score -= 5000;
                }
            }
            for (var k = 0; k < reglabels.length; k++) {
                var label = reglabels[k];
                if (v[0] - sx/2 < label.x + label.width/2 &&
                    v[0] + sx/2 > label.x - label.width/2 &&
                    v[1] - sy < label.y &&
                    v[1] > label.y - label.size) {
                    score -= 20000;
                }
            }
            if (zero[j] <= 0) {
                score -= 500;
            }
            if (v[0] + sx/2 > 0.5 * zero.mesh.config.extent.width) {
                score -= 50000;
            }
            if (v[0] - sx/2 < -0.5 * zero.mesh.config.extent.width) {
                score -= 50000;
            }
            if (v[1] > 0.5 * zero.mesh.config.extent.height) {
                score -= 50000;
            }
            if (v[1] - sy < -0.5 * zero.mesh.config.extent.height) {
                score -= 50000;
            }
            if (score > bestscore) {
                bestscore = score;
                best = j;
            }
        }
        reglabels.push({
            text: text,
            x: zero.mesh.vxs[best][0],
            y: zero.mesh.vxs[best][1],
            size:sy,
            width:sx
        });
    }
    texts = svg.selectAll('text.region').data(reglabels);
    texts.enter()
        .append('text')
        .classed('region', true);
    texts.exit()
        .remove();
    svg.selectAll('text.region')
        .attr('x', function (d) {return 1000*d.x})
        .attr('y', function (d) {return 1000*d.y})
        .style('font-size', function (d) {return 1000*d.size})
        .style('text-anchor', 'middle')
        .text(function (d) {return d.text})
        .raise();

}
