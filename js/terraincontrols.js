var points = 800;

var defaults = {
    npts: points,
    extent: defaultExtent
};

var primDiv = d3.select("#prim");

var pointsControls = d3.select("#pointsControls");
var heightControls = d3.select("#heightControls");
var viewControls = d3.select("#viewControls");

var primSVG = addSVG(primDiv);
var primZero = zero(generateGoodMesh(points));

var viewMode = "default";

var primViewCoast = true;
var primViewRivers = true;
var primViewSlope = true;

var cityViewScore = false;

var cityRender = new CityRender(primZero);

function addSVG(div) {
    return div.insert("svg", ":first-child")
           .attr("height", 800)
           .attr("width", 800)
           .attr("viewBox", "-500 -500 1000 1000");
}

function CityRender(h) {
    this.params = defaultParams;
    this.cities = [];
    this.terr = null;
    this.clearCities = function () {
        this.cities = [];
    };
}

function primDraw() {
    primSVG.selectAll().remove();
    cityRender.terr = getTerritories(cityRender, primZero);
    if (viewMode == "erodeViewErosion") {
        visualizeVoronoi(primSVG, erosionRate(primZero));
    } else if (viewMode == "cityViewScore") {
        var score = cityScore(cityRender.cities, primZero);
        visualizeVoronoi(primSVG, score, d3.max(score) - 0.5);
    } else if (viewMode == "territories") {
        visualizeVoronoi(primSVG, cityRender.terr);
    } else if (viewMode == "erodeViewRate") {
        visualizeVoronoi(primSVG, primZero, 0, 1);
    } else if (viewMode == "default") {
        visualizeVoronoi(primSVG, primZero, -1, 1);
    } else if (viewMode == "nothing") {
        primSVG.selectAll("path.field").remove();
    }
    if (primViewRivers) {
        drawPaths(primSVG, "river", getRivers(primZero, 0.01));
    } else {
        drawPaths(primSVG, "river", []);
    }
    if (primViewSlope) {
        visualizeSlopes(primSVG, primZero);
    } else {
        visualizeSlopes(primSVG, zero(primZero.mesh));
    }

    drawPaths(primSVG, 'border', getBorders(cityRender, primZero));
    visualizeCities(primSVG, cityRender, primZero);

    if (primViewCoast) {
        drawPaths(primSVG, "coast", contour(primZero, 0));
    } else {
        drawPaths(primSVG, "coast", []);
    }


}

primDraw();

pointsControls.append("button")
    .text("Generate random points")
    .on("click", function () {
        primZero = zero(generateGoodMesh(points));
        cityRender.clearCities();
        primDraw();
    });

pointsControls.append("button")
    .text("Improve points")
    .on("click", function () {
        var pts = improvePoints(
            primZero.mesh.pts, 1, primZero.mesh.extent);
        primZero = zero(makeMesh(pts, primZero.mesh.extent));
        cityRender.clearCities();
        primDraw();
    });


heightControls.append("button")
    .text("Generate random heightmap")
    .on("click", function () {
        primZero = generateRandomHeightmap(primZero);
        cityRender.clearCities();
        primDraw();
    });

heightControls.append("button")
    .text("Reset to flat")
    .on("click", function () {
        primZero = zero(primZero.mesh);
        cityRender.clearCities();
        primDraw();
    });

heightControls.append("button")
    .text("Add random slope")
    .on("click", function () {
        primZero = add(primZero, slope(primZero.mesh, randomVector(4)));
        cityRender.clearCities();
        primDraw();
    });

heightControls.append("button")
    .text("Add cone")
    .on("click", function () {
        primZero = add(primZero, cone(primZero.mesh, -0.5));
        primDraw();
    });

heightControls.append("button")
    .text("Add inverted cone")
    .on("click", function () {
        primZero = add(primZero, cone(primZero.mesh, 0.5));
        primDraw();
    });

heightControls.append("button")
    .text("Add five blobs")
    .on("click", function () {
        primZero = add(primZero, mountains(primZero.mesh, 5));
        primDraw();
    });

heightControls.append("button")
    .text("Normalize heightmap")
    .on("click", function () {
        primZero = normalize(primZero);
        primDraw();
    });

heightControls.append("button")
    .text("Round hills")
    .on("click", function () {
        primZero = peaky(primZero);
        primDraw();
    });

heightControls.append("button")
    .text("Relax")
    .on("click", function () {
        primZero = relax(primZero);
        primDraw();
    });

heightControls.append("button")
    .text("Set sea level to median")
    .on("click", function () {
        primZero = setSeaLevel(primZero, 0.5);
        primDraw();
    });

heightControls.append("button")
    .text("Erode")
    .on("click", function () {
        primZero = doErosion(primZero, 0.1);
        primDraw();
    });

heightControls.append("button")
    .text("Clean coastlines")
    .on("click", function () {
        primZero = cleanCoast(primZero, 1);
        primZero = fillSinks(primZero);
        primDraw();
    });

var modes = {
    "default": "Showing Heightmap",
    "erodeViewErosion": "Showing Erosion",
    "erodeViewRate": "Showing Erosion Rate",
    "cityViewScore": "Showing City Scoring",
    "territories": "Showing Territories",
    "nothing": "No shading."
}

var erodeBut = primDiv.append("button")
    .text(modes[viewMode])
    .on("click", function () {
        if (modes[viewMode] == undefined) {
            viewMode = "default";
        }
        var nextMode = "default";
        var useNext = false;
        for (var mode in modes) {
            if (useNext) {
                nextMode = mode;
                break;
            }
            if (mode == viewMode) {
                useNext = true;
            }
        }
        viewMode = nextMode;
        erodeBut.text(modes[viewMode]);
        primDraw();
    });

var primCoastBut = primDiv.append("button")
    .text("Showing coastline")
    .on("click", function () {
        primViewCoast = !primViewCoast;
        primCoastBut.text(primViewCoast ? "Showing coastline" : "Hiding coastline");
        primDraw();
    });

var primRiverBut = primDiv.append("button")
    .text("Showing rivers")
    .on("click", function () {
        primViewRivers = !primViewRivers;
        primRiverBut.text(primViewRivers ? "Showing rivers" : "Hiding rivers");
        primDraw();
    });


var primSlopeBut = primDiv.append("button")
    .text("Showing slope shading")
    .on("click", function () {
        primViewSlope = !primViewSlope;
        primSlopeBut.text(primViewSlope ? "Showing slope shading" : "Hiding slope shading");
        primDraw();
    });

primDiv.append("button")
    .text("Add new city")
    .on("click", function () {
        placeCity(cityRender, primZero);
        primDraw();
    });

primDiv.append("button")
    .text("Clear cities")
    .on("click", function () {
        cityRender.clearCities();
        primDraw();
    });
