var defaults = {
    pointCount: 4096,
    extent: {
        width: 1,
        height: 1
    },
    ncities: 15,
    nterrs: 5,
    fontsizes: {
        region: 40,
        city: 25,
        town: 20
    }
}

var primDiv = d3.select("#prim");

var pointsControls = d3.select("#pointsControls");
var heightControls = d3.select("#heightControls");
var viewControls = d3.select("#viewControls");

var primSVG = addSVG(primDiv);
var primZero = flatten(generateGoodMesh(defaults.pointCount, defaults.extent));

var viewMode = "default";

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

function addCheckbox(toElement, id, label, initial, onChange) {
    var holder = toElement.append("div");
    holder.append("input")
        .attr("type", "checkbox")
        .attr("checked", initial)
        .attr("id", id)
        .on("change", onChange);
    holder.append("label").attr("for", id).text(label);
}

function CityRender(h) {
    this.params = defaults;
    this.cities = [];
    this.terr = [];
    this.rivers = [];
    this.coasts = [];
    this.borders = [];
    this.clearCities = function () {
        this.cities = [];
    };
}

function primDraw() {
    primSVG.selectAll().remove();
    cityRender.terr = getTerritories(cityRender, primZero);

    var viewMode = document.getElementById("shadingMode").value;
    if (viewMode == "erodeViewErosion") {
        visualizeVoronoi(primSVG, primZero, 0, 1);
    } else if (viewMode == "downhill") {
        visualizeVoronoi(primSVG, downhill(primZero));
    } else if (viewMode == "flux") {
        visualizeVoronoi(primSVG, getFlux(primZero));
    } else if (viewMode == "slope") {
        visualizeVoronoi(primSVG, getSlope(primZero));
    } else if (viewMode == "cityViewScore") {
        var score = cityScore(cityRender.cities, primZero);
        visualizeVoronoi(primSVG, score, d3.max(score) - 0.5);
    } else if (viewMode == "territories") {
        visualizeVoronoi(primSVG, cityRender.terr);
    } else if (viewMode == "erodeViewRate") {
        visualizeVoronoi(primSVG, erosionRate(primZero));
    } else if (viewMode == "heightmap") {
        visualizeVoronoi(primSVG, primZero, -1, 1);
    } else if (viewMode == "nothing") {
        primSVG.selectAll("path.field").remove();
    }
    if (d3.select("#showRivers").property("checked")) {
        drawPaths(primSVG, "river", getRivers(primZero, 0.01));
    } else {
        drawPaths(primSVG, "river", []);
    }
    if (d3.select("#showSlopes").property("checked")) {
        visualizeSlopes(primSVG, primZero);
    } else {
        visualizeSlopes(primSVG, flatten(primZero.mesh));
    }

    drawPaths(primSVG, 'border', getBorders(cityRender, primZero));
    visualizeCities(primSVG, cityRender, primZero);

    if (d3.select("#showCoast").property("checked")) {
        drawPaths(primSVG, "coast", contour(primZero, 0));
    } else {
        drawPaths(primSVG, "coast", []);
    }

    drawLabels(primSVG, cityRender, primZero);
}

pointsControls.append("button")
    .text("Generate random points")
    .on("click", function () {
        primZero = flatten(generateGoodMesh(points));
        cityRender.clearCities();
        primDraw();
    });

pointsControls.append("button")
    .text("Improve points")
    .on("click", function () {
        var pts = improvePoints(
            primZero.mesh.pts, 1, primZero.mesh.extent);
        primZero = flatten(makeMesh(pts, primZero.mesh.extent));
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
        primZero = flatten(primZero.mesh);
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
    "nothing": "None",
    "heightmap": "Heightmap",
    "erodeViewErosion": "Erosion",
    "erodeViewRate": "Erosion Rate",
    "cityViewScore": "City Scoring",
    "territories": "Territories",
    "flux": "Flux",
    "slope": "Slope",
    "downhill": "Downhill"
}


function addSelectBox(toElement, id, label, defaultValue, options, onChange) {
    var holder = toElement.append("div");
    var select = holder.append("select")
        .attr("id", id)
        .on("change", onChange);
    var size = 0;
    for (var optionValue in options) {
        if (options.hasOwnProperty(optionValue)) {
            size += 1;
            var option = select.append("option").attr("value", optionValue);
            option.text(options[optionValue]);
            if (optionValue === defaultValue) {
                option.property("selected", "selected");
            }
        }
    }
    select.attr("size", size);
    holder.append("label").attr("for", id).text(label);

}
addSelectBox(
    viewControls, "shadingMode", "Shading", "nothing", modes, primDraw);

addCheckbox(viewControls, "showCoast", "Coastlines", true, primDraw);
addCheckbox(viewControls, "showRivers", "Rivers", true, primDraw);
addCheckbox(viewControls, "showSlopes", "Slope shading", true, primDraw);

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

primDraw();
