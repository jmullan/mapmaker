var defaults = {
    pointCount: 640,
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
var primZero = flatten(generateGoodMesh(defaults));

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
    if (viewMode == "flux") {
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
        visualizeVoronoi(primSVG, primZero);
    } else if (viewMode == "waterDepth") {
        visualizeVoronoi(primSVG, getWaterDepth(primZero), 0, 1, ["#a0d6b4", "#4169e1"]);
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
    /*
    if (d3.select("#showWater").property("checked")) {
        visualizeWater(primSVG, getWaterDepth(primZero), -1, 1);
    } else {
        primSVG.selectAll("path.water").remove();
    }
     */
    if (d3.select("#showCoast").property("checked")) {
        drawPaths(primSVG, "coast", []);
        drawPaths(primSVG, "coast", contour(primZero, 0));
    } else {
        drawPaths(primSVG, "coast", []);
    }

    drawLabels(primSVG, cityRender, primZero);
}




pointsControls.append("input")
    .attr("type", "range")
    .attr("id", "pointsCount")
    .attr("value", defaults.pointCount)
    .attr("min", 100)
    .attr("max", 64000)
    .on("change", function () {
        defaults.pointCount = this.value;
        primZero = flatten(generateGoodMesh(defaults));
        cityRender = new CityRender(primZero)
        primDraw();
    });

pointsControls.append("button")
    .text("Generate random points")
    .on("click", function () {
        primZero = flatten(generateGoodMesh(defaults));
        cityRender = new CityRender(primZero)
        primDraw();
    });

pointsControls.append("button")
    .text("Improve points")
    .on("click", function () {
        var mesh = primZero.mesh;
        var points = getImprovedPoints(mesh.points);
        primZero = flatten(makeMesh(points));
        cityRender = new CityRender(primZero)
        primDraw();
    });


heightControls.append("button")
    .text("Generate random heightmap")
    .on("click", function () {
        primZero = generateRandomHeightmap(primZero);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Reset to flat")
    .on("click", function () {
        primZero = flatten(primZero.mesh);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Add random slope")
    .on("click", function () {
        primZero = add(primZero, tiltMap(primZero.mesh, randomVector(4)));
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Add cone")
    .on("click", function () {
        primZero = add(primZero, cone(primZero.mesh, -0.5));
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Add inverted cone")
    .on("click", function () {
        primZero = add(primZero, cone(primZero.mesh, 0.5));
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Add blob")
    .on("click", function () {
        primZero = add(primZero, mountains(
            primZero.mesh, 1, randBetween(0.02, 0.10), randBetween(0.01, 1)));
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Add scoop")
    .on("click", function () {
        primZero = add(primZero, mountains(
            primZero.mesh, 1, 0.05, -randBetween(0.01, 1)));
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Normalize heightmap")
    .on("click", function () {
        primZero = normalize(primZero);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Round hills")
    .on("click", function () {
        primZero = peaky(primZero);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Relax")
    .on("click", function () {
        primZero = relax(primZero);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Set sea level to median")
    .on("click", function () {
        primZero = setSeaLevel(primZero, 0.5);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Erode")
    .on("click", function () {
        primZero = doErosion(primZero, 0.1);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Roughen")
    .on("click", function () {
        primZero = roughen(primZero, 20, 0.1);
        primZero = roughen(primZero, 100, 0.2);
        primZero = roughen(primZero, 1000, 0.4);
        primZero = doErosion(primZero, 0.1);
        cityRender = new CityRender(primZero)
        primDraw();
    });

heightControls.append("button")
    .text("Clean coastlines")
    .on("click", function () {
        primZero = cleanCoast(primZero, 1);
        primZero = fillSinks(primZero);
        cityRender = new CityRender(primZero)
        primDraw();
    });

var modes = {
    "nothing": "None",
    "heightmap": "Heightmap",
    "erodeViewRate": "Erosion Rate",
    "cityViewScore": "City Scoring",
    "territories": "Territories",
    "flux": "Flux",
    "slope": "Slope",
    "waterDepth": "Water"
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
    viewControls, "shadingMode", "Shading", "waterDepth", modes, primDraw);

addCheckbox(viewControls, "showCoast", "Coastlines", true, primDraw);
addCheckbox(viewControls, "showRivers", "Rivers", true, primDraw);
addCheckbox(viewControls, "showSlopes", "Slope shading", true, primDraw);
addCheckbox(viewControls, "showWater", "Water", false, primDraw);

primDiv.append("button")
    .text("Add new city")
    .on("click", function () {
        placeCity(cityRender, primZero);
        primDraw();
    });

primDiv.append("button")
    .text("Clear cities")
    .on("click", function () {
        cityRender = new CityRender(primZero)
        primDraw();
    });

primDraw();
