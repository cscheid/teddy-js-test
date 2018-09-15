/*global Delaunator */

var pts = [];
var downPos;
var minDistance = 30;
var isDown = false;
d3.select("#main")
    .append("rect")
    .attr("width", 600)
    .attr("height", 600)
    .attr("cursor", "pointer")
    .attr("fill", "white")
    .on("mousedown", () => {
        isDown = true;
        downPos = [d3.event.offsetX, d3.event.offsetY];
    })
    .on("mouseup", () => {
        isDown = false;
        downPos = undefined;
    })
    .on("mousemove", () => {
        if (!isDown)
            return;
        var pt = [d3.event.offsetX, d3.event.offsetY];
        if (distance(downPos, pt) < minDistance)
            return;
        downPos = pt;
        pts.push(d3.event.offsetX);
        pts.push(d3.event.offsetY);
        updateDelaunay();
    });

function updateDelaunay()
{
    try {
        var del = new Delaunator(pts);
        
        drawDelaunay(del, pts);
    } catch (e) {
        console.log(e);
    }
}

function drawDelaunay(del, pts)
{
    d3.select("#main").selectAll("g").remove();

    var triangleSides = determineTriangleSides(del, pts);

    var circumcircles = d3.range(del.triangles.length / 3).map(function(i) {
        var a = pointFromTri(del, del.triangles[3*i]);
        var b = pointFromTri(del, del.triangles[3*i+1]);
        var c = pointFromTri(del, del.triangles[3*i+2]);
        return [circumcenter(a,b,c), triangleSides[i]];
    });

    var circumCircleSVG = d3.select("#main")
        .append("g")
        .selectAll("circle")
        .data(circumcircles)
        .enter()
        .append("circle")
        .attr("cx", (d) => d[0][0])
        .attr("cy", (d) => d[0][1])
        .attr("r", (d) => d[0][2])
        .attr("stroke", d3.lab(80, 0, 0))
        .attr("pointer-events", "none")
        .attr("fill-opacity", 0.1)
        .attr("fill", (d) => d[1] ? d3.lab(80, 0, 30) : d3.lab(80, 0, -30));
    
    var edges = d3.select("#main")
        .append("g")
        .selectAll("line")
        .data(d3.range(del.triangles.length))
        .enter()
        .append("line")
        .attr("pointer-events", "none")
        .attr("x1", (i) => del.coords[2*del.triangles[i]])
        .attr("y1", (i) => del.coords[2*del.triangles[i]+1])
        .attr("x2", (i) => del.coords[2*del.triangles[nextVertexInTri(i)]])
        .attr("y2", (i) => del.coords[2*del.triangles[nextVertexInTri(i)]+1])
        .attr("stroke", "black")
        .attr("stroke-width", "2px");

    var path = d3.line()
        .x(i => pts[(2*i) % pts.length])
        .y(i => pts[(2*i+1) % pts.length]);
    path = path(d3.range(pts.length/2 + 1));
    
    var polySVG = d3.select("#main")
        .append("g")
        .append("path")
        .attr("pointer-events", "none")
        .attr("d", path)
        .attr("stroke", "red")
        .attr("fill", "none")
        .attr("stroke-width", 4);

    var pointSVG = d3.select("#main")
        .append("g")
        .selectAll("circle")
        .data(del.ids)
        .enter()
        .append("circle")
        .attr("pointer-events", "none")
        .attr("cx", (d) => del.coords[d*2])
        .attr("cy", (d) => del.coords[d*2+1])
        .attr("r", 5)
        .attr("fill", "brown");

    // update three.js scene
    meshes.forEach((m) => scene.remove(m));
    meshes = circumcircles.filter(c => c[1]).map(c => {
        console.log(c);
        var ballMesh = new THREE.Mesh( ballGeometry, material );
        
        ballMesh.position.set((c[0][0] - 300) / 20, (300 - c[0][1]) / 20, 0);
        ballMesh.scale.set(c[0][2] / 20, c[0][2] / 20, c[0][2] / 20);
        ballMesh.castShadow = true;
        return ballMesh;
    });
    meshes.forEach(m => scene.add(m));
}


//////////////////////////////////////////////////////////////////////////////
// three.JS part

var container = document.getElementById( 'three-container' );
camera = new THREE.PerspectiveCamera( 30, 1, 1, 5000 );
camera.position.set( 30, 30, 250 );
camera.lookAt( new THREE.Vector3(0, 0, 0) );
scene = new THREE.Scene();
scene.background = new THREE.Color().setHSL( 0.6, 0, 1 );
scene.fog = new THREE.Fog( scene.background, 1, 5000 );
// LIGHTS
hemiLight = new THREE.HemisphereLight( 0xffffff, 0xffffff, 0.6 );
hemiLight.color.setHSL( 0.6, 1, 0.6 );
hemiLight.groundColor.setHSL( 0.095, 1, 0.75 );
hemiLight.position.set( 0, 50, 0 );
scene.add( hemiLight );
hemiLightHelper = new THREE.HemisphereLightHelper( hemiLight, 10 );
scene.add( hemiLightHelper );
//
dirLight = new THREE.DirectionalLight( 0xffffff, 1 );
dirLight.color.setHSL( 0.1, 1, 0.95 );
dirLight.position.set( -1, 1.75, 1 );
dirLight.position.multiplyScalar( 30 );
scene.add( dirLight );
dirLight.castShadow = true;
dirLight.shadow.mapSize.width = 2048;
dirLight.shadow.mapSize.height = 2048;
var d = 50;
dirLight.shadow.camera.left = -d;
dirLight.shadow.camera.right = d;
dirLight.shadow.camera.top = d;
dirLight.shadow.camera.bottom = -d;
dirLight.shadow.camera.far = 3500;
dirLight.shadow.bias = -0.0001;
dirLightHeper = new THREE.DirectionalLightHelper( dirLight, 10 );
scene.add( dirLightHeper );
// GROUND
var groundGeo = new THREE.PlaneBufferGeometry( 10000, 10000 );
var groundMat = new THREE.MeshPhongMaterial( { color: 0xffffff, specular: 0x050505 } );
groundMat.color.setHSL( 0.095, 1, 0.75 );
var ground = new THREE.Mesh( groundGeo, groundMat );
ground.rotation.x = -Math.PI/2;
ground.position.y = -33;
scene.add( ground );
ground.receiveShadow = true;
// SKYDOME
var vertexShader = document.getElementById( 'vertexShader' ).textContent;
var fragmentShader = document.getElementById( 'fragmentShader' ).textContent;
var uniforms = {
    topColor:    { value: new THREE.Color( 0x0077ff ) },
    bottomColor: { value: new THREE.Color( 0xffffff ) },
    offset:      { value: 33 },
    exponent:    { value: 0.6 }
};
uniforms.topColor.value.copy( hemiLight.color );
scene.fog.color.copy( uniforms.bottomColor.value );
var skyGeo = new THREE.SphereBufferGeometry( 4000, 32, 15 );
var skyMat = new THREE.ShaderMaterial( { vertexShader: vertexShader, fragmentShader: fragmentShader, uniforms: uniforms, side: THREE.BackSide } );
var sky = new THREE.Mesh( skyGeo, skyMat );
scene.add( sky );
// MODEL

var ballGeometry = new THREE.SphereBufferGeometry( 1, 32, 32 );
var material = new THREE.MeshPhysicalMaterial( {
    color: new THREE.Color().setHSL( 0, 0.5, 0.25 ),
    metalness: 0,
    roughness: 0.5,
    clearCoat:  0,
    clearCoatRoughness: 0,
    reflectivity: 0,
    envMap: null
} );

var meshes = [];
var ballMesh = new THREE.Mesh( ballGeometry, material );
ballMesh.position.set(0, 0, 0);
ballMesh.rotation.y = Math.PI;
ballMesh.castShadow = true;
meshes.push(ballMesh);

meshes.forEach(m => scene.add(m));

// RENDERER
renderer = new THREE.WebGLRenderer( { antialias: true } );
renderer.setPixelRatio( window.devicePixelRatio );
renderer.setSize(600, 600);
container.appendChild( renderer.domElement );
renderer.gammaInput = true;
renderer.gammaOutput = true;
renderer.shadowMap.enabled = true;

function animate() {
    requestAnimationFrame( animate );
    render();
}

animate();

function render() {
    renderer.render( scene, camera );
}

