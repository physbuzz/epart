How do I want the API to behave?


ParticleList
placeVolume
save
load


// Run configuration for 10, 100, 1000 particles in a box
std::string prefix="run_1e3_";
int imgw=640;
int imgh=480;
float temperature=1;
vec2 simulationVolume=vec2({1.0,1.0f*imgh/imgw});
int nparticles=1000;
ImageParams p={640,480,



