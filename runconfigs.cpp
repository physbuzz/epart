//Pseudocode for running different setups.


namespace RunConfig {
//for example, basic_nparticles(1000,"run_1e3_f"); generates 
//100 640x480 frames of 1000 particles bouncing around.
void basic_nparticles(int nparticles, std::string prefix, int nframes=100) {
    //std::string prefix="run_1e3_";
    //int nparticles=1000;
    //int nframes=100;
    int imgw=640;
    int imgh=480;
    int frameskip=2;

    float temperature=1;
    //Packing fraction. This has to be small!
    float eta=0.06;

    vec2 simulationVolume=vec2({1.0,1.0f*imgh/imgw});
    //imgw,imgh,cx,cy,realsize
    ImageParams imageparams={imgw,imgh,simulationVolume[0]/2,simulationVolume[1]/2,simulationVolume[0]};

    float r=std::sqrt(eta/(nparticles*M_PI));
    float dt=Simulator::findDt(temperature,eta);
    float maxH=Simulator::findMaxH(temperature,eta,dt);

    ParticleList<float,2> plist();
    plist.placeVolume(vec2({0.0f,0.0f}),simulationVolume,nparticles,temperature);

    Simulator s(simulationVolume,maxH);
    s.ingestParticles(plist);

    for(int frame=0;frame<nframes;frame++){
        s.step(dt,frameskip);
        s.saveCircleImage(imageparams,prefix,frame,4);
    }
}

void basic_damBreak(int nparticles, std::string prefix, std::string datafile="", int nframes=1000) {
    //std::string prefix="run_1e3_";
    //int nparticles=1000;
    //int nframes=100;
    int imgw=640;
    int imgh=480;

    //Target so that 1s in the simulation is exactly 60 frames.
    int fps=60; 

    float temperature=1;
    //Packing fraction. This has to be small!
    float eta=0.06;

    vec2 simulationVolume=vec2({1.0,1.0f*imgh/imgw});
    //imgw,imgh,cx,cy,realsize
    ImageParams imageparams={imgw,imgh,simulationVolume[0]/2,simulationVolume[1]/2,simulationVolume[0]};

    float r=std::sqrt(eta/(nparticles*M_PI));
    float dt=Simulator::findDt(temperature,eta);
    float maxH=Simulator::findMaxH(temperature,eta,dt);
    
    int frameskip;
    if(dt>(1.0/fps)){
        frameskip=1;
        dt=1.0/fps;
    } else { //dt<=1/fps
        frameskip=int(std::ceil(1.0/(fps*dt)));
    }

    ParticleList<float,2> plist();
    //fill only the left half of the volume with particles.
    plist.placeVolume(vec2({0.0f,0.0f}),vec2({simulationVolume[0]/2,simulationVolume[1]}),nparticles,temperature);

    Simulator s(simulationVolume,maxH);
    s.ingestParticles(plist);

    std::unique_ptr<std::ofstream> dataout;
    if(datafile.length()>0){
        dataout=std::make_unique<std::ofstream>(datafile);
        (*dataout)<<"t, S, beta, nZero, nOne, nTwoOrMore\n";
    }
    for(int frame=0;frame<nframes;frame++){
        s.step(dt,frameskip);
        s.saveDensityImages(imageparams,prefix,frame,4);
        if(dataout){
            (*dataout)<<(dt*frameskip*(frame+1))<<", ";
            (*dataout)<<s.getEntropy()<<", "<<s.getBeta()<<", "<<s.getNZero()<<", "<<s.getNOne()<<", "<<s.getNTwoOrMore()<<"\n";
        }
    }
    if(dataout) {
        dataout->close();
    }
}

void basic_weighttest(int nparticles, std::string prefix, int nframes=100) {
    //std::string prefix="run_1e3_";
    //int nparticles=1000;
    //int nframes=100;
    int imgw=640;
    int imgh=480;
    int frameskip=2;

    float temperature=1;
    //Packing fraction. This has to be small!
    float eta=0.06;

    vec2 simulationVolume=vec2({1.0,1.0f*imgh/imgw});
    //imgw,imgh,cx,cy,realsize
    ImageParams imageparams={640,480,simulationVolume[0]/2,simulationVolume[1]/2,simulationVolume[0]};

    float r=std::sqrt(eta/(nparticles*M_PI));
    float dt=Simulator::findDt(temperature,eta);
    float maxH=Simulator::findMaxH(temperature,eta,dt);

    ParticleList<float,2> plist();
    plist.placeVolume(vec2({0.0f,0.0f}),simulationVolume,nparticles,temperature);

    Simulator s(simulationVolume,maxH);
    s.ingestParticles(plist);

    for(int frame=0;frame<nframes;frame++){
        s.step(dt,frameskip);
        s.saveCircleImage(imageparams,prefix,frame,4);
    }
}

}


int nparticles=1e6;
float L=3.0;
float H=1.0;
float T=1.0;
float area=L*H;
float nDensity=nparticles/area;

float dt,maxVel,r;
dt=Simulator::recommendDT(nparticles,T,area);
maxVel=Simulator::recommendmaxVel(nparticles,T,area);
r=Simulator::recommendmaxVel(nparticles,T,area);

ParticleList pl;
pl.placeParticleVolume(vec2({0.0f,0.0f}),vec2({L,L}),nDensity,T);

Simulator<float,2> sim;
sim.setDomainSize(vec2({L,H}));
sim.setMaxVel(maxVel);
sim.setMaxH(dt*maxVel);
sim.ingestParticleList(pl);

float throatSize=H/2;
float powModifier=0.5; 
float paddingLeft=0.2;
float paddingBottom=0.0f;
bool addCutout=true;
float cutoutWidth=0.2;
float cutoutHeight=(H-throatSize)/2;
int nverts=50;

Inflow *in=sim.addInflow();
in->setBL(vec2({0.0f,0.0f}));
in->setTR(vec2({0.1f,H}));
in->setPX(1.0f);
in->setPY(0.0f);
in->setDensity(nDensity);
in->setTemperature(T);

Outflow *out=sim.addOutflow();
out->setBL(vec2({L-0.1f,0.0f}));
out->setTR(vec2({L,H}));
out->setPX(1.0f);
out->setPY(0.0f);
out->setDensity(nDensity);

//Bottom part of the wind tunnel
Mesh2D *mesh=sim.addMesh();
mesh->beginShape();
for(int i=0;i<=nverts;i++){
    mesh->vertex(-0.1+(L+0.2)*i/nverts,
}
...
mesh->endShape();

//Top part
mesh=sim.addMesh();
mesh->beginShape();
for(int i=0;i<=nverts;i++){
    mesh->vertex(-0.1+(L+0.2)*i/nverts,
}
...
mesh->endShape();





struct MovingWallParams {
    //Velocity of the floating wall
    float wall_velocity;
    //Total simulated time
    float total_time;
    //Simulated time per frame
    float time_per_frame;
};

void movingWallTest(int nparticles, std::string prefix, 
    float wallVelocity,
    float totalTime=50.0f,
    float timePerFrame=(1.0f/60.0f),
    int imgw=640,
    int imgh=480, 
    bool verbose=true) {

    float temperature=1;
    //Packing fraction. This has to be small!
    float eta=0.06;

    vec2 simulationVolume=vec2({2.0f,1.0f});

    //imgw,imgh,cx,cy,realsize
    ImageParams imageparams={imgw,imgh,simulationVolume[0]/2,simulationVolume[1]/2,simulationVolume[0]};

    //particle radius and timestep
    float r,dt,maxH;
    int nframes,frameskip;

    nframes=int(std::ceil(totalTime/timePerFrame));
    //Satisfy packing fraction constraint (N*pi*r**2 / area == eta)
    r=std::sqrt(eta/(nparticles*M_PI));
    //dt recommended to keep a small number of collisions per timestep.
    float dt1=Simulator::findDt(temperature,eta);
    //dt recommended to evenly get the number of frames.
    float dt2=totalTime/nframes;
    if(dt2<dt1){
        //If we request a high framerate, we can just have 1 simulation step
        //per frame.
        dt=dt2;
        frameskip=1;
    } else {
        //else we have multiple collisions per frame of output, so we need
        //multiple substeps.
        frameskip=std::ceil(dt2/dt1);
        dt=dt2/frameskip;
        //The conditions being satisfied here are
        //dt2 == dt*frameskip, and dt<=dt1
    }
    //Find the grid cell size for the PGrid data structure.
    maxH=Simulator::findMaxH(temperature,eta,dt);

    ParticleList<float,2> plist();
    plist.placeVolume(vec2({0.0f,0.0f}),vec2({1.0f,1.0f}),nparticles,temperature);

    Simulator s(simulationVolume,maxH);
    s.ingestParticles(plist);

    //s::initializeMovingWall(float xpos, float v);
    s.initializeMovingWall(1.0f,wallVelocity);

    //StatisticsGrid finds n,beta,px,py by dividing space into a grid.
    s.initializeStatisticsGrid(imgw/2,imgw/4);
    s.recalculateStatisticsGrid();
    s.findParticleEntropies(); //update internal "entropy per particle" calculation

    //CSV of other stuff (energy, total entropy, etc.)
    CSVData simData("sim_data.csv");
    simData.newColumn("time");
    simData.newColumn("Px");
    simData.newColumn("Py");
    simData.newColumn("N");
    simData.newColumn("N_true");
    simData.newColumn("T");
    simData.newColumn("S");
    simData.newColumn("E");
    //Histograms of parameters with y integrated out. px,py,n,t.
    CSVData pxntSliceHistogram("pxnt_slice_hist.csv");
    //Particle data in the center of the screen
    CSVData pxCenterHistogram("px_hist.csv");
    float pxHistogramMin=-3;
    float pxHistogramMax=-3;
    float pxHistogramBinSize=0.05;
    float pxHistogramSize=0.2;

    #include <cstdarg> 
    CSVData(std::string filename) { }
    void newColumn(std::string name) { }
    void newRow(std::vector<float> arg) { }
    void newRow(std::vector<std::string> arg) { }
    void newRow(float arg, ...) { }
    void save() { }

    //CollisionSimulator
    step(float dt,int nsteps)
    void initializeMovingWall(float xpos, float velocity) { }
    void initializeStatisticsGrid(int nx, int ny){ }
    void recalculateStatisticsGrid() { }
    void findParticleEntropies() { }
    void saveStatsImages(imageparams,prefix,frame,n)
    void saveNewEntropyImage(imageparams,prefix,frame,n)
    float getTime() { }
    float getTotalPx() { }
    float getTotalPy() { }
    float getTotalN() { }
    float getAverageT() { }
    float getTotalS() { }
    float getTotalE() { }
    std::vector<float> getpxntSliceHistogram()
    std::vector<float> getpxntSliceHistogramRowLabels()
    std::vector<float> getpxntSliceHistogramRow1()
    std::vector<float> getpxHistogram(
            VectorND<float,2> bl,
            VectorND<float,2> tr,
            float histmin,
            float histmax,
            float binsize) { }

    for(int frame=0;frame<nframes;frame++){
        s.step(dt,frameskip);
        
        s.recalculateStatisticsGrid();
        s.findParticleEntropies();

        //Plots of px,py,n,T,s
        s.saveStatsImages(imageparams,prefix,frame,4);

        //Plots of the entropy generated
        s.saveEntropyGeneratedImage(imageparams,prefix,frame,4);

        //Save simulation parameters like total energy, 
        simData.newRow(s.getTime(),s.getTotalPx(),s.getTotalPy(),s.getTotalN(),nparticles,s.getAverageT(),s.getTotalS(), s.getTotalE());
        simData.save();

        //Save momentum histogram
        pxntSliceHistogram.newRow(s.getpxntSliceHistogram());
        pxntSliceHistogram.save();

        //Save momentum histogram
        pxCenterHistogram.newRow(s.getpxHistogram(
                    vec2({simulationVolume[0]/2-pxHistogramSize/2,simulationVolume[1]/2-pxHistogramSize/2}),
                    vec2({simulationVolume[0]/2+pxHistogramSize/2,simulationVolume[1]/2+pxHistogramSize/2}),
                    pxHistogramMin,
                    pxHistogramMax,
                    pxHistogramBinSize));
        pxCenterHistogram.save();
    }
}





/*
Initializing the simulation:
Size of the domain, position of any obstacles or movers,
Setup of boundary conditions, initial fill with particleList
Setup of any control regions (eg thermostats) or probes
Definition of dt, maxVel, and r.
Calculation of the relation between physical parameters (eta, mean free path, pressure)

Running the simulation, internals:
Update of any movable objects
Optimized inner loop: update particles, run interactions with boundaries,
run interactions with static and movable meshes.
Handle control regions (inflow/outflow).
Update easy particle stats.
Update easy timing stats.

Running the simulation, externals:
Calculate global stats only if needed, and only the parameters needed.
Saving the stats to a table that can be accessed later
Checkpointing and saving to external files
Calculation of physical params (to understand what should be happening in the program)
Logging to check if there's any craziness happening (high density, high rate of missed collisions)
Drawing:
- Drawing different stats.
- Drawing the particles as circles.
- Camera movement and updating
- Choosing color functions


What things are already totally separate?
- Image, DoubleImage,
- PGrid,
- ParticleList
- VectorND


*/







