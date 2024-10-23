#include <iostream>
#include <vector>
#include <cmath>

#include "utils.h"
#include "VectorND.h"
#include "phystructs.h"
#include "ParticleList.h"
#include "PGrid.h"
#include "ImageUtil.h"
#include "easytime.h"
#include "CSVData.h"
#include "CollisionSimulator.h"




using namespace std;

void movingWallTest(int nparticles, std::string prefix, 
    float wallVelocity,
    float totalTime=50.0f,
    float timePerFrame=(1.0f/30.0f),
    int imgw=640,
    int imgh=480, 
    bool verbose=true) {

    float temperature=0.5f;
    //Packing fraction. This has to be small!
    float eta=0.06;
    float L=2.0f;

    VectorND<float,2> simulationVolume=VectorND<float,2>({1.2*L,L});

    //imgw,imgh,cx,cy,realsize
    CollisionSimulator::ImageParams imageparams={imgw,imgh,simulationVolume[0],simulationVolume[0]/2,simulationVolume[1]/2};

    //particle radius and timestep
    float radius,dt,maxH;
    int nframes,frameskip;

    nframes=int(std::ceil(totalTime/timePerFrame));
    //Satisfy packing fraction constraint (N*pi*r**2 / area == eta)
    radius=std::sqrt(eta/(nparticles*M_PI));

    //Find the grid cell size for the PGrid data structure.
    maxH=CollisionSimulator::findMaxH(radius);

    //dt recommended to keep a small number of collisions per timestep.
    float dt1=CollisionSimulator::findDt(temperature,maxH);

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


    ParticleList<float,2> pl;


    //This really should be its own method, but usage of PGrid is awkward here.
    //plist.placeVolume(VectorND<float,2>({0.0f,0.0f}),VectorND<float,2>({1.0f,1.0f}),nparticles,temperature);
    //Simulator s(simulationVolume,maxH);
    //s.ingestParticles(plist);
    for(int i=0;i<nparticles;i++) {
        float vmag=std::sqrt(2*temperature);
        float theta=(rand()*2.0f*M_PI)/RAND_MAX;
        VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
        VectorND<float,2> vnew({vmag*std::cos(theta),vmag*std::sin(theta)});
        pl.plist.push_back(Particle<float,2>{pnew,vnew,pnew,vnew,-1});
    }
    CollisionSimulator cl(pl,simulationVolume,maxH,radius,dt);
    PGrid<float,2> &s=cl.s;
    s.rebuildGrid();
    for(int passes=0;passes<5;passes++){
        for(Particle<float,2> *p1 : s.updateLoop()){
            bool collisionFree=true;
            for(Particle<float,2> *p2 : s.nearbyLoop(p1->pos,2*radius)){
                if(p1==p2)
                    continue;
                if((p1->pos-p2->pos).length2()<4*radius*radius){
                    collisionFree=false;
                    break;
                }
            }
            if(!collisionFree){
                VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
                p1->pos=pnew;
                p1->posnew=pnew;
            }
        }
    }

    //s::initializeMovingWall(float xpos, float v);
    //s.initializeMovingWall(1.0f,wallVelocity);

    //StatisticsGrid finds n,beta,px,py by dividing space into a grid.
    //statsGridN controls the supersampling of cells. Larger statsGridN = bigger cells.
    //The statistics get better the more particles are in a single statsgrid cell, which
    //will be about num_particles_in_cell = nparticles*(statsGridN/imgw)**2
    float aspect=simulationVolume[1]/simulationVolume[0];
    int statsGridN=6; 
    cl.initializeStatsGrid(imgw/statsGridN,int(imgw*aspect/statsGridN), nparticles/(L*L),1.0f);
    cl.recalculateStatsGrid();
    //s.findParticleEntropies(); //update internal "entropy per particle" calculation

    //CSV of other stuff (energy, total entropy, etc.)
    CSVData simData("sim_data.csv");
    simData.newColumn("time");
    simData.newColumn("Px");
    simData.newColumn("Py");
    simData.newColumn("averagen");
    simData.newColumn("N");
    simData.newColumn("N_true");
    simData.newColumn("T");
    simData.newColumn("S");
    simData.newColumn("E");

    /*
    //Histograms of parameters with y integrated out. px,py,n,t.
    CSVData pxntSliceHistogram("pxnt_slice_hist.csv");
    //Particle data in the center of the screen
    CSVData pxCenterHistogram("px_hist.csv");
    float pxHistogramMin=-3;
    float pxHistogramMax=-3;
    float pxHistogramBinSize=0.05;
    float pxHistogramSize=0.2;
    */

    std::cout<<"Rendering "<<nframes<<" frames with "<<frameskip<<" substeps per frame."<<std::endl;

    for(int frame=0;frame<nframes;frame++){
        cl.step(frameskip);
        
        cl.recalculateStatsGrid();
        //s.findParticleEntropies();

        //Plots of px,py,n,T,s
        //s.saveStatsImages(imageparams,prefix,frame,4);

        //Plots of the entropy generated
        //s.saveEntropyGeneratedImage(imageparams,prefix,frame,4);

        //Save simulation parameters like total energy, 
        //simData.newRow(s.getTime(),s.getTotalPx(),s.getTotalPy(),s.getTotalN(),nparticles,s.getAverageT(),s.getTotalS(), s.getTotalE());
        simData.newRow();
        simData.newColumn(cl.stats.currentTime);
        simData.newColumn(cl.stats.totalPx);
        simData.newColumn(cl.stats.totalPy);
        simData.newColumn(cl.stats.averagen);
        simData.newColumn(cl.stats.totalN);
        simData.newColumn(nparticles);
        simData.newColumn(cl.stats.averageT);
        simData.newColumn(cl.stats.totalS);
        simData.newColumn(cl.stats.totalE);
        simData.save();

        float pp=0.05f;
        float rr=0.08f;
        //cl.saveDensityImages(rr,pp,imageparams,"run_",frame,5,dt*frame);
        cl.saveStatsgridImage(imageparams,"statsgrid_",frame,5);

        //Save momentum histogram
        //pxntSliceHistogram.newRow(s.getpxntSliceHistogram());
        //pxntSliceHistogram.save();

        //Save momentum histogram
        /*
        pxCenterHistogram.newRow(s.getpxHistogram(
                    VectorND<float,2>({simulationVolume[0]/2-pxHistogramSize/2,simulationVolume[1]/2-pxHistogramSize/2}),
                    VectorND<float,2>({simulationVolume[0]/2+pxHistogramSize/2,simulationVolume[1]/2+pxHistogramSize/2}),
                    pxHistogramMin,
                    pxHistogramMax,
                    pxHistogramBinSize));
        pxCenterHistogram.save();
        */
    }
}



    /*
void placeVolumeInPlist(ParticleList<float,2> &plist, VectorND<float,2> bl, VectorND<float,2> tr,
        float radius,int nparticles, float temperature){

    for(int i=0;i<nparticles;i++){
        float vmag=std::sqrt(2*temperature);
        float theta=(rand()*2.0f*M_PI)/RAND_MAX;
        VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
        VectorND<float,2> vnew({vmag*std::cos(theta),vmag*std::sin(theta)});
        plist.push_back(Particle<float,2>{pnew,vnew,pnew,vnew,-1});
    }

    VectorND<float,2> diff=tr-bl;

    PGrid<float,2> pg

    PGrid(std::vector<Particle<Float,DIM> > *plist,VectorND<Float,DIM> domainSize,Float maxH) : 
    CollisionSimulator cl(pl,domainSize,maxH);
    PGrid<float,2> &s=cl.s;
    s.rebuildGrid();
    for(int passes=0;passes<5;passes++){
        for(Particle<float,2> *p1 : s.updateLoop()){
            bool collisionFree=true;
            for(Particle<float,2> *p2 : s.nearbyLoop(p1->pos,2*radius)){
                if(p1==p2)
                    continue;
                if((p1->pos-p2->pos).length2()<4*radius*radius){
                    collisionFree=false;
                    break;
                }
            }
            if(!collisionFree){
                VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
                p1->pos=pnew;
                p1->posnew=pnew;
            }
        }
    }
}*/
int main() {

    movingWallTest(1000000,"test",0.0f,50.0f);
    return 0;
}

