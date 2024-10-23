#ifndef COLLISIONSIMULATOR_H
#define COLLISIONSIMULATOR_H

#define COLLISION_EPS2 0.0000001

#include "utils.h"
#include "VectorND.h"
#include "phystructs.h"
#include "ParticleList.h"
#include "PGrid.h"
#include "ImageUtil.h"
#include "easytime.h"



class CollisionSimulator {
public:

    struct ImageParams {
        int imgw;
        int imgh;
        float realsize;
        float cx;
        float cy;
    };
    struct CollisionStats {
        //time spent calculating collisions
        float collisionTime; 
        //time spent in drawing and saving functions
        float drawingTime;        

        //fraction of particles not colliding at all
        int nZero; 
        //fraction with one collision
        int nOne; 
        //fraction with two collisions detected - this should stay low!
        int nTwoOrMore; 

        float currentTime;
        float totalPx;
        float totalPy;
        float totalN;
        float averagen;
        float averageT;
        float totalS;
        float totalE;
        bool statsAreCurrent;
    } stats;

    struct PhysicsQueryStruct {
        float s, n, px, py, beta, e;
    };

    struct StatsGridStruct { 
        float area;
        int sgw,sgh;
        int nref;
        float betaref;
        DoubleImage px,py,n,beta,e,s,ds;
        bool statsAreCurrent;
    } statsgrid;
    void initializeStatsGrid(int nx, int ny,int nref, float betaref){ 
        statsgrid.sgw=nx;
        statsgrid.sgh=ny;
        statsgrid.px=DoubleImage(nx,ny);
        statsgrid.py=DoubleImage(nx,ny);
        statsgrid.n=DoubleImage(nx,ny);
        statsgrid.beta=DoubleImage(nx,ny);
        statsgrid.e=DoubleImage(nx,ny);
        statsgrid.s=DoubleImage(nx,ny);
        statsgrid.ds=DoubleImage(nx,ny);
        statsgrid.statsAreCurrent=false;
        statsgrid.area=s.domainSize[0]*s.domainSize[1]/(statsgrid.sgw*statsgrid.sgh);
        statsgrid.nref=nref;
        statsgrid.betaref=betaref;
        assert(statsgrid.area>0);//,"statsgrid area less than zero. PGrid not initialized?");
    }
    VectorND<int,2> getStatsGridIndex(VectorND<float,2> pos){
        VectorND<int,2> ret({int(std::floor(statsgrid.sgw*pos.x[0]/s.domainSize[0])),
            int(std::floor(statsgrid.sgh*pos.x[1]/s.domainSize[1]))});
        ret.clampCube(VectorND<int,2>({0,0}),VectorND<int,2>({statsgrid.sgw-1,statsgrid.sgh-1}));
        return ret;
    }
    //fract defined as pos == (getStatsGridIndex(pos)+getStatsGridFract(pos))*s.domainSize/statsgrid.sgwh
    VectorND<float,2> getStatsGridFract(VectorND<float,2> pos){
        VectorND<int,2> ind=getStatsGridIndex(pos);
        VectorND<float,2> ret({statsgrid.sgw*pos[0]/s.domainSize[0]-ind[0],
            statsgrid.sgh*pos[1]/s.domainSize[1]-ind[1]});

        return ret;
    }
    float lerp(float a, float b, float t){ return (b-a)*t+a; }
    float statsGridInterpDensity(VectorND<float,2> pos) {
        VectorND<int,2> ind=getStatsGridIndex(pos);
        VectorND<float,2> fract=getStatsGridFract(pos);
        float nll=statsgrid.n.get(ind[0],ind[1]);
        float nlu=statsgrid.n.get(ind[0],std::min(ind[1]+1,statsgrid.sgh-1));
        float nul=statsgrid.n.get(std::min(ind[0]+1,statsgrid.sgw-1),ind[1]);
        float nuu=statsgrid.n.get(std::min(ind[0]+1,statsgrid.sgw-1),std::min(ind[1]+1,statsgrid.sgh-1));
        return lerp(lerp(nll,nul,fract[0]),lerp(nlu,nuu,fract[0]),fract[1]);
    }
    float statsGridInterpDS(VectorND<float,2> pos) {
        VectorND<int,2> ind=getStatsGridIndex(pos);
        VectorND<float,2> fract=getStatsGridFract(pos);
        float nll=statsgrid.ds.get(ind[0],ind[1]);
        float nlu=statsgrid.ds.get(ind[0],std::min(ind[1]+1,statsgrid.sgh-1));
        float nul=statsgrid.ds.get(std::min(ind[0]+1,statsgrid.sgw-1),ind[1]);
        float nuu=statsgrid.ds.get(std::min(ind[0]+1,statsgrid.sgw-1),std::min(ind[1]+1,statsgrid.sgh-1));
        return lerp(lerp(nll,nul,fract[0]),lerp(nlu,nuu,fract[0]),fract[1]);
    }
    float statsGridInterpTemperature(VectorND<float,2> pos) {
        VectorND<int,2> ind=getStatsGridIndex(pos);
        VectorND<float,2> fract=getStatsGridFract(pos);
        float nll=1.0f/statsgrid.beta.get(ind[0],ind[1]);
        float nlu=1.0f/statsgrid.beta.get(ind[0],std::min(ind[1]+1,statsgrid.sgh-1));
        float nul=1.0f/statsgrid.beta.get(std::min(ind[0]+1,statsgrid.sgw-1),ind[1]);
        float nuu=1.0f/statsgrid.beta.get(std::min(ind[0]+1,statsgrid.sgw-1),std::min(ind[1]+1,statsgrid.sgh-1));
        return lerp(lerp(nll,nul,fract[0]),lerp(nlu,nuu,fract[0]),fract[1]);
    }
    void recalculateStatsGrid() {

        statsgrid.px.zeroImage();
        statsgrid.py.zeroImage();
        statsgrid.n.zeroImage();
        statsgrid.e.zeroImage();
        statsgrid.ds.zeroImage();

        stats.totalPx=0;
        stats.totalPy=0;
        stats.totalN=0;
        stats.averagen=0;
        stats.averageT=0;
        stats.totalS=0;
        stats.totalE=0;
        //First, get the totals that we can calculate in the first pass (total momentum, particle number,
        //total energy)
        for(Particle<float,2> &p : *s.plist){
            VectorND<int,2> ind=getStatsGridIndex(p.pos);
            int x=ind[0];
            int y=ind[1];
            statsgrid.n.increase(x,y,1);
            statsgrid.px.increase(x,y,p.vel.x[0]);
            statsgrid.py.increase(x,y,p.vel.x[1]);
            statsgrid.e.increase(x,y,(p.vel.x[0]*p.vel.x[0]+p.vel.x[1]*p.vel.x[1])*0.5f);
        }
        for(int x=0;x<statsgrid.sgw;x++){
            for(int y=0;y<statsgrid.sgh;y++){
                auto &e=statsgrid.e.at(x,y);
                auto &n=statsgrid.n.at(x,y);
                auto &px=statsgrid.px.at(x,y);
                auto &py=statsgrid.py.at(x,y);
                float betadenom=e-0.5*(px*px+py*py)/n; // energy minus energy due to COM motion
                                                       
                auto &beta=statsgrid.beta.at(x,y);
                if(n<=1 || betadenom<=0){
                    beta=std::numeric_limits<float>::infinity();
                    statsgrid.s.put(x,y,0.0f);
                } else {
                    beta=n/betadenom;
                    statsgrid.s.put(x,y,-(n/statsgrid.area)*std::log((n/statsgrid.area)*beta/(statsgrid.nref*statsgrid.betaref)));
                }

                stats.totalE+=e;
                stats.totalPx+=px;
                stats.totalPy+=py;
                stats.totalN+=n;

                if(n>=1){
                    px=px/n;
                    py=py/n;
                }
                e=e/statsgrid.area;
                n=n/statsgrid.area;

                //average density is (sum_i n_i)/N ~= int n(x)^2 d^dx
                //Think of it differently: each particle in this cell has density n.
                //There are  n*area particles in this cell. So the total is going to be Sum(n*n*area)/N
                stats.averagen+=n*n*statsgrid.area;
                stats.averageT+=n*statsgrid.area/beta;
                stats.totalS+=statsgrid.s.get(x,y)*statsgrid.area;
            }
        }
        stats.averagen/=stats.totalN;
        stats.averageT/=stats.totalN;

        float minds=10000,maxds=-10000,avgds=0;

        for(Particle<float,2> &p : *s.plist){
            p.lasts=p.thiss;
            p.thiss=std::log(statsGridInterpTemperature(p.pos)*statsgrid.nref*statsgrid.betaref/(statsGridInterpDensity(p.pos)));

            float ds=(p.thiss-p.lasts);
            if(std::isfinite(ds)){
                avgds+=ds;
                if(ds>maxds)
                    maxds=ds;
                if(ds<minds)
                    minds=ds;
                VectorND<int,2> ind=getStatsGridIndex(p.pos);
                int x=ind[0];
                int y=ind[1];
                if(statsgrid.n.at(x,y)>0){
                    statsgrid.ds.at(x,y)+=ds/(statsgrid.n.at(x,y)*statsgrid.area);
                }
            }
        }
        std::cout<<"(min,max,avg) ds = ("<<minds<<", "<<maxds<<", "<<(avgds/stats.totalN)<<")"<<std::endl;

        stats.statsAreCurrent=true;
        statsgrid.statsAreCurrent=true;
    }



    PGrid<float,2> s;
    float maxH;
    float radius;
    float dt;

    //Conditions for dt and maxh.
    //dt*velocity<maxh
    //Expected number of collisions < some critical number (I won't 
    //guarantee this on this first pass, I guess) 
    //A good rule of thumb is maxh ~= (safety factor)*(diameter)
    //The factor of 6.0f is a safety factor in findDt. It means that we have to be
    //6 standard deviations above the mean in velocity, in order to travel more 
    //than maxH in one timestep.
    static float findMaxH(float radius){
        return 5*radius;
    }
    static float findDt(float temperature,float maxH){
        return maxH/(6.0f*std::sqrt(2.0f*temperature));
    }

    void initializeStats(){
        stats.collisionTime=0; 
        stats.drawingTime=0;        
        stats.nZero=0; 
        stats.nOne=0; 
        stats.nTwoOrMore=0; 

        stats.currentTime=0;
        stats.totalPx=0;
        stats.totalPy=0;
        stats.totalN=0;
        stats.averagen=0;
        stats.averageT=0;
        stats.totalS=0;
        stats.totalE=0;
        stats.statsAreCurrent=false;
    }


    CollisionSimulator(ParticleList<float,2> &pl, VectorND<float,2> domainSize, float maxH, float radius, float dt) :
        stats{},
        s(&pl.plist,domainSize,maxH),
        maxH(maxH),
        radius(radius),
        dt(dt)
    { 
        initializeStats();
    }


    void stepOnce(){
        stats.statsAreCurrent=false;
        statsgrid.statsAreCurrent=false;
        stats.nZero=0;
        stats.nOne=0;
        stats.nTwoOrMore=0;
        stats.currentTime+=dt;

        for(Particle<float,2> &p : *s.plist){
            p.collision=-1;
            p.posnew=VectorND<float,2>({0.0f,0.0f});
            p.velnew=VectorND<float,2>({0.0f,0.0f});
        }

        float radius2=radius*radius;

        //TODO: make this generic
        for(Particle<float,2> *p1 : s.updateLoop()){
            if(p1->collision<0) {
                for(Particle<float,2> *p2 : s.nearbyLoop(p1->pos,maxH)){
                    assert(p2!=nullptr);
                    if(p1==p2 || p2->collision==0)
                        continue;
                    auto x1=p1->pos;
                    auto x2=p2->pos;
                    auto v1=p1->vel;
                    auto v2=p2->vel;
                    auto dx=x2-x1;
                    auto dv=v2-v1;
                    float inner=dx.dot(dv);
                    if(inner>=0)
                        continue;
                    float dv2=dv.length2();
                    float d=inner*inner-dv2*(dx.length2()-4.0f*radius2);
                    if(d<=0)
                        continue;
                    float t1=(-inner-sqrt(d))/dv2;
                    if(t1<COLLISION_EPS2 || t1>dt)
                        continue;
                    //If we get to this point, there's a collision within time dt.
                    //If the other particle has already undergone a collision, also ignore
                    //but set a flag for it.
                    //
                    //Basically, we only ever handle collisions for two pairs of particles
                    //with their collision flags set to -1.
                    if(p2->collision>0){
                        //This is not a perfect count of the number of double collisions. 
                        //Need to double check that it has some bearing to ground truth!
                        stats.nTwoOrMore+=1;
                        stats.nOne-=1;
                        continue;
                    }

                    float t2=dt-t1;
                    p1->collision=1;
                    p2->collision=1;
                    stats.nOne+=2;

                    p1->posnew=p1->pos+p1->vel*t1;
                    p2->posnew=p2->pos+p2->vel*t1;
                    //Collision of equal masses; we reverse each relative velocity along the direction
                    //of their collision vector.
                    //
                    //Start with velocity v1. Put it in CM frame:
                    //v1_CM=v1-(v1+v2)/2
                    //reverse it along the rhat direction:
                    //v1_CM -> v1_CM-2*rhat(rhat.dot(v1_CM))
                    //add (v1+v2)/2 to put it back in non-CM and simplify:
                    //v1_new = v1+2*rhat*rhat.dot((v2-v1)/2)
                    //       = v1+rhat*rhat.dot(v2-v1)
                    //cout<<(p2->posnew-p1->posnew).length()<<endl;
                    auto dx2=(p2->posnew-p1->posnew).normalized();
                    p1->velnew=p1->vel+dx2*dx2.dot(dv);
                    p2->velnew=p2->vel-dx2*dx2.dot(dv);

                    //time evolve the rest of the way.
                    p1->posnew+=p1->velnew*t2;
                    p2->posnew+=p2->velnew*t2;
                    break;
                }
            }

            if(p1->collision<0){
                p1->posnew=p1->pos+dt*p1->vel;
                p1->velnew=p1->vel;
                p1->collision=0;
                stats.nZero+=1;
            }
            p1->pos=p1->posnew;
            p1->vel=p1->velnew;
            if(p1->posnew.x[0]<0){
                p1->vel.x[0]=-p1->vel.x[0];
                p1->pos.x[0]=-p1->pos.x[0];
            }
            if(p1->posnew.x[0]>s.domainSize[0]){
                p1->vel.x[0]=-p1->vel.x[0];
                p1->pos.x[0]=2*s.domainSize[0]-p1->pos.x[0];
            }
            if(p1->posnew.x[1]<0){
                p1->vel.x[1]=-p1->vel.x[1];
                p1->pos.x[1]=-p1->pos.x[1];
            }
            if(p1->posnew.x[1]>s.domainSize[1]){
                p1->vel.x[1]=-p1->vel.x[1];
                p1->pos.x[1]=2*s.domainSize[1]-p1->pos.x[1];
            }
            /* It's still possible to wind up with p1->pos outside the boundaries
             * after these checks, but the PGrid updater will clamp the position to the
             * boundaries
             * */
        } 
    }
    void step(int nsteps){
        for(int i=0;i<nsteps;i++){
            stepOnce();
        }

    }

    /*
void initializeMovingWall(float xpos, float velocity) { }
void findParticleEntropies() { }
void saveStatsImages(imageparams,prefix,frame,n)
void saveNewEntropyImage(imageparams,prefix,frame,n)
std::vector<float> getpxntSliceHistogram()
std::vector<float> getpxntSliceHistogramRowLabels()
std::vector<float> getpxntSliceHistogramRow1()
std::vector<float> getpxHistogram(
        VectorND<float,2> bl,
        VectorND<float,2> tr,
        float histmin,
        float histmax,
        float binsize) { }*/
    //Return physical values after sampling particles within rmax of pos.
    //c is the calculated value from radialwc(rmax,p)
    PhysicsQueryStruct querySimulation(VectorND<float,2> pos,float rmax, float p,float c=-1.0f){
        if(c<=0.0f)
            c=radialwc(rmax,p);

        PhysicsQueryStruct ret{0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};

        for(Particle<float,2> *p2 : s.nearbyLoop(pos,rmax)){
            float r=(p2->pos-pos).length();
            float w=radialw(r,p,rmax);
            ret.n+=w;
            ret.px+=p2->vel.x[0]*w;
            ret.py+=p2->vel.x[1]*w;
            ret.e+=0.5f*p2->vel.length2()*w;
        }
        //If we didn't pick up any particles, just return zero.
        if(ret.n<=COLLISION_EPS2)
            return ret;
        //calculate the expected momentum and expected energy
        ret.px/=ret.n;
        ret.py/=ret.n;
        ret.e/=ret.n;
        float h=0.0f;
        int nparticles=0;
        for(Particle<float,2> *p2 : s.nearbyLoop(pos,rmax)){
            float r=(p2->pos-pos).length();
            float w=radialw(r,p,rmax);
            if(w>0){
                nparticles++;
            }
            float p1x=p2->vel.x[0]-ret.px;
            float p1y=p2->vel.x[1]-ret.py;
            h+=0.5f*(p1x*p1x+p1y*p1y)*w;
        }
        ret.beta=ret.n/h;
        ret.n/=c;
        float z11=100000.0f;
        ret.s=ret.n*(2.0f-log(ret.n*ret.beta/z11));

        if(nparticles<=1){
            //Can't estimate beta and s if there's only one particle
            ret.beta=0.0f;
            ret.s=0.0f;
        }
        return ret;
    }
    void saveStatsgridImage(ImageParams ip,
            std::string prefix, int fnamei, int padcount){

        int imw=ip.imgw;
        int imh=ip.imgh;

        Image outimgDensity(imw,imh);
        Image outimgDS(imw,imh);
        float realsize=ip.realsize;
        float cx=ip.cx;
        float cy=ip.cy;
        float aspect=float(imh)/imw;

        //reference number density.
        //Average should be s.domainSize.product()/s.plist->size().
        float nref=(3.0f*s.plist->size())/s.domainSize.product();
        float totals=0.0f;

        for(int a=0;a<imw;a++){
            for(int b=0;b<imh;b++){
                float x=cx+(float(a)/imw-0.5f)*realsize;
                float y=cy+(float(b)/imh-0.5f)*realsize*aspect;
                VectorND<float,2> pos({x,y});

                float density=statsGridInterpDensity(pos);

                float sc=density/nref;
                outimgDensity.put(a,b,intToRGB(sc*255,sc*255,sc*255));

                if(density>0){
                    //60 is a temporary magic number; it's really frameskip/dt.
                    //to get the rate of change of entropy per particle per frame.
                    float ds=statsGridInterpDS(pos)*30;
                    //float ds=statsGridInterpDS(pos)/(density*statsgrid.area)*30;
                    int red=ds>0?int(ds*10):0;
                    int blue=ds<0?int(-ds*10):0;
                    outimgDS.put(a,b,intToRGB(red,0,blue));
                } else {
                    outimgDS.put(a,b,intToRGB(0,0,0));
                }
            }
        }
        outimgDensity.save(getFilename(prefix,fnamei,padcount,".bmp"));
        outimgDS.save(getFilename(prefix+"ds_",fnamei,padcount,".bmp"));
    }
    static float radialw(float r, float p,float rmax){
        return (r>=rmax)?0.0f:(1.0f/(p+(r/rmax))-1.0f/(p+1.0f));
    }
    static float radialwc(float rmax,float p){
        return rmax*rmax*M_PI*(2.0f-1.0f/(1.0f+p)+2.0f*p*std::log(p/(1.0f+p)));
    }
    void saveDensityImage(float radiusPrime,float p,
            ImageParams ip,
            std::string prefix, int fnamei, int padcount){

        int imw=ip.imgw;
        int imh=ip.imgh;

        Image outimg(imw,imh);
        float realsize=ip.realsize;
        float cx=ip.cx;
        float cy=ip.cy;
        float aspect=float(imh)/imw;
        float cc=radialwc(radiusPrime,p);


        //reference number density.
        //Average should be s.domainSize.product()/s.plist->size().
        float nref=(3.0f*s.plist->size())/s.domainSize.product();
        float totals=0.0f;

        for(int a=0;a<imw;a++){
            for(int b=0;b<imh;b++){
                float x=cx+(float(a)/imw-0.5f)*realsize;
                float y=cy+(float(b)/imh-0.5f)*realsize*aspect;
                VectorND<float,2> pos({x,y});

                auto q=querySimulation(pos,radiusPrime,p,cc);
                float sc=q.n/nref;
                outimg.put(a,b,intToRGB(sc*255,sc*255,sc*255));
                totals+=q.s*realsize*realsize*aspect/(imw*imh);
                /*
                if(!accept){
                    float m=0.9f;
                    float c=std::cos(particleIndex);
                    float s=std::cos(particleIndex);
                    auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                    outimg.put(a,b,intToRGB(rgb.r,rgb.g,rgb.b));
                }
                else
                    outimg.put(a,b,intToRGB(0,0,0));*/
            }
        }
        //auto q=querySimulation(s,VectorND<float,2>({0.5f,0.5f}),radiusPrime,p,cc);
        //cout<<"Entropy density at some pos: "<<q.s<<endl;
        //cout<<"Entropy: "<<totals<<endl;
        outimg.save(getFilename(prefix,fnamei,padcount,".bmp"));
    }

    void saveDensityImages(float radiusPrime, float p,
            ImageParams ip,
            std::string prefix, int fnamei, int padcount, float timevalue){
        int imw=ip.imgw;
        int imh=ip.imgh;

        Image imgDensity(imw,imh);
        Image imgpx(imw,imh);
        Image imgpy(imw,imh);
        Image imgs(imw,imh);
        Image imge(imw,imh);


        float realsize=ip.realsize;
        float cx=ip.cx;
        float cy=ip.cy;
        float aspect=float(imh)/imw;
        float cc=radialwc(radiusPrime,p);


        //reference number density.
        //Average should be s.domainSize.product()/s.plist->size().
        float nref=(3.0f*s.plist->size())/s.domainSize.product();
        float totals=0.0f;



        float nref2=(1.0f*s.plist->size())/s.domainSize.product();
        float z11=10000.0f;
        float sref=nref2*(2.0f-std::log(nref2*1.0f/z11));
        //float srefmax=sref*1.5f;
        //float srefmin=sref/2.0f;
        float srefmax=410000.0f/4.0f;
        float srefmin=340000.0f/32.0f;;
        //cout<<srefmin<<" : "<<srefmax<<endl;
        
        float totaln=0.0f;
        float totale=0.0f;


        for(int a=0;a<imw;a++){
            for(int b=0;b<imh;b++){
                float x=cx+(float(a)/imw-0.5f)*realsize;
                float y=cy+(float(b)/imh-0.5f)*realsize*aspect;
                VectorND<float,2> pos({x,y});

                auto q=querySimulation(pos,radiusPrime,p,cc);
                totals+=q.s*realsize*realsize*aspect/(imw*imh);
                totaln+=q.n*realsize*realsize*aspect/(imw*imh);
                totale+=q.e*q.n*realsize*realsize*aspect/(imw*imh);

                //Density map:
                float sc=q.n/nref;
                imgDensity.put(a,b,intToRGB(sc*255,sc*255,sc*255));

                //momentum map:
                sc=q.px;
                int red=sc>0?int(std::log(1+sc)*90):0;
                int blue=sc<0?int(std::log(1-sc)*90):0;
                imgpx.put(a,b,intToRGB(red,0,blue));
                //y momentum
                sc=q.py/5.0;
                //red=sc>0?int(sc*255):0;
                //blue=sc<0?int(-sc*255):0;
                red=sc>0?int(std::log(1+sc)*255):0;
                blue=sc<0?int(std::log(1-sc)*255):0;
                imgpy.put(a,b,intToRGB(red,0,blue));

                //entropy picture
                sc=(q.s-srefmin)/(srefmax-srefmin);
                imgs.put(a,b,intToRGB(sc*255,sc*255,sc*255));

                //energy picture
                sc=std::log(1+q.e)/5.0;
                imge.put(a,b,intToRGB(sc*255,sc*255,sc*255));

                /*
                if(!accept){
                    float m=0.9f;
                    float c=std::cos(particleIndex);
                    float s=std::cos(particleIndex);
                    auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                    outimg.put(a,b,intToRGB(rgb.r,rgb.g,rgb.b));
                }
                else
                    outimg.put(a,b,intToRGB(0,0,0));*/
            }
        }
        //auto q=querySimulation(s,VectorND<float,2>({0.5f,0.5f}),radiusPrime,p,cc);
        //cout<<"Entropy density at some pos: "<<q.s<<endl;
        std::cout<<timevalue<<", "<<totals<<", "<<totaln<<", "<<totale<<", ";
        imgDensity.save(getFilename(prefix+"density",fnamei,padcount,".bmp"));
        imgpx.save(getFilename(prefix+"px",fnamei,padcount,".bmp"));
        imgpy.save(getFilename(prefix+"py",fnamei,padcount,".bmp"));
        imgs.save(getFilename(prefix+"s",fnamei,padcount,".bmp"));
        imge.save(getFilename(prefix+"e",fnamei,padcount,".bmp"));
    }

    void saveImage(float radius, 
            ImageParams ip,
            std::string prefix, int fnamei, int padcount){
        float drawR=radius;


        int imw=ip.imgw;
        int imh=ip.imgh;

        Image outimg(imw,imh);
        float realsize=ip.realsize;
        float cx=ip.cx;
        float cy=ip.cy;
        float aspect=float(imh)/imw;

        float drawRSquared=drawR*drawR;
        int nParticlesChecked=0;
        for(int a=0;a<imw;a++){
            for(int b=0;b<imh;b++){
                float x=cx+(float(a)/imw-0.5f)*realsize;
                float y=cy+(float(b)/imh-0.5f)*realsize*aspect;
                VectorND<float,2> pos({x,y});

                float x2=cx+(float(a+1)/imw-0.5f)*realsize;
                float y2=cy+(float(b+1)/imh-0.5f)*realsize*aspect;
                VectorND<float,2> pos2({x2,y2});
                if(!(s.positionToIntvec(pos)==s.positionToIntvec(pos2))){
                    outimg.put(a,b,intToRGB(255,255,255));
                    continue;
                }

                bool accept=true;
                for(Particle<float,2> *p2 : s.nearbyLoop(VectorND<float,2>({x,y}),radius)){
                    nParticlesChecked++;

                    if((pos-p2->pos).length2()<drawRSquared) {
                        accept=false;
                        break;
                    }
                }
                if(!accept){
                    /*float m=0.9f;
                    float c=std::cos(particleIndex);
                    float s=std::cos(particleIndex);
                    auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                    outimg.put(a,b,intToRGB(rgb.r,rgb.g,rgb.b));*/
                    outimg.put(a,b,intToRGB(255,255,255));
                }
                else
                    outimg.put(a,b,intToRGB(0,0,0));
            }
        }
        //cout<<"Checked "<<(float(nParticlesChecked)/(imw*imh))<<" particles per pixel"<<endl;
        outimg.save(getFilename(prefix,fnamei,padcount,".bmp"));
    }
};

#endif

