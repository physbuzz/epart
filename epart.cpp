#include <iostream>
#include <vector>
#include <cmath>

#include "VectorND.h"
#include "phystructs.h"
#include "ParticleList.h"
#include "PGrid.h"
#include "ImageUtil.h"

#define EPS2 0.00000

//Some string manipulation functions for saving files. pad_int(1234,5) returns "01234".
std::string pad_int(int arg, int padcount) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(padcount) << arg;
    return ss.str();
}

//Returns a file name in the form of "prefix00###suffix". For example "image0032.bmp"
std::string getFilename(std::string prefix, int num, int padcount, std::string suffix) {
    return prefix + pad_int(num, padcount) + suffix;
}


using namespace std;

struct ImageParams {
    int imgw;
    int imgh;
    float realsize;
    float cx;
    float cy;
};

void saveImage(PGrid<float,2> &s,ParticleList<float,2> &pl,float radius, 
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
    for(int a=0;a<imw;a++){
        for(int b=0;b<imh;b++){
            float x=cx+(float(a)/imw-0.5f)*realsize;
            float y=cy+(float(b)/imh-0.5f)*realsize*aspect;

            bool accept=true;
            int particleIndex=0;
            VectorND<float,2> pos({x,y});
            VectorND<int,2> avec=s.positionToIntvec(pos);
            //loop over all adjacent cells
            for(int dx=-1;dx<=1 && accept;dx++){
                for(int dy=-1;dy<=1 && accept;dy++){
                    if(avec[0]+dx<0
                            ||avec[0]+dx>=s.numCells[0]
                            ||avec[1]+dy<0
                            ||avec[1]+dy>=s.numCells[1])
                        continue;

                    //loop over the particles in the adjacent cells
                    int bind=s.intvectorToIndex(avec+VectorND<int,2>({dx,dy}));
                    for(size_t p2ind=0;p2ind<s.idarr[bind].size();p2ind++){
                        Particle<float,2> *p2=&((*s.plist)[s.idarr[bind][p2ind]]);
                        if((pos-p2->pos).length2()<drawRSquared) {
                            accept=false;
                            particleIndex=s.idarr[bind][p2ind];
                            break;
                        }
                    }
                }
            }
            if(!accept){
                float m=0.9f;
                float c=std::cos(particleIndex);
                float s=std::cos(particleIndex);
                auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                outimg.put(a,b,intToRGB(rgb.r,rgb.g,rgb.b));
            }
            else
                outimg.put(a,b,intToRGB(0,0,0));
        }
    }
    outimg.save(getFilename(prefix,fnamei,padcount,".bmp"));
}
void drawImage()
{
    /*
    int imw=3840;
    int imh=2160;
    float drawR=radius*0.9;
    

    DoubleImage dimg(imw,imh);
    float realsize=0.9;
    float cx=0.5;
    float cy=0.5;
    float aspect=float(imh)/imw;

    float drawRSquared=drawR*drawR;
    for(int a=0;a<imw;a++){
        for(int b=0;b<imh;b++){
            float x=cx+(float(a)/imw-0.5f)*realsize;
            float y=cy+(float(b)/imh-0.5f)*realsize*aspect;

            bool accept=true;
            VectorND<float,2> pos({x,y});
            VectorND<int,2> avec=s.positionToIntvec(pos);
            //loop over all adjacent cells
            for(int dx=-1;dx<=1 && accept;dx++){
                for(int dy=-1;dy<=1 && accept;dy++){
                    if(avec[0]+dx<0
                            ||avec[0]+dx>=s.numCells[0]
                            ||avec[1]+dy<0
                            ||avec[1]+dy>=s.numCells[1])
                        continue;

                    //loop over the particles in the adjacent cells
                    int bind=s.intvectorToIndex(avec+VectorND<int,2>({dx,dy}));
                    for(size_t p2ind=0;p2ind<s.idarr[bind].size();p2ind++){
                        Particle<float,2> *p2=&((*s.plist)[s.idarr[bind][p2ind]]);
                        if((pos-p2->pos).length2()<drawRSquared) {
                            accept=false;
                            break;
                        }
                    }
                }
            }
            if(accept)
                dimg.put(a,b,0.0);
            else
                dimg.put(a,b,1.0);
        }
    }
    std::cout<<"Done, saving."<<std::endl;
    //dimg.unitStretch();
    Image outimg(dimg.getData(),imw,imh);
    outimg.save("tmp.bmp");
    std::cout<<"Done."<<std::endl;*/
}


void updateOnce(PGrid<float,2> &s,float radius,float dt){
    int nZero=0;
    int nOne=0;
    int nTwo=0;

    for(Particle<float,2> &p : *s.plist){
        p.collision=-1;
        p.posnew=VectorND<float,2>({0.0f,0.0f});
        p.velnew=VectorND<float,2>({0.0f,0.0f});
    }

    float radius2=radius*radius;
    for(Particle<float,2> *p1 : s.updateLoop()){
        int aind=s.getParticleIndex(*p1);
        VectorND<int,2> avec=s.indexToIntvector(aind);

        //loop over all adjacent cells
        for(int dx=-1;dx<=1 && (p1->collision<0);dx++){
            for(int dy=-1;dy<=1 && (p1->collision<0);dy++){
                //ignore cells that are out of bounds
                if(avec[0]+dx<0
                        ||avec[0]+dx>=s.numCells[0]
                        ||avec[1]+dy<0
                        ||avec[1]+dy>=s.numCells[1])
                    continue;
                
                //loop over the particles in the adjacent cells
                int bind=s.intvectorToIndex(s.indexToIntvector(aind)+VectorND<int,2>({dx,dy}));
                for(size_t p2ind=0;p2ind<s.idarr[bind].size();p2ind++){
                    Particle<float,2> *p2=&((*s.plist)[s.idarr[bind][p2ind]]);
                    if(p1==p2 || p2->collision>=0)
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
                    if(t1<EPS2 || t1>dt)
                        continue;
                    //If we get to this point, there's a collision within time dt.
                    //If the other particle has already undergone a collision, also ignore
                    //but set a flag for it.
                    //
                    //Basically, we only ever handle collisions for two pairs of particles
                    //with their collision flags set to -1.
                    if(p2->collision>0){
                        nTwo++;
                        continue;
                    }

                    float t2=dt-t1;
                    p1->collision=1;
                    p2->collision=1;
                    nOne++;

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
                    //cout<<(p1->vel.length2()+p2->vel.length2()-p1->velnew.length2()-p2->velnew.length2())<<endl;
                    break;
                }
            }
        }

        if(p1->collision<0){
            p1->posnew=p1->pos+dt*p1->vel;
            p1->velnew=p1->vel;
            p1->collision=0;
        }
        p1->pos=p1->posnew;
        p1->vel=p1->velnew;
        if(p1->posnew.x[0]<0){
            p1->vel.x[0]=-p1->vel.x[0];
            p1->pos.x[0]=-p1->pos.x[0];
        }
        if(p1->posnew.x[0]>s.domainMax[0]){
            p1->vel.x[0]=-p1->vel.x[0];
            p1->pos.x[0]=2*s.domainMax[0]-p1->pos.x[0];
        }
        if(p1->posnew.x[1]<0){
            p1->vel.x[1]=-p1->vel.x[1];
            p1->pos.x[1]=-p1->pos.x[1];
        }
        if(p1->posnew.x[1]>s.domainMax[1]){
            p1->vel.x[1]=-p1->vel.x[1];
            p1->pos.x[1]=2*s.domainMax[1]-p1->pos.x[1];
        }
        /* It's still possible to wind up with p1->pos outside the boundaries
         * after these checks, but the PGrid updater will clamp the position to the
         * boundaries
         * */
    } 
}


int main(){
    float temperature=1.0f;
    //Expected velocities are sqrt(2T/m)
    //time to cross a boundary ~= dx/sqrt(2T/m)
    int nparticles=1000000;




    float radius=0.0002f;
    float L=2.0f;
    VectorND<float,2> domainSize({2.0f*L,L});
    float maxH=0.001f;
    float dt=maxH/(6.0f*std::sqrt(2.0f*temperature));
    

    ParticleList<float,2> pl;
    
    const int maxtries=500;
    for(int i=0;i<nparticles;i++){
        int j=0;
        for(;j<maxtries;j++){
            VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
            int k=0;
            for(;k<pl.plist.size();k++){
                if((pnew-pl.plist[k].pos).length2()<4*radius*radius){
                    break;
                }
            }
            if(k==pl.plist.size()){
                float vmag=std::sqrt(2*temperature);
                float theta=(rand()*2.0f*M_PI)/RAND_MAX;
                VectorND<float,2> vnew({vmag*std::cos(theta),vmag*std::sin(theta)});
                pl.plist.push_back(Particle<float,2>{pnew,vnew,pnew,vnew,-1});
                break;
            } 
        }
        if(j==maxtries) {
            cout<<"Tried maxtries times and couldn't place new particle."<<endl;
        }
    }

    /*
            VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
    VectorND<float,2> pos1({1.0f-0.1,1.0f});
    VectorND<float,2> vel1({1.0f,0.0f});
    VectorND<float,2> pos2({1.0f+0.1,1.0f});
    VectorND<float,2> vel2({-1.0f,0.0f});
    pl.plist.push_back(Particle<float,2>{pos1,vel1,pos1,vel1,-1});
    pl.plist.push_back(Particle<float,2>{pos2,vel2,pos2,vel2,-1});*/

    PGrid<float,2> s(&pl.plist,domainSize,maxH);
    s.rebuildGrid();


    int nframes=200;
    int frameskip=1;
    ImageParams ip{};
    ip.imgw=640;
    ip.imgh=480;
    ip.realsize=0.05f;
    ip.cx=2.0f;
    ip.cy=1.0f;
    for(int i=0;i<=nframes*frameskip;i++){
        updateOnce(s,radius,dt);
        if(i%frameskip==0){
            float e=0;
            for(int j=0;j<s.plist->size();j++){
                e+=0.5*s.plist->at(j).vel.length2();
            }
            cout<<"On step "<<i<<". Energy is: "<<e<<endl;
            saveImage(s,pl,radius,ip,"first",(i/frameskip),5);
        }
    }


    /*
    
    float ratio=0;
    int nframes=600;
    int frameskip=10;
    for(int i=0;i<=nframes*frameskip;i++){
        ratio+=updateOnce(s,fract*radius,0.1*radius);
        if(i%frameskip==0){
            cout<<"On step "<<i<<endl;
            saveImage(s,pl,radius,fract,(i/frameskip));
        }
    }
    cout<<"Acceptance ratio: "<<(ratio/(frameskip*nframes))<<endl;*/



    return 0;
}





