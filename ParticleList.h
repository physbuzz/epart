#ifndef PARTICLELIST_H
#define PARTICLELIST_H

#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "phystructs.h"

template<typename Float,int DIM>
class ParticleList { 
public:

    std::vector<Particle<Float,DIM> > plist;
    ParticleList() : plist() { }
    //initialize a numCells[0] by numCells[1] by ... grid, each point a distance w away 
    //from the next.
    void initializeGrid(VectorND<int,DIM> numCells, Float w){
        int maxindex=numCells.product();
        plist.reserve(maxindex);
        for(int index=0;index<maxindex;index++){

            int indexprime=index;
            Particle<Float,DIM> p;
            for(int i=0;i<DIM;i++){
                p.pos[i]=w*(indexprime%numCells[i]);
                indexprime=indexprime/numCells[i];
            }
            plist.push_back(p);
        }
    }
    //initialize an N by N by ... grid, each a distance w away.
    void initializeGrid(int N, Float w){
        VectorND<int,DIM> numCells;
        for(int i=0;i<DIM;i++)
            numCells[i]=N;
        initializeGrid(numCells,w);
    }
    void initializeRandom(int N, Float L){
        for(int n=0;n<N;n++){
            Particle<Float,DIM> p;
            for(int i=0;i<DIM;i++){
                p.pos[i]=(L*rand())/RAND_MAX;
            }
            plist.push_back(p);
        }
    }

    void save(std::string filename){
        std::ofstream pfile;
        pfile<<std::setprecision(9);
        pfile.open (filename.c_str());
        if(pfile.is_open()){
            for(int i=0;i<plist.size();i++){
                auto p=plist.at(i);
                for(int d=0;d<DIM;d++){
                    pfile<<p.pos[d];
                    if(d<DIM-1)
                        pfile<<" ";
                    else 
                        pfile<<"\n";
                }
                for(int d=0;d<DIM;d++){
                    pfile<<p.vel[d];
                    if(d<DIM-1)
                        pfile<<" ";
                    else 
                        pfile<<"\n";
                }
            }
            pfile.close();
        } else {std::cerr<<"Could not open file "<<filename<<" for particle data writing";}
    }

    bool load(std::string filename){
        std::ifstream pfile;
        pfile.open (filename.c_str());
        if(pfile.is_open()){
            plist=std::vector<Particle<Float,DIM> >();
            for(std::string line; std::getline(pfile, line); )   //read stream line by line
            {
                std::istringstream in(line);      //make a stream for the line itself
                Particle<Float,DIM> p;
                for(int d=0;d<DIM;d++){
                    in>>p.pos[d];
                }
                if(!in.fail())
                    plist.push_back(p);
                else {
                    std::cerr<<"ParticleList::Load failed to read line: "<<line<<std::endl;
                }
            }
            pfile.close();
        } else {
            std::cerr<<"Could not open file "<<filename<<" for particle data writing";
            return false;
        }
        return true;
    }
};
#endif
