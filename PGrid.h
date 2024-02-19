#ifndef PGRID_H
#define PGRID_H
#include <utility>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "phystructs.h"


template<typename Float,int DIM>
class SlowParticleManager {
public:
    std::vector<Particle<Float,DIM> > *plist;
    SlowParticleManager() : plist(nullptr) { }


    class SlowParticlePairIterator {
        SlowParticleManager<Float,DIM> *p;
        public:
        size_t i;
        size_t j;
        SlowParticlePairIterator(SlowParticleManager<Float,DIM> *k,size_t position1=0,size_t position2=1) : p(k), i(position1), j(position2){ }
        SlowParticlePairIterator& operator++() {
            i++;
            if(i>=j) {
                i=0;
                j++;
            }
            return *this;
        }
        std::pair<Particle<Float,DIM>*,Particle<Float,DIM>*> operator*() const {
            return std::make_pair(&(*p->plist)[i],&(*p->plist)[j]);
        }
        bool operator!=(const SlowParticlePairIterator &other){
            return (i!=other.i || j!=other.j);
        }
    };

    class SlowParticlePairRange{
        SlowParticleManager<Float,DIM> *p;
        public:
        SlowParticlePairRange(SlowParticleManager<Float,DIM> *p) : p(p){ }

        SlowParticlePairIterator begin() { 
            return SlowParticlePairIterator(p,0,1);
        }
        SlowParticlePairIterator end() { 
            return SlowParticlePairIterator(p,0,p->plist->size());
        }
    };

    SlowParticlePairRange pairs(){
        return SlowParticlePairRange(this);
    }

};

template<typename Float,int DIM>
class PGrid {
public:
    static_assert(0<DIM,"PGrid DIM must be positive.");

    std::vector<std::vector<size_t> > idarr;
    VectorND<int,DIM> numCells;
    VectorND<Float,DIM> boxWidth;
    VectorND<Float,DIM> domainSize;
    VectorND<Float,DIM> domainMax;//iota smaller than domainSize.
    bool needsRebuild;
    size_t productOfSizes;
    Float maxH;

    //sanity check: for DIM=3...
    //ret=arg[2]
    //ret=ret*numCells[1]+arg[1]
    //ret=ret*numCells[0]+arg[0]
    //so ret=arg[2]*numCells[1]*numCells[0]+arg[1]*numCells[0]+arg[0]
    // the index is k*I*J+j*I+i
    int intvectorToIndex(VectorND<int,DIM> arg){
        int ret=0;
        for(int i=DIM-1;i>=0;i--)
            ret=ret*numCells[i]+arg[i];

        return ret;
    }
    //sanity check: for DIM=3...
    //the index is ret[2]*numCells[1]*numCells[0]+ret[1]*numCells[0]+ret[0]
    //so the two functions are inverses of each other.
    VectorND<int,DIM> indexToIntvector(int index){
        VectorND<int,DIM> ret;
        for(int i=0;i<DIM;i++){
            ret[i]=index%numCells[i];
            index=index/numCells[i];
        }
        return ret;
    }
    //Update the particle at (*plist)[pind] from grid position aind to position bind.
    void particleChangeGrid(size_t pind, size_t aind, size_t bind){
    //std::vector<std::vector<size_t> > idarr;
        assert(aind<idarr.size()&&bind<idarr.size());
        std::vector<size_t>::iterator pos=std::find(idarr[aind].begin(),idarr[aind].end(),pind);
        assert(pos!=idarr[aind].end());
        idarr[aind].erase(pos);
        idarr[bind].push_back(pind);
    }
    std::vector<Particle<Float,DIM> > *plist;

    PGrid &operator=(const PGrid&) = delete;
    PGrid(const PGrid&) = delete;
    PGrid() : plist(nullptr) { }
    PGrid(std::vector<Particle<Float,DIM> > *plist) : idarr(), numCells(), boxWidth(), domainSize(), 
        domainMax(), needsRebuild(true), productOfSizes(1), maxH(1), plist(nullptr) { }

    PGrid(std::vector<Particle<Float,DIM> > *plist,VectorND<Float,DIM> domainSize,Float maxH) : 
        idarr(), numCells(), boxWidth(), domainSize(), 
        domainMax(), needsRebuild(true), productOfSizes(0), 
        maxH(maxH), plist(plist) {
        setParams(domainSize,maxH);
    }

    void setParams(VectorND<Float,DIM> domainSize,Float maxH){
        assert(maxH>0 &&domainSize.min()>Float(0));
        this->maxH=maxH;
        this->productOfSizes=1;
        //The guarantee that we want is that:
        //  1. With clamped boundary conditions, we apply "if(particle.x[i]>domainMax[i]) particle.x[i]=domainMax[i];", 
        //  and for periodic boundary conditions, "if(particle.x[i]>domainMax[i]) particle.x[i]=std::fmod(particle.x[i],domainSize[i])"
        //  2. After applying boundary conditions, 0<=size_t(particle.x[i]/boxWidth[i])<numCells[i].
        for(size_t i=0;i<DIM;i++){
            this->numCells[i]=size_t(domainSize[i]/maxH);
            this->numCells[i]=std::max(this->numCells[i],1);
            this->boxWidth[i]=domainSize[i]/this->numCells[i];
            this->domainSize[i]=this->boxWidth[i]*this->numCells[i]; //condition 2
            this->domainMax[i]=std::nextafter(0.9999*domainSize[i],Float(0)); // condition 1
            productOfSizes*=this->numCells[i];
        }
        needsRebuild=true;
    }

    void ingestParticles(std::vector<Particle<Float,DIM> > *plist){
        this->plist=plist;
        needsRebuild=true;
    }


    VectorND<int,DIM> positionToIntvec(const VectorND<Float,DIM> &arg){
        VectorND<int,DIM> ivec;
        for(int i=0;i<DIM;i++){
            ivec[i]=std::floor(arg[i]/boxWidth[i]);
        }
        return ivec;
    }
    size_t getParticleIndex(Particle<Float,DIM> &arg){
        VectorND<int,DIM> ivec;
        for(int i=0;i<DIM;i++){
            if(!((arg.pos[i]/boxWidth[i])<numCells[i])){
                std::cout<<"Pesky out of bounds bug."<<std::endl;
                std::cout<<"i: "<<i<<std::endl;
                std::cout<<"pos: "<<arg.pos<<std::endl;
                std::cout<<"ratio[i]: "<<int(std::floor(arg.pos[i]/boxWidth[i]))<<std::endl;
                std::cout<<"numCells[i]: "<<numCells[i]<<std::endl;
            
                arg.pos.clampCube(VectorND<Float,DIM>(),domainMax);
                std::cout<<"domainmax: "<<domainMax[i]<<std::endl;
                std::cout<<"ratio[i]: "<<int(std::floor(arg.pos[i]/boxWidth[i]))<<std::endl;
            }
            assert((arg.pos[i]/boxWidth[i])>=0);
            assert((arg.pos[i]/boxWidth[i])<numCells[i]);
            ivec[i]=std::floor(arg.pos[i]/boxWidth[i]);
        }
        return intvectorToIndex(ivec);
    }

    void rebuildGrid(){
        assert(productOfSizes>0);
        assert(plist!=nullptr);

        //std::cout<<"Product of Sizes: "<<productOfSizes<<std::endl;
        
        idarr=std::vector<std::vector<size_t> >(productOfSizes);
        for(size_t n=0;n<plist->size();n++) {
            (*plist)[n].pos.clampCube(VectorND<Float,DIM>(),domainMax);
            VectorND<int,DIM> ivec;
            for(int i=0;i<DIM;i++){
                assert(((*plist)[n].pos[i]/boxWidth[i])>=0);
                //std::cout<<((*plist)[n].pos[i]/boxWidth[i])<<std::endl;
                //std::cout<<domainSize[i]<<std::endl;
                assert(((*plist)[n].pos[i]/boxWidth[i])<numCells[i]);
                ivec[i]=std::floor((*plist)[n].pos[i]/boxWidth[i]);
            }
            //size_t indx=intvectorToIndex(ivec);
            //std::cout<<"Index of particle "<<n<<" is "<<indx<<std::endl;

            //idarr[].push_back(n);
            //in 2d: idarr[idy*numCellsX+idx].push_back(n);
            int idx=intvectorToIndex(ivec);
            assert(0<=idx && idx<idarr.size());
            idarr[idx].push_back(n);
        }
        needsRebuild=false;
    }



    //Do this after the iterators!
    //void timestepGrid(Float dt); 

    class PGridPairIt {
        PGrid<Float,DIM> *pg;
        VectorND<int,DIM> a;
        size_t i;
        VectorND<int,DIM> da;
        size_t j;
        size_t aind;
        size_t bind;

        bool isValid(){
            if(isEmpty())
                return false;
            if(aind==bind)
                return (i<j)&&(j<pg->idarr[aind].size());
            auto b=a+da;
            for(int k=0;k<DIM;k++){
                if(b[k]>=pg->numCells[k])
                    return false;
            }
            return (i<pg->idarr[aind].size())&&(j<pg->idarr[bind].size());
        }



        //step b index
        //make sure that 0<=a<numCells and 0<=a+da<numCells,
        //as well as aind=intvectorToIndex(a); bind=intvectorToIndex(a+da);
        //If we reach all the way to the end, increase aind by 1.
        void stepABOnce(){
            //if there are no particles at aind, no point in iterating over
            //2^DIM b index cells.
            if(pg->idarr[aind].size()==0){
                aind++;
                bind=aind;
                a=pg->indexToIntvector(aind);
                da=VectorND<int,DIM>();
                i=0;
                j=1;
                return;
            }
            int k=0;
            for(;k<DIM;k++) {
                da[k]++;
                if(da[k]==2)
                    da[k]=0;
                else 
                    break;
            }
            if(k==DIM){
                aind++;
                bind=aind;
                a=pg->indexToIntvector(aind);
                da=VectorND<int,DIM>();
                i=0;
                j=1;
            } else {
                i=0;
                j=0;
                bind=pg->intvectorToIndex(a+da);
            }
        }
        void stepABNextValid(){
            do {
                stepABOnce();
            } while(!isEmpty()&&!isValid());
        }
        
        //Each call to popFront either moves aind+=1 or aind+=0. 
        //It leaves 0<=a<numCells and 0<=a+da<numCells valid, 
        //as well as aind=intvectorToIndex(a); bind=intvectorToIndex(a+da);
        //however it might be that !(0<=i<idarr[aind].size()) or !(0<=j<idarr[bind].size())
        //so you still need to call isValid()
        void pop_front(){
            if(isEmpty())
                return;
            assert(!(pg->needsRebuild));

            if(aind==bind){
                i++;
                if(i>=j){
                    i=0;
                    j++;
                }
                if(j>=pg->idarr[aind].size()){
                    //We have exhausted all aind <-> aind pairs.
                    stepABNextValid();
                }
            } else {
                i++;
                if(i>=pg->idarr[aind].size()){
                    i=0;
                    j++;
                    if(j>=pg->idarr[bind].size()){
                        //We have exhausted all aind <-> bind pairs.
                        j=0;
                        stepABNextValid();
                    }
                }
            }
            if(isValid()){
                first=&((*(pg->plist))[(pg->idarr[aind])[i]]);
                second=&((*(pg->plist))[(pg->idarr[bind])[j]]);
            } else {
                first=nullptr;
                second=nullptr;
            }
        }
        public:
        Particle<Float,DIM> *first;
        Particle<Float,DIM> *second;
        explicit PGridPairIt(PGrid<Float,DIM> *pg,size_t aind=0) : 
            pg(pg), a(), i(0), da(), j(1), aind(aind), bind(aind), first(nullptr),second(nullptr) { 
            a=pg->indexToIntvector(aind);
            while(!isEmpty()&&!isValid()){
                pop_front();
            }
            if(isValid()){
                first=&((*(pg->plist))[(pg->idarr[aind])[i]]);
                second=&((*(pg->plist))[(pg->idarr[bind])[j]]);
            }
        }
        bool isEmpty(){
            return aind>=pg->productOfSizes;
        }
        PGridPairIt& operator++() {
            pop_front();
            return *this;
        }
        std::pair<Particle<Float,DIM>*,Particle<Float,DIM>*> operator*() const {
            return std::make_pair(first,second);
        }
        bool operator==(const PGridPairIt &other){
            assert(pg==other.pg);
            return (aind>=pg->idarr.size() && other.aind>=pg->idarr.size())||(
                    aind==other.aind&&bind==other.bind&&i==other.i&&j==other.j);
        }
        bool operator!=(const PGridPairIt &other){
            return !(operator==(other));
        }
    };

    class PGridPairRange{
        PGrid<Float,DIM> *pg;
        public:
        PGridPairRange(PGrid<Float,DIM> *pg) : pg(pg){ }

        PGridPairIt begin() { 
            return PGridPairIt(pg,0);
        }
        PGridPairIt end() { 
            return PGridPairIt(pg,pg->plist->size());
        }
    };

    class PGridUpdateIt {
        PGrid<Float,DIM> *pg;
        size_t pind; // index into pg->plist
        size_t aind; // index into pg->idarr
        void updateAIndex(){
            if(pg!=nullptr&&pind<pg->plist->size()){
                aind=pg->getParticleIndex((*pg->plist)[pind]);
                part=&((*pg->plist)[pind]);
            } else {
                aind=0;
                part=nullptr;
            }
        }
    public:
        Particle<Float,DIM> *part;
        PGridUpdateIt(PGrid<Float,DIM> *pg,size_t pind) : pg(pg),pind(pind),aind(0),part(nullptr) { 
            updateAIndex();
        }

        PGridUpdateIt& operator++() {
            (*pg->plist)[pind].pos.clampCube(VectorND<Float,DIM>(),pg->domainMax);
            size_t aind_new=pg->getParticleIndex((*pg->plist)[pind]);
            if(aind!=aind_new){
                pg->particleChangeGrid(pind,aind,aind_new);
            }
            pind++;
            updateAIndex();
            return *this;
        }
        Particle<Float,DIM>* operator*() const {
            return part;
        }
        bool operator!=(const PGridUpdateIt &other) const {
            assert(pg==other.pg);
            return (pind!=other.pind);
        }


    };
    class PGridUpdateRange {
        PGrid<Float,DIM> *pg;
    public:
        PGridUpdateRange(PGrid<Float,DIM> *pg) : pg(pg){ }

        PGridUpdateIt begin() { 
            return PGridUpdateIt(pg,0);
        }
        PGridUpdateIt end() { 
            return PGridUpdateIt(pg,pg->plist->size());
        }
    };

    PGridPairRange pairs(){
        return PGridPairRange(this);
    }
    PGridUpdateRange updateLoop(){
        return PGridUpdateRange(this);
    }

};
#endif
