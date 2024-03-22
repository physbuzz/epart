#ifndef PGRID_H
#define PGRID_H
#include <utility>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "phystructs.h"

template<typename Float,int DIM>
class PGrid {
public:
    static_assert(0<DIM,"PGrid DIM must be positive.");

    std::vector<std::vector<size_t> > idarr;

    //integer number of cells in each dimension
    VectorND<int,DIM> numCells;

    //Box dimensions of each cell
    VectorND<Float,DIM> boxWidth;

    //Upper bounds on the particle locations
    VectorND<Float,DIM> domainSize;

    bool needsRebuild;
    size_t productOfSizes;

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
        needsRebuild(true), productOfSizes(1), plist(nullptr) { }

    PGrid(std::vector<Particle<Float,DIM> > *plist,VectorND<Float,DIM> domainSize,Float maxH) : 
        idarr(), numCells(), boxWidth(), domainSize(), 
        needsRebuild(true), productOfSizes(0), 
        plist(plist) {
        setParams(domainSize,maxH);
    }

    void setParams(VectorND<Float,DIM> domainSize,Float maxH){
        assert(maxH>0 &&domainSize.min()>Float(0));
        this->productOfSizes=1;
        this->domainSize=domainSize;
        for(size_t i=0;i<DIM;i++){
            this->numCells[i]=std::max(int(domainSize[i]/maxH),1);
            this->boxWidth[i]=domainSize[i]/this->numCells[i];
        }
        this->productOfSizes=this->numCells.product();
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
    VectorND<int,DIM> positionToIntvecClamped(const VectorND<Float,DIM> &arg){
        VectorND<int,DIM> ivec;
        for(int i=0;i<DIM;i++){
            ivec[i]=std::clamp(int(arg[i]/boxWidth[i]),0,numCells[i]-1);
        }
        return ivec;
    }

    size_t getParticleIndex(Particle<Float,DIM> &arg){
        VectorND<int,DIM> ivec=positionToIntvecClamped(arg.pos);
        return intvectorToIndex(ivec);
    }

    void rebuildGrid(){
        assert(productOfSizes>0);
        assert(plist!=nullptr);
        idarr=std::vector<std::vector<size_t> >(productOfSizes);
        for(size_t pindx=0;pindx<plist->size();pindx++) {
            (*plist)[pindx].pos.clampCube(VectorND<Float,DIM>(),domainSize);
            VectorND<int,DIM> ivec=positionToIntvecClamped((*plist)[pindx].pos);
            int idx=intvectorToIndex(ivec);
            idarr.at(idx).push_back(pindx);
        }
        needsRebuild=false;
    }



    //Do this after the iterators!
    //void timestepGrid(Float dt); 
    /*
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

    PGridPairRange pairs(){
        return PGridPairRange(this);
    }
    
    */

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
            (*pg->plist)[pind].pos.clampCube(VectorND<Float,DIM>(),pg->domainSize);
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
    PGridUpdateRange updateLoop(){
        return PGridUpdateRange(this);
    }




    class PGridRectIt {
        PGrid<Float,DIM> *pg;

        //Lowest coordinate in the search rectangle (inclusive)
        VectorND<int,DIM> bl;
        //Highest coordinate in the search rectangle (exclusive)
        VectorND<int,DIM> tr;
        //We look at all grid coordinates with bl[k]<=loc[k]<tr[k].
        VectorND<int,DIM> loc;

        // index into pg->idarr
        size_t aind; 
        // index into pg->plist
        size_t pind;


        void stepAind(){
            //Increase loc[0]. If we reach the bound tr[0],
            //set loc[0]=bl[0] and increase loc[1], etc.
            for(int k=0;k<DIM;k++) {
                loc[k]++;
                if(loc[k]>=tr[k])
                    loc[k]=bl[k];
                else {
                    aind=pg->intvectorToIndex(loc);
                    return;
                }
            }
            //If we reach here then we've just set loc[DIM-1]=bl[DIM-1],
            //but we want the iterator to be empty now.
            loc=tr;
            aind=pg->intvectorToIndex(loc);
        }
        bool isEmpty(){
            return loc==tr;
        }
        bool isValid(){
            return (!isEmpty())&&pind<pg->idarr[aind].size();
        }

        void pop_front(){
            if(isEmpty())
                return;
            assert(!(pg->needsRebuild));
            assert(aind>=0 && aind<pg->idarr.size());

            //If there are particles left in the current cell, get the next one.
            pind++;
            if(pind<pg->idarr[aind].size()){
                partptr=&((*(pg->plist))[(pg->idarr[aind])[pind]]);
                return;
            }

            //Else find the next nonempty cell, if there is one.
            pind=0;
            stepAind();
            while(!isEmpty() && pg->idarr[aind].size()==0){
                stepAind();
            }
            if(isEmpty()){
                partptr=nullptr;
            } else {
                partptr=&((*(pg->plist))[(pg->idarr[aind])[pind]]);
            }
        }
        public:

        Particle<Float,DIM> *partptr;

        explicit PGridRectIt(PGrid<Float,DIM> *pg,VectorND<int,DIM> bl,VectorND<int,DIM> tr) : 
            pg(pg), bl(bl),tr(tr), loc(bl), aind(0),pind(0), partptr(nullptr){ 

            for(int k=0;k<DIM;k++){
                //Iterator expects that the passed in coordinates are inside the simulation domain
                //assert(bl[k]>=0 && bl[k]<pg->numCells[k]);
                //assert(tr[k]>=1 && tr[k]<=pg->numCells[k]);
            }
            aind=pg->intvectorToIndex(loc);
            if(isEmpty())
                return;
            if(isValid())
                partptr=&((*(pg->plist))[(pg->idarr[aind])[pind]]);
            else
                pop_front();
        }
        //Used to construct the end iterator with PGridRectIt(pg,bl,tr,tr);
        explicit PGridRectIt(PGrid<Float,DIM> *pg,VectorND<int,DIM> bl,VectorND<int,DIM> tr, VectorND<int,DIM> loc) : 
            pg(pg), bl(bl),tr(tr), loc(loc), aind(0),pind(0), partptr(nullptr){ 

            for(int k=0;k<DIM;k++){
                //Iterator expects that the passed in coordinates are inside the simulation domain
                //assert(bl[k]>=0 && tr[k]<=pg->numCells[k]);
            }
            aind=pg->intvectorToIndex(loc);
            pop_front();
        }
        PGridRectIt& operator++() {
            pop_front();
            return *this;
        }
        Particle<Float,DIM>* operator*() const {
            return partptr;
        }
        bool operator!=(const PGridRectIt &other){
            return !(this->isEmpty());
        }
    };
    class PGridNearbyRange {
        PGrid<Float,DIM> *pg;
        VectorND<Float,DIM> pos;
        Float r;
        VectorND<int,DIM> bl;
        VectorND<int,DIM> tr;
    public:
        PGridNearbyRange(PGrid<Float,DIM> *pg, VectorND<Float,DIM> pos, Float r) : 
            pg(pg), pos(pos), r(r),bl(),tr(){ 

            VectorND<Float,DIM> blcoord;
            VectorND<Float,DIM> trcoord;
            for(int k=0;k<DIM;k++){
                blcoord[k]=pos[k]-r;
                trcoord[k]=pos[k]+r;
            }
            bl=pg->positionToIntvec(blcoord);
            tr=pg->positionToIntvec(trcoord);
            for(int k=0;k<DIM;k++){
                tr[k]=tr[k]+1;
            }
            bl.clampCube(VectorND<int,DIM>(),pg->numCells);
            tr.clampCube(VectorND<int,DIM>(),pg->numCells);

            bool empty=false;
            for(int k=0;k<DIM;k++){
                if(bl[k]==tr[k])
                    empty=true;
            }
            if(empty)
                bl=tr;
        }
        PGridRectIt begin() { 
            return PGridRectIt(pg,bl,tr);
        }
        PGridRectIt end() { 
            return PGridRectIt(pg,bl,tr,tr);
        }
    };
    PGridNearbyRange nearbyLoop(VectorND<Float,DIM> pos, Float r){
        return PGridNearbyRange(this, pos, r);
    }
};
#endif
