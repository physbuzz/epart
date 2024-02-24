CPP=g++


epart: epart.cpp PGrid.h VectorND.h phystructs.h ParticleList.h ImageUtil.h
	$(CPP) -g epart.cpp -o epart

run: epart
	./epart

epart-fast: epart.cpp PGrid.h VectorND.h phystructs.h ParticleList.h ImageUtil.h
	$(CPP) -O3 epart.cpp -o epart-fast

#circlepack: circlepack.cpp PGrid.h VectorND.h phystructs.h ParticleList.h
#	$(CPP) -g circlepack.cpp -o circlepack
#
#run: circlepack
#	./circlepack
#
#
#tidy:
#	clang-tidy circlepack.cpp --
#
#clean: 
#	-rm ./circlepack
#	-rm ./fast
#
#FORCE: 
