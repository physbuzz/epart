CPP=g++


epart: epart.cpp PGrid.h VectorND.h phystructs.h ParticleList.h ImageUtil.h
	$(CPP) -g epart.cpp -o epart

run: epart
	./epart
#circlepack: circlepack.cpp PGrid.h VectorND.h phystructs.h ParticleList.h
#	$(CPP) -g circlepack.cpp -o circlepack
#
#run: circlepack
#	./circlepack
#
#fast: circlepack.cpp PGrid.h VectorND.h phystructs.h ParticleList.h
#	$(CPP) -DNDEBUG -O3 circlepack.cpp -o fast
#	./fast
#
#tidy:
#	clang-tidy circlepack.cpp --
#
#clean: 
#	-rm ./circlepack
#	-rm ./fast
#
#FORCE: 
