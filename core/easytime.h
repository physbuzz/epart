#ifndef EASYTIME_H
#define EASYTIME_H

#include <iostream>
#include <ctime>

/* Example usage:
EasyTimer x=EasyTimer();
x.tick();

some_operation_that_takes_a_ton_of_time();

std::cout<<"Time Elapsed Since Last Tick: "<<x.tick()<<std::endl;

some_operation_that_takes_a_ton_of_time();

std::cout<<"Time Elapsed Since Last Tick: "<<x.tick()<<std::endl;

*/
class EasyTimer {
	clock_t timer;
public:
	EasyTimer();
	double tick();
};

inline 
EasyTimer::EasyTimer(){
	timer=std::clock();
}

inline
double EasyTimer::tick(){
	clock_t timer2=clock();
	double duration=double(timer2-timer)/CLOCKS_PER_SEC;
	timer=timer2;
	return duration; 
}

#endif