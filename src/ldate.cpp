/*
 * ldate.cpp
 *
 *  Created on: 30 d�c. 2016
 *      Author: hnnguyen
 */


#include "LCalendar.h"

int main(int argc, char **argv) {
	LCalendar *objLunarCalendar = new LCalendar(7);
	objLunarCalendar->toString();
	return 0;
}
