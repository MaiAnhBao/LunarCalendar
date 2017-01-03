/*
 * LCalendar.h
 *
 *  Created on: 30 déc. 2016
 *      Author: hnnguyen
 */

#ifndef LUNARCALENDAR_LCALENDAR_H_
#define LUNARCALENDAR_LCALENDAR_H_

#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include <cmath>

struct LunarDate {
	int mLunarDay;
	int mLunarMonth;
	int mLunarYear;
};

class LCalendar {
private:
	const double PI = 3.14159265359;
	const int dayChangeFromJulianToGregorian = 2299160;
	const int dayOfYear = 365;
	const float fJuliusDateFirstJanuary1900 = 2415021.076998695;
	const int  iJuliusDateFirstJanuary1900 = 2415021;
	const float fSynodicMonth = 29.530588853;
	const int iJulianDayOfCentury = 36525;

	const std::vector<std::string> vecCan{""};
	const std::vector<std::string> vecChi{""};
	const std::vector<std::string> vecSolarTerms{""};

	LCalendar(const LCalendar&) = delete;
	LCalendar &operator=(const LCalendar&) = delete;

	float mTimeZone;
	LunarDate mLunarDate;
	bool fullMonth = true;


	float getTimeZone();

	int convertFromDateToJuliusDate(int year, int month, int day);
	std::tm convertFromJuliusDateToDate(int juliusDate);
	int getNewMoonDay(int k, int timeZone);
	int getSunLongtitude(int jdn, int timeZone);
	int getLunarMonthElevent(int year, int timeZone);
	int getLeapMonthOffset(int a11, int timeZone);
	void convertSolarToLunar(int dd, int mm, int yy, int timeZone);
	std::tm convertLunarToSolar(int lunarDay, int lunarMonth, int lunarYear,
			int lunarLeap, int timeZone);
public:
//	LCalendar() = default;
	LCalendar(float fTimeZone = 0.0);
	LCalendar(std::tm iDate, float fTimeZone = 0.0);
	~LCalendar() = default;
	void setTimeZone(float fTimeZone);
	LunarDate getLunarDate();
	void toString();
};

#endif /* LUNARCALENDAR_LCALENDAR_H_ */
