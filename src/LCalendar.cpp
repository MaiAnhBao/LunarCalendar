/*
 * LCalendar.cpp
 *
 *  Created on: 30 déc. 2016
 *      Author: hnnguyen
 */

#include "LCalendar.h"

/**
 Description: Convert to Julius Date from Date
 @param the date
 @return the julius date
 */
int LCalendar::convertFromDateToJuliusDate(int year, int month, int day) {
#ifdef LOG
	std::cout << "LOG " << "convertFromDateToJuliusDate <<<<<<" << std::endl;
#endif
#ifdef DEBUG
	std::cout << "INPUT: " << " year " << year << " month " << month << " day " << day << std::endl;
#endif
	if ((day > 31) or (month > 12)) {
		return -1;
	}
	int a = (14 - month) / 12;
	int y = year + 4800 - a;
	int m = month + 12 * a - 3;
#ifdef DEBUG
	std::cout << "a = " << a << " y = " << y << " m = " << m << std::endl;
#endif
	int juliusDate = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100
			+ y / 400 - 32045;
	if (juliusDate < 2299161) {
		juliusDate = day + (153 * m + 2) / 5 + 365 * y + y / 4 - 32083;
	}
#ifdef LOG
	std::cout << "LOG " << "convertFromDateToJuliusDate >>>>>>" << std::endl;
#endif
	return juliusDate;
}

/**
 Description: Convert to Date from Julius Date
 @param the julius date
 @return the date
 */
std::tm LCalendar::convertFromJuliusDateToDate(int juliusDate) {
#ifdef LOG
	std::cout << "LOG " << "convertFromJuliusDateToDate >>>>>>" << std::endl;
#endif
	int a, b, c, d, e, m, day, month, year;
	if (juliusDate > dayChangeFromJulianToGregorian) {
		a = juliusDate + 32044;
		b = (a * 4 + 3) / 146097;
		c = a - (b * 146097) / 4;
#ifdef DEBUG
		std::cout << "\033[1;35m" << "DEBUG: " << "a = " << a << "b = " << b << "c = " << c << "\033[0m" << std::endl;
#endif
	} else {
		b = 0;
		c = juliusDate + 32082;
	}
	d = (4 * c + 3) / 1461;
	e = c - 1461 * d / 4;
	m = (5 * e + 2) / 153;
#ifdef DEBUG
	std::cout << "\033[1;35m" << "DEBUG: " << "d = " << d << " e = " << e << " m = " << m << "\033[0m\n" << std::endl;
#endif

	day = e - (153 * m + 2) / 5 + 1;
	month = m + 3 - 12 * (m / 10);
	year = b * 100 + d - 4800 + m / 10;
#ifdef DEBUG
	std::cout << "\033[1;35m" << "DEBUG: " << "day = " << day << " month = " << month << " year = " << year << "\033[0m\n" << std::endl;
#endif

	std::tm retVal;
	retVal.tm_year = year - 1900;
	retVal.tm_mon = month - 1;
	retVal.tm_mday = day;
#ifdef LOG
	std::cout << "LOG " << "convertFromJuliusDateToDate <<<<<<" << std::endl;
#endif
	return retVal;
}

/**
 Description: get the time of new moon k from 1/1/1900
 @param k-th of new moon,
 time zone
 @return the julius day of new moon
 */
int LCalendar::getNewMoonDay(int k, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "getNewMoonDay >>>>>>" << std::endl;
#endif
	int retVal;
	double T, T2, T3, dr, Jd1, M, Mpr, F, C1, deltat, JdNew;
	T = k / 1236.85;
	T2 = T * T;
	T3 = T2 * T;
	dr = PI / 180;
	Jd1 = 2415020.75933 + 29.53058868 * k + 0.0001178 * T2 - 0.000000155 * T3;
	Jd1 = Jd1 + 0.00033 * std::sin((166.56 + 132.87 * T - 0.009173 * T2) * dr); // Mean new moon
	M = 359.2242 + 29.10535608 * k - 0.0000333 * T2 - 0.00000347 * T3; //Sun's mean anomaly
	Mpr = 306.0253 + 385.81691806 * k + 0.0107306 * T2 + 0.00001236 * T3;
	F = 21.2964 + 390.67050646 * k - 0.0016528 * T2 - 0.00000239 * T3;
	C1 = (0.1734 - 0.000393 * T) * std::sin(M * dr)
			+ 0.0021 * std::sin(2 * dr * M);
	C1 = C1 - 0.4068 * std::sin(Mpr * dr) + 0.0161 * std::sin(dr * 2 * Mpr);
	C1 = C1 - 0.0004 * std::sin(dr * 3 * Mpr);
	C1 = C1 + 0.0104 * std::sin(dr * 2 * F) - 0.0051 * std::sin(dr * (M + Mpr));
	C1 = C1 - 0.0074 * std::sin(dr * (M - Mpr))
			+ 0.0004 * std::sin(dr * (2 * F + M));
	C1 = C1 - 0.0004 * std::sin(dr * (2 * F - M))
			- 0.0006 * std::sin(dr * (2 * F + Mpr));
	C1 = C1 + 0.0010 * std::sin(dr * (2 * F - Mpr))
			+ 0.0005 * std::sin(dr * (2 * Mpr + M));
	if (T < -11) {
		deltat = 0.001 + 0.000839 * T + 0.0002261 * T2 - 0.00000845 * T3
				- 0.000000081 * T * T3;
	} else {
		deltat = -0.000278 + 0.000265 * T + 0.000262 * T2;
	};
	JdNew = Jd1 + C1 - deltat;
	retVal = static_cast<int>(JdNew + 0.5 + timeZone / 24);
#ifdef DEBUG
	std::cout << "DEBUG: " << "new moon day " << retVal << std::endl;
#endif
#ifdef LOG
	std::cout << "LOG " << "getNewMoonDay <<<<<<" << std::endl;
#endif
	return retVal;
}

/**
 Description: calculate the Sun's longtitude
 @param: the julius date
 time zone
 @return the Sun's longtitude of the corresponding date
 */
int LCalendar::getSunLongtitude(int jdn, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "getSunLongtitude >>>>>>" << std::endl;
#endif
	float T, T2, dr, M, L0, DL, L;
	T = (jdn - 2451545.5 - timeZone / 24) / 36525; // Time in Julian centuries from 2000-01-01 12:00:00 GMT
	T2 = T * T;
	dr = PI / 180; // degree to radian
	M = 357.52910 + 35999.05030 * T - 0.0001559 * T2 - 0.00000048 * T * T2; // mean anomaly, degree
	L0 = 280.46645 + 36000.76983 * T + 0.0003032 * T2; // mean longitude, degree
	DL = (1.914600 - 0.004817 * T - 0.000014 * T2) * std::sin(dr * M);
	DL = DL + (0.019993 - 0.000101 * T) * std::sin(dr * 2 * M)
			+ 0.000290 * std::sin(dr * 3 * M);
	L = L0 + DL; // true longitude, degree
#ifdef DEBUG
	std::cout << "DEBUG: " << "sun longtitude degree = "  << L << std::endl;
#endif
	L = L * dr;
	L = L - PI * 2 * (static_cast<int>(L / (PI * 2))); // Normalize to (0, 2*PI)

#ifdef DEBUG
					std::cout << "DEBUG: " << "sun longtitude = " << static_cast<int>(L / PI * 6) << " degree "  << L << std::endl;
#endif

#ifdef LOG
	std::cout << "LOG " << "getSunLongtitude <<<<<<" << std::endl;
#endif
	return static_cast<int>(L / PI * 6);
}

/**
 Description: get the new moon of the november of year
 @param: year, time zone
 @return: the new moon of the november of the year
 */
int LCalendar::getLunarMonthElevent(int year, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "getLunarMonthElevent >>>>>>" << std::endl;
#endif
#ifdef DEBUG
	std::cout << "DEBUG: " << "year " << year << " time zone: " << timeZone << std::endl;
#endif
	int k, off, nm, sunLong;
	off = convertFromDateToJuliusDate(year, 12, 31) - 2415021;

#ifdef DEBUG
	std::cout << "DEBUG: " << "different from " << "2415021 " << "to 31/12/" << year << ": " << off << std::endl;
#endif

	k = static_cast<int>(off / 29.530588853);

#ifdef DEBUG
	std::cout << "DEBUG: " << "the order of new moon day " << k << std::endl;
#endif

	nm = getNewMoonDay(k, timeZone);

#ifdef DEBUG
	std::cout << "DEBUG: " << "new monday of November " << nm << std::endl;
#endif

	sunLong = getSunLongtitude(nm, timeZone); // sun longitude at local midnight

#ifdef DEBUG
	std::cout << "DEBUG: " << "Sun longtitude of November " << sunLong << std::endl;
#endif

	if (sunLong >= 9) {
#ifdef DEBUG
		std::cout << "DEBUG: " << "sun longtitude > 9" << std::endl;
#endif
		nm = getNewMoonDay(k - 1, timeZone);
	}

#ifdef LOG
	std::cout << "LOG " << "getLunarMonthElevent <<<<<<" << std::endl;
#endif
	return nm;
}

/**
 Description: Get the Leap month
 @param: the julius date of the november, time zone
 @return: the leap month
 */
int LCalendar::getLeapMonthOffset(int a11, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "getLeapMonthOffset >>>>>>" << std::endl;
#endif
	int k, last, arc, i;
	k = static_cast<int>((a11 - 2415021.076998695) / 29.530588853 + 0.5);
	last = 0;
	i = 1; // We start with the month following lunar month 11
	arc = getSunLongtitude(getNewMoonDay(k + i, timeZone), timeZone);
	do {
		last = arc;
		i++;
		arc = getSunLongtitude(getNewMoonDay(k + i, timeZone), timeZone);
	} while (arc != last && i < 14);
#ifdef LOG
	std::cout << "LOG " << "getLeapMonthOffset <<<<<<" << std::endl;
#endif
	return i - 1;
}

/**
 Description: convert from solar date to lunar date
 @param: day, month, year, time zone
 @return: none
 */
void LCalendar::convertSolarToLunar(int dd, int mm, int yy, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "convertSolarToLunar >>>>>>" << std::endl;
#endif
	int k, dayNumber, monthStart, a11, b11, lunarDay, lunarMonth, lunarYear,
			lunarLeap, diff, nextMonthStart;
	int leapMonthDiff;
	dayNumber = convertFromDateToJuliusDate(yy, mm, dd);
#ifdef DEBUG
	std::cout << "DEBUG: " << " julius day  = " << dayNumber << std::endl;
#endif
	k = static_cast<int>((dayNumber - 2415021.076998695) / 29.530588853);
#ifdef DEBUG
	std::cout << "DEBUG: " << "the k-th new moon from 1/1/1900 " << k << std::endl;
#endif
	monthStart = getNewMoonDay(k + 1, timeZone);

#ifdef DEBUG
	std::cout << "DEBUG: " << "month Start = " << monthStart << std::endl;
#endif

	if (monthStart > dayNumber) {
		nextMonthStart = monthStart;
		monthStart = getNewMoonDay(k, timeZone);
	} else {
		nextMonthStart = getNewMoonDay(k+2,timeZone);
	}

	//check full month
	if (nextMonthStart - monthStart < 30) {
		fullMonth = false;
	}

#ifdef DEBUG
	std::cout << "DEBUG: " << "Month start after compare = " << monthStart
	<< std::endl;
#endif
	a11 = getLunarMonthElevent(yy, timeZone);
#ifdef DEBUG
	std::cout << "DEBUG: " << "November of this year" << a11 << std::endl;
#endif
	b11 = a11;
	if (a11 >= monthStart) {
		lunarYear = yy;
		a11 = getLunarMonthElevent(yy - 1, timeZone);
	} else {
		lunarYear = yy + 1;
		b11 = getLunarMonthElevent(yy + 1, timeZone);
	}
#ifdef DEBUG
	std::cout << "DEBUG: " << "November after comparing " << a11 << " " << b11 << std::endl;
#endif
	lunarDay = dayNumber - monthStart + 1;
#ifdef DEBUG
	std::cout << "DEBUG: " << "Lunar days " << lunarDay << std::endl;
#endif

	diff = static_cast<int>((monthStart - a11) / 29);
#ifdef DEBUG
	std::cout << "DEBUG: " << "different between month start and the last november " << diff << std::endl;
#endif

	lunarLeap = 0;
	lunarMonth = diff + 11;

	if (b11 - a11 > 365) {
#ifdef DEBUG
		std::cout << "DEBUG: " << " there is a leap month" << std::endl;
#endif
		leapMonthDiff = getLeapMonthOffset(a11, timeZone);
		if (diff >= leapMonthDiff) {
			lunarMonth = diff + 10;
			if (diff == leapMonthDiff) {
				lunarLeap = 1;
			}
		}
	}
	if (lunarMonth > 12) {
		lunarMonth = lunarMonth - 12;
	}
	if (lunarMonth >= 11 && diff < 4) {
		lunarYear -= 1;
	}
#ifdef DEBUG
	std::cout << "\033[1;35m" << "DEBUG: " << lunarDay << " " << lunarMonth << " "
			<< lunarYear << "\033[0m\n" << std::endl;
#endif


	mLunarDate.mLunarDay = lunarDay;
	mLunarDate.mLunarMonth = lunarMonth;
	mLunarDate.mLunarYear = lunarYear;

#ifdef LOG
	std::cout << "LOG " << "convertSolarToLunar <<<<<<" << std::endl;
#endif
}

/**
 Description: convert from lunar date to the solar date
 @param: day, month, year of lunar date, lunar leap, time zone
 @return: day, month, year of solar date
 */
std::tm LCalendar::convertLunarToSolar(int lunarDay, int lunarMonth,
		int lunarYear, int lunarLeap, int timeZone) {
#ifdef LOG
	std::cout << "LOG " << "convertLunarToSolar >>>>>>" << std::endl;
#endif
	std::tm retVal;
	int k, a11, b11, off, leapOff, leapMonth, monthStart;
	if (lunarMonth < 11) {
		a11 = getLunarMonthElevent(lunarYear - 1, timeZone);
		b11 = getLunarMonthElevent(lunarYear, timeZone);
	} else {
		a11 = getLunarMonthElevent(lunarYear, timeZone);
		b11 = getLunarMonthElevent(lunarYear + 1, timeZone);
	}
	off = lunarMonth - 11;
	if (off < 0) {
		off += 12;
	}
	if (b11 - a11 > 365) {
		leapOff = getLeapMonthOffset(a11, timeZone);
		leapMonth = leapOff - 2;
		if (leapMonth < 0) {
			leapMonth += 12;
		}
		if (lunarLeap != 0 && lunarMonth != leapMonth) {
			return retVal;
		} else if (lunarLeap != 0 || off >= leapOff) {
			off += 1;
		}
	}
	k = static_cast<int>(0.5 + (a11 - 2415021.076998695) / 29.530588853);
	monthStart = getNewMoonDay(k + off, timeZone);
#ifdef DEBUGSOLAR
	std::cout << "\033[1;35m" << "DEBUG: " << "month start: " << monthStart << "\033[0m\n" << std::endl;
#endif
#ifdef LOG
	std::cout << "LOG " << "convertLunarToSolar <<<<<<" << std::endl;
#endif
	return convertFromJuliusDateToDate(monthStart + lunarDay - 1);
}

/**
 * @Description: set Time Zone
 * @param: float Time Zone
 * @return none
 */
void LCalendar::setTimeZone(float fTimeZone) {
	mTimeZone = fTimeZone;
}

/**
 * @Description: get Time Zone
 * @param: Void
 * @return Time Zone
 */
float LCalendar::getTimeZone() {
	return mTimeZone;
}

/**
 * @Description: get Lunar date
 * @param Void
 * @return: lunar date
 */
LunarDate LCalendar::getLunarDate() {
	return mLunarDate;
}

/**
 * @Description: constructor default with time zone
 * @param: time zone
 * @return: lunar date of today
 */
LCalendar::LCalendar(float fTimeZone) :
		mTimeZone(fTimeZone) {

	std::chrono::time_point<std::chrono::system_clock> currentTime;
	currentTime = std::chrono::system_clock::now();
	std::time_t current_time = std::chrono::system_clock::to_time_t(
			currentTime);

	std::tm local_time = *localtime(&current_time);
	int year = local_time.tm_year + 1900;
	int month = local_time.tm_mon + 1;
	int day = local_time.tm_mday;

	convertSolarToLunar(day, month, year, mTimeZone);
}

/**
 * @Description: constructor with the date and time zone
 * @param: time zone, date
 * @return: lunar date of the corresponding date
 */
LCalendar::LCalendar(std::tm iDate, float fTimeZone) :
		mTimeZone(fTimeZone) {
	int year = iDate.tm_year + 1900;
	int month = iDate.tm_mon + 1;
	int day = iDate.tm_mday;
	convertSolarToLunar(day, month, year, mTimeZone);

}

/**
 * print lunar date
 * @param none
 * @return none
 */
void LCalendar::toString() {
	if (fullMonth) {
		std::cout << "Full Month" << std::endl;
	} else {
		std::cout << "Not Full Month" << std::endl;
	}
	std::cout << "Day " << mLunarDate.mLunarDay << " Month " << mLunarDate.mLunarMonth << " Year " << mLunarDate.mLunarYear << std::endl;
}
