import datetime as dt
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None 
class calendar_convert:
    def __init__(self,date_input,time_zone=7):
        self.date_input = date_input.reset_index(drop = True)
        self.time_zone = time_zone
        self.bgn = dt.date(1800,1,1)

    #Number of days from 1st Julian day
    def julian_day(self, dd):
        x = pd.to_datetime(dd).dt.date
        gre = (x - self.bgn)//np.timedelta64(1,'D') + 2378497
        return gre

    #Convert that back into Gregorian date
    def julian_day_to_date(self,gre):
        a = gre - 2378497
        date = [self.bgn + dt.timedelta(days=d) for d in a]
        return date
    
    def new_moon(self, k_th):
        # Time in Julian centuries from 1900 January 0.5
        time_julian = k_th / 1236.85
        time_julian_2 = time_julian * time_julian
        time_julian_3 = time_julian_2 * time_julian
        degree_to_radian = np.pi / 180
        julian_day_1 = (2415020.75933 + 29.53058868 * k_th + 0.0001178 * time_julian_2 - 0.000000155 * time_julian_3)
        julian_day_1 = (julian_day_1 + 0.00033*np.sin((166.56 + 132.87*time_julian - 0.009173 * time_julian_2) *
                                         degree_to_radian))
        # Mean new moon
        mean_new_moon = (359.2242 + 29.10535608*k_th -
                         0.0000333*time_julian_2 - 0.00000347*time_julian_3)
        # Sun's mean anomaly
        sun_mean_anomaly = (306.0253 + 385.81691806*k_th +
                            0.0107306*time_julian_2 + 0.00001236*time_julian_3)
        # Moon's mean anomaly
        moon_mean_anomaly = (21.2964 + 390.67050646*k_th -
                             0.0016528*time_julian_2 - 0.00000239*time_julian_3)
        # Moon's argument of latitude
        moon_arg_lat = ((0.1734 - 0.000393*time_julian) *
                        np.sin(mean_new_moon*degree_to_radian) +
                        0.0021*np.sin(2*degree_to_radian*mean_new_moon))
        moon_arg_lat = (moon_arg_lat -
                        0.4068*np.sin(sun_mean_anomaly*degree_to_radian)
                        + 0.0161*np.sin(degree_to_radian*2*sun_mean_anomaly))
        moon_arg_lat = (moon_arg_lat -
                        0.0004*np.sin(degree_to_radian*3*sun_mean_anomaly))
        moon_arg_lat = (moon_arg_lat +
                        0.0104*np.sin(degree_to_radian*2*moon_mean_anomaly)
                        - 0.0051 * np.sin(degree_to_radian *
                                            (mean_new_moon + sun_mean_anomaly)))
        moon_arg_lat = (moon_arg_lat -
                        0.0074*np.sin(degree_to_radian *
                                        (mean_new_moon - sun_mean_anomaly))
                        + 0.0004*np.sin(degree_to_radian *
                                          (2*moon_mean_anomaly + mean_new_moon)))
        moon_arg_lat = (moon_arg_lat - 0.0004*np.sin(degree_to_radian *
                                                       (2*moon_mean_anomaly -
                                                        mean_new_moon))
                        - 0.0006 * np.sin(degree_to_radian *
                                            (2*moon_mean_anomaly
                                             + sun_mean_anomaly)))
        moon_arg_lat = (moon_arg_lat + 0.0010*np.sin(degree_to_radian *
                                                       (2*moon_mean_anomaly -
                                                        sun_mean_anomaly))
                        + 0.0005*np.sin(degree_to_radian *
                                          (2*sun_mean_anomaly + mean_new_moon))
                        )

        if k_th.shape==():
            if time_julian < -11:
                deltat = (0.001 + 0.000839*time_julian + 0.0002261*time_julian_2
                      - 0.00000845*time_julian_3 - 0.000000081*time_julian*time_julian_3)
            else:
                deltat = -0.000278 + 0.000265*time_julian + 0.000262*time_julian_2
        else:
            deltat = -0.000278 + 0.000265*time_julian + 0.000262*time_julian_2

            deltat[time_julian < -11] = (0.001 + 0.000839*time_julian + 0.0002261*time_julian_2
                          - 0.00000845*time_julian_3 -  0.000000081*time_julian*time_julian_3)

        new_julian_day = julian_day_1 + moon_arg_lat - deltat  

        return np.floor(new_julian_day + 0.5 + self.time_zone / 24.)
    
    def sun_longitude(self, dayNumber):
        '''
        ' Compute the longitude of the sun at any time.
        ' Parameter: floating number jdn, the number of days since 1/1/4713 BC noon
        '''
        jdn = dayNumber - 0.5 - self.time_zone / 24
        time_in_julian = (jdn - 2451545.0) / 36525.
        # Time in Julian centuries
        # from 2000-01-01 12:00:00 GMT
        time_in_julian_2 = time_in_julian * time_in_julian
        degree_to_radian = np.pi / 180.  # degree to radian
        mean_time = (357.52910 + 35999.05030*time_in_julian
                     - 0.0001559*time_in_julian_2 -
                     0.00000048 * time_in_julian*time_in_julian_2)
        # mean anomaly, degree
        mean_degree = (280.46645 + 36000.76983*time_in_julian +
                       0.0003032*time_in_julian_2)
        # mean longitude, degree
        mean_long_degree = ((1.914600 - 0.004817*time_in_julian -
                             0.000014*time_in_julian_2)
                            * np.sin(degree_to_radian*mean_time))
        mean_long_degree += ((0.019993 - 0.000101*time_in_julian) *
                             np.sin(degree_to_radian*2*mean_time) +
                             0.000290*np.sin(degree_to_radian*3*mean_time))
        long_degree = mean_degree + mean_long_degree  # true longitude, degree
        long_degree = long_degree * degree_to_radian
        long_degree = long_degree - np.pi*2*(long_degree // (np.pi*2))
        # Normalize to (0, 2*np.pi)
        return (long_degree*6)// np.pi
    
    def lunar_month_nov(self,x):
        '''
        ' Find the 1st day of Lunar November for the given year and time zone
        '''
        gt = x.dt.strftime('%y-12-31')
        gt = pd.to_datetime(gt, yearfirst = True)
        off = self.julian_day(gt) - 2415021.
        k = off // 29.530588853
        lunar_month = self.new_moon(k)
        sun_long = self.sun_longitude(lunar_month)
        # sun longitude at local midnight
        lunar_month[sun_long>= 9] = self.new_moon(k - 1)
        return lunar_month
    
    def leap_offset(self, a11):
        '''
        ' Find the index of the leap month.
        '''
        k = np.floor((a11 - 2415021.076998695) / 29.530588853 + 0.5)
        last = 0
        arc = self.sun_longitude(self.new_moon(k))
        for j in range(len(a11)):
            i = 1  # start with month following lunar month 11
            arc[j] = self.sun_longitude(self.new_moon(k[j] + i))
            while True:
                last = arc[j]
                i += 1
                arc[j] = self.sun_longitude(self.new_moon(k[j] + i))
                if not (arc[j] != last and i < 14):
                    break
            arc[j] = i-1
        return arc
    
    def solar_to_lunar(self):
        date_table = pd.DataFrame(self.date_input)
        date_table.columns = ['solar_date']
        date_table['lunar_year'] = date_table['solar_date'].dt.strftime('%Y')
        date_table['lunar_year'] = date_table['lunar_year'].astype('int')
        day_number = self.julian_day(self.date_input)

        k = (day_number - 2415021.076998695) // 29.530588853

        month_start = self.new_moon(k + 1)
        month_start[month_start > day_number] = self.new_moon(k)

        a11 = self.lunar_month_nov(self.date_input)
        b11 = a11.copy()

        b11[a11 < month_start] = self.lunar_month_nov(self.date_input + pd.DateOffset(years=1))
        date_table['lunar_year'][a11 < month_start] = date_table['lunar_year'] + 1

        a11[a11 >= month_start] = self.lunar_month_nov(self.date_input- pd.DateOffset(years=1))

        lunar_day = day_number - month_start + 1
        date_table['diff'] = (month_start - a11) // 29
        date_table['lunar_month'] = date_table['diff'] + 11
        date_table['leap'] = False
        leap_month_diff = date_table['leap'].copy()  
        ss = b11-a11
        leap_month_diff[ss > 365] = self.leap_offset(a11)

        date_table['lunar_month'][(ss > 365)&(date_table['diff'] >= leap_month_diff)] = date_table['diff'] + 10
        date_table['leap'][(ss > 365)&(date_table['diff'] == leap_month_diff)] = True    
        date_table['lunar_month'][date_table['lunar_month'] > 12] = date_table['lunar_month'] - 12
        date_table['lunar_year'][(date_table['lunar_month'] >= 11)&(date_table['diff'] < 4)] = date_table['lunar_year']-1
        date_table['lunar_day'] = lunar_day
        date_table['lunar_date'] = pd.to_datetime(date_table[['lunar_year', 'lunar_month', 'lunar_day']].rename(columns={'lunar_year': 'year', 'lunar_month': 'month', 'lunar_day': 'day'}))
        
        return date_table[['solar_date','lunar_date','leap']]