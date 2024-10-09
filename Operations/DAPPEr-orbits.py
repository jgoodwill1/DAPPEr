#####
# created by Jarod Dagney on 2024 Oct. 09
# script to simulate amount of time a satellite will be above the horizon 
# for the University of Delaware CubeSat Ground Station, located at the
# Mount Cuba Astronomical Observatory
#
# use pip install beyond to install the library needed for making the simulation
# Things to change as needed:
#   Use a different ground station
#       station = create_station()
#           arg1: Name of ground station
#           arg2: list (latitude, longitude, altitude)
#   Change the satellite to be tracked
#       change the link in requests.get(...)
#           the number following satellite/tle/#####, replace the #s with the NORAD ID
#           of the satellite you want to track
#   Change the step size of the rough search (done to optimize run time)
#       change roughSearchStep to a number in seconds (see comments around the variable definition)
#####

from beyond.dates import Date, timedelta
from beyond.io.tle import Tle
from beyond.frames import create_station
import numpy as np
import requests
import matplotlib.pyplot as plt


station = create_station('Mt Cuba', (39.7845, -75.6327, 0)) 

#USE an api to get the TLE of the ISS: NORAD ID=25544
ret = requests.get(f"https://api.n2yo.com/rest/v1/satellite/tle/25544&apiKey=PA96AC-YFF725-V49R26-59PI").json()
tle = ret["tle"]
lines = tle.splitlines()

# Assign the lines to variables
line1 = lines[0]
line2 = lines[1]


tle = Tle(line1 + "\n" + line2)

start_time = Date.now()
simulation_end = timedelta(days = int(input("Enter the number of days you want to simulate in this script: ")))
minEl = int(input("What is the minimum elevation (in degrees) that you want to use?: "))

goodDates = []
roughSearchStep = 300 #number of seconds to do the rough search step. This should be GREATER than 60, otherwise you'll be here forever
#roughSearchStep should be less than the length of a pass (5-10 minutes), otherwise you may miss a pass entirely
# 300 is a good number: if a pass is less than 5 minutes (300 seconds) it's a poor pass. But that's assuming iss orbit

#do a rough search to find the times when the satellite is in view of the ground station.
for orb in station.visibility(tle.orbit(), start=start_time, stop=simulation_end, step=timedelta(seconds=roughSearchStep), events=True):
    goodDates.append(orb.date)


i = 0
eventsList = []
#loop over the times when the satellite is visible and group them into unique passes
#we just want one time for each pass to use to generate the full list of passes
while(i < len(goodDates) - 2):
    eventsList.append(goodDates[i])
    #if the difference is less than 1000 seconds this is 
    while(i < len(goodDates) - 1 and goodDates[i+1] - goodDates[i] < timedelta(seconds=roughSearchStep*4)):
        i += 1
        #print(i)
    i += 1

#eventsList now contains one entry for each pass: each entry is a date, use this date to generate the high resolutions passes

#start the pass generation search at least twice the amount of roughSearchStep before the time gathered from the large search until you reach LOS
passesList = []
for ev in eventsList:
    start_time = ev - timedelta(seconds=roughSearchStep * 2)
    #a pass won't last longer than 30 minutes so this will ensure we generate the full pass
    end_time = timedelta(minutes=30)
    gap = 3 #seconds: chosen as a pretty good resolution, can go higher for faster run speed (depends on your need for precision)
    myPass = []
    for orb in station.visibility(tle.orbit(), start=start_time, stop=end_time, step=timedelta(seconds = gap), events=True):
        myPass.append(orb)

        if orb.event and orb.event.info == "LOS":
            #print("found end of this pass")
            passesList.append(myPass)
            break
    
print(f"Found {len(passesList)} passes over {simulation_end.days} days")
#print(passesList)
#passesList now contains lists of the orbital event returned from Beyond, but only for the passes previously located

passDurations = []
for i in range(len(passesList)):
    azims = []
    elevs = []
    j = 0
    foundStart = False
    foundEnd = False
    #passStart = 0
    #passEnd = 0
    passStart = Date.now()
    passEnd = Date.now()
    while(j < len(passesList[i])):
        if(not foundStart and np.degrees(passesList[i][j].phi) > minEl):
            #print("Found the Start")
            foundStart = True
            passStart = passesList[i][j].date
        elif(not foundEnd and foundStart and np.degrees(passesList[i][j].phi) < minEl):
            #print("Found the end")
            foundEnd = True
            passEnd = passesList[i][j].date
        j+=1
    #print(f"{passEnd} - {passStart} = {passEnd - passStart}")
    #print(f'This pass has {(passEnd - passStart).seconds} seconds above {minEl} degrees')
    duration = (passEnd - passStart).seconds
    passDurations.append(duration)

print(passDurations)
#remove the 0s where there are no times above the minimum elevations
passDurations1 = []
for t in passDurations:
    if(t == 0):
        continue
    passDurations1.append(t)

print(len(passDurations1))
print(f'{np.average(passDurations1)} seconds per pass on average ')
print(f"Found {len(passDurations1)} good passes over {simulation_end.days} days")





    # for orb in passesList[i]:
    #     azims.append(-np.degrees(orb.theta) % 360)
    #     elevs.append(90 - np.degrees(orb.phi))




    # plt.figure()
    # ax = plt.subplot(111, projection='polar')
    # ax.set_theta_direction(-1)
    # ax.set_theta_zero_location('N')
    # plt.plot(np.radians(azims), elevs, '.')
    # ax.set_yticks(range(0, 90, 20))
    # ax.set_yticklabels(map(str, range(90, 0, -20)))
    # ax.set_rmax(90)
    # plt.show()


# for orb in station.visibility(tle.orbit(), start=start_time, stop=end_time, step=timedelta(seconds = gap), events=True):
#     azim = -np.degrees(orb.theta) % 360
#     elev = np.degrees(orb.phi)
#     r = orb.r / 1000.
#     print("{event:10} {tle.name}  {date:%Y-%m-%dT%H:%M:%S.%f} {azim:7.2f} {elev:7.2f} {r:10.2f}".format(
#         date=orb.date, r=r, azim=azim, elev=elev,
#         tle=tle, event=orb.event if orb.event is not None else ""
#     ))