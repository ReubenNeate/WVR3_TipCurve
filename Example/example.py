import TipCurves_v2 as tp #import the function

date, time, angle, t23, t31 = tp.ProcessWVR3("T240312.csv") #Run the ProcessWVR3 function 

index = 10 #The index of the data we want to preview. This will change the time of the data we preview

#print a preview of the data
print("date: ", date)
print("time: ", time[index])
print("        Brightness temperature")
print("angle   23.8 GHz                31.6 GHz ")
for a in range(8):
    print(angle[index][a], "° ", t23[index][a], ' K, ', t31[index][a], ' K' )
