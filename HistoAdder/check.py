data_crab_info=[]

with open('Data_UL2018_crabjob_pass_fractions.txt') as f:
    for line in f:
        x = line.split()
        data_crab_info.append(dict({"name":x[0], "pass":float(x[1]))
