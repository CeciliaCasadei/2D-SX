# -*- coding: utf-8 -*-

def get_shannonSamplings(sampling, qMax):
    samplings = []
    for l in range(0, 1000):
        q_rod_sampling = l * sampling
        if q_rod_sampling <= qMax:
            samplings.append(q_rod_sampling)
        else:
            samplings.append(q_rod_sampling)  # ADD ONE MORE
            break        
    n = len(samplings)     # ONLY NON_NEGATIVE
    for i in range(1, n):
        samplings.append(-samplings[i])            
    samplings.sort()  
    return samplings    