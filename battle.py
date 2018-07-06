
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from numpy.random import randn
import statistics
import math
import time

arr=[]

N = input('Enter the matrix size(20-70k): ')
rho_d = input('Enter the value of rho-a(0.2-0.8): ')
rho_a = input('Enter the value of rho-b(0.2-0.8): ')
strategy = raw_input('Enter A for aggressive strategy and D for defensive strategy: ')


for ii in range(1,50,1):
    #declaring matrix size and initial conditions
    #N = 50    #matrix size
    R = ii      #value of range
    a0 = 50000   #initial reserve for attackers
    d0 = 50000   #initial reserve for defenders
    A0 = a0/2000  #parameter for colour settings
    plt.pause(1)
    data = np.zeros((N+(2*R),N+(2*R)))      #2d matrix for Battleground
    #rho_a = 0.2             #initial deployment probability for attackers
    #rho_d = 0.45            #deployment probability for defenders
    b1 = 1                  #attacking/defensive power of value 1
    b2 = 2                  #attacking/defensive power of value 2
    itr = 200               #number of iteretions
    flag_eradicated = 0
    #code for simulation of cellular automata
    fig, ax = plt.subplots()
    cax = ax.imshow(data[R:N+R,R:N+R], interpolation='nearest', cmap=cm.bwr, vmin=-A0, vmax=A0)
    ax.set_title('Battlefield')
    cbar = fig.colorbar(cax, ticks=[-100, 0, 100])
    plt.pause(1)

    #Initialize Defenders randomly with rho_d density in lattice
    a_reserve = a0
    d_reserve = d0
    for i in range(R,N+R):
        for j in range(R,N+R):
            r = np.random.rand()
            if r<rho_d:
                r = np.random.rand()
                if r<0.5:
                    data[i][j] = data[i][j] + b2
                    d_reserve = d_reserve-b2
                else:
                    data[i][j] = data[i][j] + b1
                    d_reserve = d_reserve-b1

    #Initialize Attackers randomly with rho_a density in lattice
    for i in range(R,N+R):
        for j in range(R,N+R):
            r = np.random.rand()
            if r<rho_a:
                r = np.random.rand()
                if r<0.5:
                    data[i][j] = data[i][j] - b2
                    a_reserve = a_reserve-b2
                else:
                    data[i][j] = data[i][j] - b1
                    a_reserve = a_reserve-b1

    plt.pause(1)
    cax = ax.imshow(data[R:N+R,R:N+R], interpolation='nearest', cmap=cm.coolwarm, vmin=-A0, vmax=A0)


    flagA = 0
    flagB = 0
    pA = 0          #Total no. of attackers in ground
    pD = 0          #Total no. of defenders in ground

    for it in range(itr):
        #Deployment of attackers
        for i in range(R,N+R):
            for j in range(R,N+R):
                #sum of defenders on 4 quadrants surrounding the attack site
                sum1 = 0;
                sum2 = 0;
                sum3 = 0;
                sum4 = 0;
                r = np.random.rand()
                #If attacker site then only deploy
                if data[i][j]<0:
                    #calculating sum in all four qudrants
                    for k in range(i-R,i+1):
                        for l in range(j-R,j+1):
                            sum2 = sum2 + data[k][l];
                    for k in range(i-R,i+1):
                        for l in range(j,j+R+1):
                            sum1 = sum1 + data[k][l];
                    for k in range(i,i+R+1):
                        for l in range(j-R,j+1):
                            sum3 = sum3 + data[k][l];
                    for k in range(i,i+R+1):
                        for l in range(j,j+R+1):
                            sum4 = sum4 + data[k][l];

                    sums = [sum1, sum2, sum3, sum4]
                    if strategy=='A':
                        max_sum = max(sums)
                    elif strategy=='D':
                        max_sum = min(sums)
                    b=b2
                    #Deploy attackers
                    if sum1 == max_sum:
                        if ((sum3 == max_sum) or (sum2==max_sum and sum4==max_sum)):
                            data[i][j] = data[i][j]
                        elif (sum2 == max_sum):
                            data[i-1][j] = data[i-1][j] - b
                        elif (sum4 == max_sum):
                            data[i][j+1] = data[i][j+1] - b
                        else:
                            data[i-1][j+1] = data[i-1][j+1] - b
                    elif sum2 == max_sum:
                        if (sum4 == max_sum) or (sum1==max_sum and sum3==max_sum):
                            data[i][j] = data[i][j]
                        elif (sum1 == max_sum):
                            data[i-1][j] = data[i-1][j] - b
                        elif (sum3 == max_sum):
                            data[i][j-1] = data[i][j-1] - b
                        else:
                            data[i-1][j-1] = data[i-1][j-1] - b
                    elif sum3 == max_sum:
                        if (sum1 == max_sum) or (sum2==max_sum and sum4==max_sum):
                            data[i][j] = data[i][j]
                        elif (sum2 == max_sum):
                            data[i][j-1] = data[i][j-1] - b
                        elif (sum4 == max_sum):
                            data[i+1][j] = data[i+1][j] - b
                        else:
                            data[i+1][j-1] = data[i+1][j-1] - b
                    elif sum4 == max_sum:
                        if (sum2 == max_sum) or (sum1==max_sum and sum3==max_sum):
                            data[i][j] = data[i][j]
                        elif (sum1 == max_sum):
                            data[i][j+1] = data[i][j+1] - b
                        elif (sum3 == max_sum):
                            data[i+1][j] = data[i+1][j] - b
                        else:
                            data[i+1][j+1] = data[i+1][j+1] - b

                    a_reserve = a_reserve-b
                    if a_reserve<0:
                        flagA=1
                if flagA==1:
                    break
            if flagA==1:
                break



        #deploy defenders
        for i in range(R,N+R):
            for j in range(R,N+R):
                r = np.random.rand()
                if r<rho_d:
                    r = np.random.rand()
                    if r<0.5:
                        data[i][j] = data[i][j] + b2
                        d_reserve = d_reserve-b2
                        if d_reserve<0:
                            flagB=1
                    else:
                        data[i][j] = data[i][j] + b1
                        d_reserve = d_reserve-b1
                        if d_reserve<0:
                            flagB=1
                if flagB==1:
                    break
            if flagB==1:
                break
        if flagB==1 or flagA==1:
            break

        cax = ax.imshow(data[R:N+R,R:N+R], interpolation='nearest', cmap=cm.coolwarm, vmin=-A0, vmax=A0)
        plt.pause(0.05)

        #calculation total number of sites occupied by defenders and attackers
        siteA = 0
        siteD = 0
        for i in range(R,N+R):
            for j in range(R,N+R):
                if data[i][j]<0:
                    siteA = siteA+1
                elif data[i][j]>0:
                    siteD = siteD+1
        if siteD==0:
            print("defenders are eradicated - ATTACKERS WON!")
            flag_eradicated = 1
            break
        if siteA==0:
            print("attackers are eradicated - DEFENDERS WON!")
            flag_eradicated = 1
            break
        #PA.append(a0)
        # PD.append(pD)
        # SITEA.append(siteA)
        # SITED.append(siteD)


    siteA = 0
    siteD = 0
    total_reserve = 0
    for i in range(R,N+R):
        for j in range(R,N+R):
            total_reserve = total_reserve + data[i][j]
            if data[i][j]<0:
                pA = pA-data[i][j]
                siteA = siteA + 1
            elif data[i][j]>0:
                pD = pD+data[i][j]
                siteD = siteD + 1
    pA = pA + a_reserve
    pD = pD + d_reserve
    if flag_eradicated==0:
        if a_reserve>0 and d_reserve<=0:
            print("Defenders ran out of resources : Attackers will win")
        elif a_reserve<=0 and d_reserve>0:
            print("Attackers ran out of resources : Defenders will win")
        elif a_reserve<=0 and d_reserve<=0:
            print("Both ran out of resourses : Match Draw")
    # print siteA, siteD, pA, pD
    total_reserve = total_reserve + a_reserve + d_reserve
    x=float(siteA*100)/(N*N)    #percentage sites occupied by attackers
    print(ii,pA,pD,a_reserve, d_reserve,siteA,siteD, total_reserve,x)
    arr.append(x)
plt.figure(2)
plt.plot(arr)
plt.ylabel('Percentage sites occupied by attackers')
plt.xlabel('Iteration')
plt.show()
