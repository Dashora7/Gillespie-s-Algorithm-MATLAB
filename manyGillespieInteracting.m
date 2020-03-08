kATS = 0.5; %Initialization of frequencies
kSTA = 5;
kSTW = 0.01;
kWTS = 0.01;
kMOVE = 1;
kENTER = 0.02;
p = 1;

time = 0; %Initialization of times
secCount = 0;
timeMax = 50000;
 
ATScount = 0; %Reaction counter initialization
STAcount = 0;
STWcount = 0;
WTScount = 0;
MOVEcount = 0;
ENTERcount = 0;
PUSHOUTcount = 0;
 
wState = [];
aState = [];
sState = [];
boxNumber = 10;
 
timeCount = 0;

wState = zeros(1, boxNumber);
aState = zeros(1, boxNumber + 1);
sState = zeros(1, boxNumber);
    
wState(1) = 0;
aState(1) = 1;
sState(1) = 0;

sState(4) = 1;

aState(7) = 1;

wState(5) = 2;

sState(10) = 1;
 
a0Array = [];
b0 = 0;
currentBox = 0;
 
savedTime = [];
aSaveLeft = [];
a = 0;

count = 0;
stateSum = 0;
 
while time <= timeMax
    a0Array = [];
    for c = 1:boxNumber
        stateSum = aState(c) + sState(c) + wState(c);
        kMOVE = 1/(1 + stateSum);
        a1 = kATS * aState(c); %Propensities for each reaction in vesicles/sec
        a2 = kSTA * sState(c);
        a3 = kSTW * sState(c);
        a4 = kWTS * wState(c);
        a5 = kMOVE * aState(c);
        a6 = 0;
        
        if (sState(c) == 1 & aState(c) == 1)
            a4 = 0;
        end
        
        if (c == 1 & aState(c) == 0 & sState(c) == 0)
            a6 = kENTER; 
        end
        
        a0 = a1 + a2 + a3 + a4 + a5 + a6; %Randomized calculation of time for each reaction
        a0Array = [a0Array, a0];
    end	
    
    
    b0 = sum(a0Array);
    r1 = rand;
    r2 = rand;   
    tau = -(1/b0)*(log(r1));
     
    for c=1:boxNumber
        if ((sum(a0Array(1:c))/b0 >= r2))
            currentBox = c;
            break;
        end
    end
        
        stateSum = aState(currentBox) + sState(currentBox) + wState(currentBox);
        kMOVE = 1/stateSum;

        a1 = kATS * aState(currentBox); %Propensities for each reaction in vesicles/sec
        a2 = kSTA * sState(currentBox);
        a3 = kSTW * sState(currentBox);
        a4 = kWTS * wState(currentBox);
        a5 = kMOVE * aState(currentBox);
        a6 = kENTER;
        
        if (sState(currentBox) == 1 & aState(currentBox) == 1)
            a4 = 0;
        end
        
        if (currentBox == 1 & aState(currentBox) == 0 & sState(currentBox) == 0)
            a0 = a1 + a2 + a3 + a4 + a5 + a6;
        else 
            a0 = a1 + a2 + a3 + a4 + a5;
        end
        sumA0 = sum(a0Array(1:(currentBox-1)));


        if((sumA0 + a1)/b0 >= r2)

            aState(currentBox) = aState(currentBox) - 1;
            sState(currentBox) = sState(currentBox) + 1;	
            ATScount = ATScount + 1;         

        elseif((sumA0 + a1 + a2)/b0 >= r2)  

            aState(currentBox) = aState(currentBox) + 1;
            sState(currentBox) = sState(currentBox) - 1;
            STAcount = STAcount + 1;

        elseif((sumA0 + a1 + a2 + a3)/b0 >= r2)

            sState(currentBox) = sState(currentBox) - 1;
            wState(currentBox) = wState(currentBox) + 1;
            STWcount = STWcount + 1;

        elseif ((sumA0 + a1 + a2 + a3 + a4)/b0 >= r2)

            wState(currentBox) = wState(currentBox) - 1;
            sState(currentBox) = sState(currentBox) + 1;
            WTScount = WTScount + 1;

        elseif ((sumA0 + a1 + a2 + a3 + a4 + a5)/b0 >= r2)
            if (not(currentBox == 10) & sState(currentBox + 1) == 1)
                if ((sumA0 + a1 + a2 + a3 + a4 + p) /b0 >= r2)
                    aState(currentBox) = aState(currentBox) - 1;
                    aState(currentBox + 1) = aState(currentBox + 1) + 1;
                    sState(currentBox + 1) = sState(currentBox + 1) - 1;
                    wState(currentBox + 1) = wState(currentBox + 1) + 1;
                else
                    aState(currentBox) = aState(currentBox) - 1;
                    sState(currentBox) = sState(currentBox) + 1;
                end
            
            elseif (not(currentBox == 10) & aState(currentBox + 1) == 1)
                aState(currentBox) = aState(currentBox) - 1;
                sState(currentBox) = sState(currentBox) + 1;
            else
            aState(currentBox) = aState(currentBox) -1;
            aState(currentBox + 1) = aState(currentBox + 1) + 1;
            MOVEcount = MOVEcount + 1;
            end
            
        else 
            aState(1) = aState(1) + 1;
            ENTERcount = ENTERcount +  1;
        end
            
    
    time = time + tau;
    secCount = secCount + tau;
    
    disp("SECOND: " + time);
    disp(aState);
    disp(sState);
    disp(wState);
    secCount = 0;  
    
    savedTime = [savedTime, time];
    %a = sum(aState)-sum(aState(1:boxNumber));
    aSaveLeft = [aSaveLeft, aState(boxNumber+1)];

end
%table(savedTime, aState(2), sState(2), wState(2));
figure;
plot(savedTime, aSaveLeft);
