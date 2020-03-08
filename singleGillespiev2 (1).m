kATS = 0.5; %Initialization of frequencies 
kSTA = 5; 
kSTW = 0.01;
kWTS = 0.01;

time = 0; %Initialization of times
timeMax = 10;

Wnum = 0; %Number of vesicles in each state
Anum = 1;
Snum = 0;

ATScount = 0; %Reaction counter initialization
STAcount = 0;
STWcount = 0;
WTScount = 0;

Wtau = 0; %Total times in each state
Atau = 0;
Stau = 0;

savedState = [];
savedTime = [];
savedDist = [];

position = 0;

while time <= timeMax
    
    a1 = kATS * Anum; %Propensities for each reaction in vesicles/sec
    a2 = kSTA * Snum;
    a3 = kSTW * Snum;
    a4 = kWTS * Wnum;
    
    a0 = a1 + a2 + a3 + a4; %Randomized calculation of time for each reaction
    r1 = rand;
    tau = -(1/a0)*(log(r1));
    position = position + tau*(Anum == 1);
    %disp(tau);
    r2 = rand;
    
    Atau = (tau * Anum) + Atau; %Change of time in each states
    Wtau = (tau * Wnum) + Wtau;
    Stau = (tau * Snum) + Stau; 
    
    if((a1/a0) >= r2)
        
        Anum = Anum - 1;
        Snum = Snum + 1;    
        ATScount = ATScount + 1;
        savedState = [savedState, "S"];
        
    elseif((a1 + a2)/a0 >= r2)      
        Anum = Anum + 1;
        Snum = Snum - 1;
        STAcount = STAcount + 1;
        savedState = [savedState, "A"];
        
    elseif((a1 + a2 + a3)/a0 >= r2)
        
        Snum = Snum - 1;
        Wnum = Wnum + 1;
        STWcount = STWcount + 1;
        savedState = [savedState, "W"];
        
    else
        
        Wnum = Wnum - 1;
        Snum = Snum + 1;
        WTScount = WTScount + 1;
        savedState = [savedState, "S"];
    end
    
    time = time + tau;
    
    savedTime = [savedTime, time];
    savedDist = [savedDist, position];
    %disp(time);  
end

array = [Wtau; Atau; Stau];
disp(array);
figure; 
bar(array);
figure;
plot(savedTime, savedDist);

