% SHYAM'S TOLERANCE MODEL
clear all 

%% STEP 1) Define parameter ranges to look across
% these rules are necessary for the model to make biological sense:
%   all parameters must be positive
%   c <= 1
%   y <= 1
%   q <= 1
%   r >= 1
%   rp >= 1

% Below are some ranges that follow all the rules **I have no idea if they
% are broad enough (or too broad) for you to see the interesting behavior.
% You will need to adjust these ranges as you explore the model to see what
% the most informative ranges are.**
% You may need to use fewer values if this takes too long to run, or more
% values if you aren't getting enough resolution to see patterns.  This
% will take trial and error to figure out.
% For each parameter, you will specify
%   a minimum value : the step size : a maximum value
% For example, "rRange = 1:0.5:10" will have matlab run the model for r=1,
% r=1.5, r=2, ... r=10
rRange = 1:2:7;
    %% ***SEE NOTE 1***
yRange = 0.1:0.2:1;
pRange = 0.5:1:4.0;
cRange = 0.01:0.02:0.1;
qRange = 0.1:0.2:1;
rpRange = 1:0.5:3.5;
kRange = 0.5:1:4.0;
KRange = 50:25:150;a
bRange = 0.05:0.05:0.5;

%% STEP 2) Specify how long to run, initial conditions, etc.
time=200; % number of years the model is run

intS = 100; % initial above ground biomass
intV1=intS; 
intH =20; % initial herbivore density

S=zeros(time,1);
S1=zeros(time,1);
S2=zeros(time,1);
V1=zeros(time,1);
V2=zeros(time,1);
V3=zeros(time,1);
V=zeros(time,1);
S(1,1)=intS;
H=zeros(time,1);
H(1,1)=intH;
V1(1,1)=intV1;
%% ***SEE NOTE 2***


%% STEP 3) Make a file to store your results
resultsfile = fopen('ModelBehavior.txt','w'); %creates a text file to store results
% eventually, we'll put results in using the following columns:
% 1 - r value used
% 2 - y
% 3 - p
% 4 - c
% 5 - q
% 6 - rp
% 7 - k
% 8 - K
% 9 - b
% 10 - Final V value
% 11 - Final S value
% 12 - Final H value
% 13 - does S appear to be cycling? (0=no, 1=yes)
%       There are several different ways we could do this, but I have it
%       set up to check for cycling by seeing if the coefficient of
%       variation in S over the last 50 years is above some small threshold 
%       (10^-6).  If it is, we'll say that S is not holding steady at its equilibrium but
%       rather that it's showing sustained fluctuations.  You might have to
%       change both of these (the 10^-6 threshold and using the last 50
%       years) to capture the model's behavior properly -- you will have to
%       try and see
%       *** REMEMBER: if it does look like S is cycling, then you must
%       ignore what is in columns 10-12.  If the populations are stable,
%       these final values are the equlibrium values.  If the populations
%       are fluctuating, then they're just some point in the cycle where
%       you happened to stop the simulations, and they don't mean much.  If
%       you find a lot of cycling, you may want to change the code so that
%       you're saving the average densities (over the last 50 or so years)
%       instead of just the final densities.
% 14 - This will be all zeros for now, but it's a placeholder in case you
% want to use this last column to store some index of biocontrol
% effectiveness


%% STEP 4) Loop through all combinations of the parameters
for rIndex = 1:length(rRange)
    r = rRange(rIndex);
    for yIndex = 1:length(yRange)
        y = yRange(yIndex);
        for pIndex = 1:length(pRange)
            p = pRange(pIndex);
            for cIndex = 1:length(cRange)
                c = cRange(cIndex);
                for qIndex = 1:length(qRange)
                    q = qRange(qIndex);
                    for rpIndex = 1:length(rpRange)
                        rp = rpRange(rpIndex);
                        for kIndex = 1:length(kRange)
                            k = kRange(kIndex);
                            for KIndex = 1:length(KRange)
                                K = KRange(KIndex);
                                for bIndex = 1:length(bRange)
                                    b = bRange(bIndex);
[rIndex yIndex pIndex cIndex qIndex rpIndex kIndex KIndex bIndex]
% simulate the model with these parameter values
for t=2:time
        
       V1(t,1) = r*y*(S(t-1,1)/(1+b*S(t-1,1))); %early-spring growth
       S1(t,1)=S(t-1,1)-((y)*((S(t-1,1))/(1+b*S(t-1,1)))); %spring growth effect on Stored biomass
       H(t,1)=rp*H(t-1,1)*(V1(t-1,1))/(k*H(t-1,1)+V1(t-1,1)); %herbivore growth
       
       %% ***SEE NOTE 3***
       V2(t,1)=(1-H(t,1)/(H(t,1)+p*V1(t,1)))*V1(t,1); %spring growth minus loss of tissue due to herbivory
       V3(t,1)=V2(t,1)+(c*(H(t,1)/(H(t,1)+p*V1(t,1))*V1(t,1))*S1(t,1)) ; %late regrowth of above ground tissue after herbivory
       S2(t,1)=max([0,S1(t,1)-(c*(H(t,1)/(H(t,1)+p*V1(t,1))*V1(t,1))*S1(t,1))]);
       V(t,1)= r*(V3(t,1)/(1+(H(t,1)/(H(t,1)+p*V1(t,1))*V1(t,1))));
       S(t,1)=S2(t,1)+q*(V(t,1)*((1-(S2(t,1)/K))));
       %%
       
end

% calculate anything you want to know about how this simulation went
cvS = sqrt(var(S(time-50:time,1)))/mean(S(time-50:time,1)); % coefficient of variation in S over last 50 yrs
diffS = length(unique(sign(diff(S(time-50:time,1)))));
cycling = cvS>0.1&&cvS<1.0&&diffS>1; % will =1 if the CV of S is greater than the threshold of 10^-6; will eqaul 0 otherwise
Veffect = ((r*y)/b)/V(time,1); % change this to whatever quantity you want to use as an index of biocontrol effectiveness
Seffect = K/S(time,1);
Peffect = ((r*y/b) + K)/(V(time,1)+S(200,1));
% store the output before moving to the next parameter set
% the '%f\t ...' part just primes matlab to put numbers into columns
% after that, you give matlab the quantities you want to store, in the
% appropriate order.  It will add a new row to the saved file each time it
% gets to this line of code (so, each line is one distinct parameter set)
fprintf(resultsfile,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',r,y,p,c,q,rp,k,K,b,V(time,1),S(time,1),H(time,1),cycling,Veffect,Seffect,Peffect);


                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
fclose(resultsfile);
