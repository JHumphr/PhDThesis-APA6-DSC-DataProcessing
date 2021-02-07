%% Takes in the Results "Vector" from LoadSeparation.m to Compute Kamal-Sourour Model

%% Manual Exclusion of Data:
Exc = true(10,5,3);

% Temperature index:
% 130 - 1; 140 -2; 150 -3; 160-4; 170-5; 180 -6;   Index = (T-120)/10;
%
% Exclusion List is done manually, this is for clearly irrelvant data
% (negative, too small, etc.)
% Folder 1:
% T = 140: Exclude 2:
Exc(1,2,2) = false;
% T = 150: Exclude 3:
Exc(1,3,3) = false; %(this one is marginal)
% T = 160: Exclude 1:
Exc(1,4,1) = false;
% T = 170: Exclude 1:
Exc(1,5,3) = false;

% Folder 2:

%% 
load('Vector.mat')

for f = 1:1;%length(Vector) %for each folder
    V = Vector{f};
    Tspan = unique(Vector{f}(:,1));
    for i = 1:1:length(Tspan)
        T = Tspan(i);
        Tlogic = (V(:,1)==T);
        A=V(Tlogic,:); 
        Aavg = mean(A);
        Aavg(1) = T;
        Aavg(2) = NaN;
        Aavg
%           Exc(f,:,:)
        
        
   
    end
    
end