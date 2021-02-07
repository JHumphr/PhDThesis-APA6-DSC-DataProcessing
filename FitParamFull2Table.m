% FitParamFull to table
load FitParamFull.mat

gaussTable = []; variC1Table = []; variC20Table = [];
tempAKS = [0 0 0 0 0 0 0 0 0 0];
c1 = [0.6 0.8 1.0 1.2 1.4 1.2 1.2 1.2 1.2 1.2 0.0];
c20p = [0.6 0.6 0.6 0.6 0.6 0.4 0.6 0.8 1.0 1.2 0.0];

test= [c1; c20p]; %test concetenate


FitParamFull{11}.kParam =zeros(8,8);
for i = 1:11 % for each folder + 1 for table spacing
% Gauss    
    if i == 11
    else
        tempAGauss = FitParamFull{i}.gParam; % create temporary array for gauss parameters 5x11 (8 gaus std Error Hmelt Htotal)
        gaussTable = [gaussTable; zeros(1,length(tempAGauss)); tempAGauss]; 
    end
% Kamal-Sourour
    kParam = FitParamFull{i}.kParam(:,:);

    for j = 1:length(kParam(:,1)) % for each set of fits
%         disp(kParam(j,:))
        tempAKS(i+(j-1).*(11),:) = [c1(i) c20p(i) kParam(j,:)];         
    end
% Kim-Avrami
    if i == 11
    else
    kAParam = FitParamFull{i}.kAParam;
    if i <= 5
        variC1 = squeeze(kAParam(3,:,:));
        rowNum = size(variC1,1);
        while rowNum < 5
        variC1 = [variC1; zeros(1,size(variC1,2))];
        rowNum = size(variC1,1);
        end
        
        variC1Table = [variC1Table; variC1]; 
    end
    if i >= 6 && i < 11
        variC20 = squeeze(kAParam(3,:,:));
        rowNum = size(variC20,1);
        while rowNum < 5
        variC20 = [variC20; zeros(1,size(variC20,2))];
        rowNum = size(variC20,1);
        end
        
        variC20Table = [variC20Table; variC20];
    end
    end
end

fileID = fopen('gaussParameters.txt','w');
fprintf(fileID,'%4.3E %4.3E %4.3E %4.3E %4.3E %4.3E %4.3E %4.3E %4.3E %4.3E %4.3E\r\n',gaussTable');
fclose(fileID);

% disp(tempAKS)
fileID = fopen('kamalParameters.txt','w');
fprintf(fileID,'%2.1f & %2.1f & %4.3f & %4.3f & %4.3f & %4.3f & %3.2f & %3.2f & %5.4f & %5.4f\\\\ \r\n',tempAKS'.*[1 1 1/1E5 1/1000 1/1E5 1/1000 1 1 1 1]'); % .*[1/1E5 1/1000 1/1E5 1/1000 1 1 1 1]'
fclose(fileID);

% disp(variC1Table)
fileID = fopen('kimParametersC1.txt','w');
fprintf(fileID,'%4.3f %3.2f %4.5f %3.5f %4.5f %4.5f \r\n',variC1Table');
fclose(fileID);

% disp(variC20Table)
fileID = fopen('kimParametersC20.txt','w');
fprintf(fileID,'%4.3f %3.2f %4.5f %3.5f %4.5f %4.5f \r\n',variC20Table');
fclose(fileID);
