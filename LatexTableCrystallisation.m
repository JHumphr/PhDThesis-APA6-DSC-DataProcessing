load FitParamFull

c1 = [0.6 0.8 1.0 1.2 1.4 1.2 1.2 1.2 1.2 1.2 0.0];
c20p = [0.6 0.6 0.6 0.6 0.6 0.4 0.6 0.8 1.0 1.2 0.0];

fileID = fopen('KimAvramiTable.txt','w');

    
%     fprintf(fileID,'n = %d & & & & &\\\\ \\hline\r\n',i);
%     fprintf(fileID,'%s & %s & %s & %s & %s & %s\\\\ \r\n','T','aeq','K','$\theta$','nc','StdErr');
    fprintf(fileID,'%s & %s & %s & %s & %s \\\\ \r\n','c','$\psi$','$n_c$','$\theta$','StdErr');
    fprintf(fileID,'%s & %s & %s & %s & %s \\\\ \\hline\r\n','[s]','[-]','[-]','[s]','[-]');
    
for j = 1:10 %each parametet set
    fprintf(fileID,'C1 = %6.2f; C20P = %6.2f & & & & \\\\\r\n',c1(j),c20p(j));
    for i = 1:1:4 %1:1:4 %for each value of nc    
        
        kA = FitParamFull{j,1}.kAParam;  % kA for each folder
        
        kAatNc = squeeze(kA(i,:,:)); % kA of the folder for the given nc
        
        fprintf(fileID,'%6.2d & %6.0f & %6.0f & %6.2d & %6.3f\\\\\r\n',kAatNc');
    end
  
end

fclose(fileID);