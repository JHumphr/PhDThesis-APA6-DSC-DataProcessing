KSTable = KSTable;


for i = 1:10
   
    fileID = fopen('kamalParametersODE.txt','w');
    fprintf(fileID,'%2.1f & %2.1f & %4.3f & %4.3f & %4.3f & %4.3f & %3.2f & %3.2f & %5.4f & %5.4f\\\\ \r\n',KSTable'.*[1 1 1/1E5 1/1000 1/1E5 1/1000 1 1 1 1]'); % .*[1/1E5 1/1000 1/1E5 1/1000 1 1 1 1]'
    fclose(fileID);
    
end

for i = 1:length(KSTable(1,:))
KS_Average(i) = mean(KSTable(:,i));

end

KS_Average