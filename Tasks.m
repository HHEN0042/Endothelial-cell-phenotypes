
initCobraToolbox (0)
changeCobraSolver('gurobi') 

model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 
Model1_I = readCbModel('Model1_I.mat');
Model2_I = readCbModel('Model2_I.mat');
Model3_I = readCbModel('Model3_I.mat');
Model4_I = readCbModel('Model4_I.mat');


TasksTable=readtable('Supp_Tables_2_ENERGY METABOLISM_redox.xls','sheet','Supp_Table_1_Recon1_Filtered');
Tasks=unique(TasksTable(:,1));
% Table, Headers, and colums:
headers=split(unique(join(TasksTable{:,2:4},'_'),'stable'),'_');
data = cell(height(Tasks)+2,length(model.rxns)+3);
for i=1:length(headers)
    data(i+2,1)=headers(i,1);
    data(i+2,2)=headers(i,2);
    data(i+2,3)=headers(i,3);
end
for i=1:length(model.rxns)
    data(1,i+3)=model.subSystems(i);
    data(2,i+3)=model.rxns(i);
end
FunctionSummaryTrauma1=transpose(data);
FunctionSummaryTrauma2=transpose(data);
FunctionSummaryTrauma3=transpose(data);
FunctionSummaryTrauma4=transpose(data);
for i=1:height(Tasks)
    %Ith Function
    IthIndex=find(ismember(TasksTable(:,1),Tasks(i,1)));
    IthTask=join(string(TasksTable{IthIndex(1),2:4}),newline);
%   IthTaskShort=TasksTable{IthIndex(1),4};
    IthTaskShort=join(string(TasksTable{IthIndex(1),2:4}),'_');    
    IthSubstrate=table2cell(TasksTable(IthIndex,5));
    IthProduct=table2cell(TasksTable(IthIndex,8));
    %Testing trauma groups
    [Flux_1, FBAsolution_1,~]=testPathway(Model1_I,IthSubstrate,IthProduct);
    [Flux_2, FBAsolution_2,~]=testPathway(Model2_I,IthSubstrate,IthProduct);
    [Flux_3, FBAsolution_3,~]=testPathway(Model3_I,IthSubstrate,IthProduct);
    [Flux_4, FBAsolution_4,~]=testPathway(Model4_I,IthSubstrate,IthProduct);
    for j=1:length(Model1_I.rxns)
        FunctionSummaryTrauma1{j+3,i+2}=FBAsolution_1.x(j);
        FunctionSummaryTrauma2{j+3,i+2}=FBAsolution_2.x(j);
        FunctionSummaryTrauma3{j+3,i+2}=FBAsolution_3.x(j);
        FunctionSummaryTrauma4{j+3,i+2}=FBAsolution_4.x(j);
    end
    subplot(ceil(sqrt(height(Tasks))),floor(sqrt(height(Tasks))),i)
    hold on;    
    t=bar([Flux_1,Flux_2,Flux_3,Flux_4]);
    ylabel('Flux Size');
    xlabel('Trauma groups');
    title(IthTask);
    ax=gca;
    %exportgraphics(ax,strcat('Path_An-',string(IthTaskShort),'.pdf')); 
end
