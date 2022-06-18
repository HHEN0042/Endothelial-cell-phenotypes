initCobraToolbox (0)
changeCobraSolver('gurobi') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Characterize the Baseline model   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: %Random sampling of hpmecgen model a level 3 sized model with level 4a constraints !!!!!!I would say that this and the above three lines are all the 1st step random sampling
%model = readCbModel('iEC2997_preclinicalmodel_annsurg');%%%% The Ccat-model should be used for the clinical models 
load('iEC3006')
model = ExpandedModelNewBoundaries
FBA = optimizeCbModel(model);%%%%This FBA is part of the random sampling protocol it is used to establish the set the boundary of the possible flux predictions to consider
model.lb(find(model.c)) = 0.5*FBA.f;%%%%Here we find half the values of the maximal predicted fluxes
[sampleMetaOutC, mixedFraction] = gpSampler(model, length(model.rxns), [], 8*3600,10000)%,4);%%%%Here we take the half maximal model and set up the random sampling
save mixedFraction_EndoRecon1NewBounds mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_EndoRecon1NewBounds FBA
pointsC = sampleMetaOutC.points;
save sampleMetaOut_EndoRecon1NewBounds.mat sampleMetaOutC
% STEP 2: List of each metabolits after random sampling 
%%% Loaded the model after the random sampling 
model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 
%%%Set up a place to put the next bit
FluxesPreclinical = zeros(length(model.rxns),7);
%%%% Get the stats about the reaction fluxes the 25th and 75th percentiles
%%%% (columns 6 and 7) would be enough but the rest may be interesting
for idx=1:length(model.rxns)    
    flux = model_points(idx,:);        
    minf = min(flux);
    meanf = mean(flux);
    maxf = max(flux);
    stdf = std(flux);
    UPstd = meanf + (stdf*2);
    LOWstd = meanf - (stdf*2);
    percent25 = prctile(flux,25);
    percent75 = prctile(flux,75);    
    FluxesPreclinical(idx,1) = minf;
    FluxesPreclinical(idx,2) = meanf;
    FluxesPreclinical(idx,3) = maxf;
    FluxesPreclinical(idx,4) = LOWstd;
    FluxesPreclinical(idx,5) = UPstd;
    FluxesPreclinical(idx,6) = percent25;
    FluxesPreclinical(idx,7) = percent75;        
end
clear('ans' , 'idx' , 'minf' , 'meanf' , 'mean(flux)' , 'maxf' , 'stdf' , 'UPstd' , 'LOWstd' , 'percent25' , 'percent75');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Define Boundaries in Trauma models    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL 1
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC model1.xlsx', 1, 'A1:B53'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
modelTrauma1_I=modelTrauma;
testFBA_model1_I = optimizeCbModel(modelTrauma1_I); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save Trauma-1_I modelTrauma1_I

%% MODEL 2
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC model2.xlsx', 1, 'A1:B53'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end 
% STEP 5: Use FBA to see if the model with the trauma constraints works
modelTrauma2_I=modelTrauma;
testFBA_model2_I = optimizeCbModel(modelTrauma2_I); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save Trauma-2_I modelTrauma2_I

%% MODEL 3
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC model3.xlsx', 1, 'A1:B53'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
modelTrauma3_I=modelTrauma;
testFBA_model3_I = optimizeCbModel(modelTrauma3_I); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save Trauma-3_I modelTrauma3_I

%% MODEL 4
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC model4.xlsx', 1, 'A1:B53'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
modelTrauma4_I=modelTrauma;
testFBA_model4_I = optimizeCbModel(modelTrauma4_I); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save Trauma-4_I modelTrauma4_I


%% MODEL 1
% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(modelTrauma1_I);
model.lb(find(modelTrauma1_I.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(modelTrauma1_I, length(modelTrauma1_I.rxns), [], 8*3600,10000);
save mixedFraction_modelTrauma1_I mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_modelTrauma1_I fba
pointsC = sampleMetaOutC.points;
save points_modelTrauma1_I pointsC
save Model1_I.mat sampleMetaOutC

%% MODEL 2
% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(modelTrauma2_I);
model.lb(find(modelTrauma2_I.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(modelTrauma2_I, length(modelTrauma2_I.rxns), [], 8*3600,10000);
save mixedFraction_modelTrauma2_I mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_modelTrauma2_I fba
pointsC = sampleMetaOutC.points;
save points_modelTrauma2_I pointsC
save Model2_I.mat sampleMetaOutC

%% MODEL 3
% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(modelTrauma3_I);
model.lb(find(modelTrauma3_I.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(modelTrauma3_I, length(modelTrauma3_I.rxns), [], 8*3600,10000);
save mixedFraction_modelTrauma3_I mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_modelTrauma3_I fba
pointsC = sampleMetaOutC.points;
save points_modelTrauma3_I pointsC
save Model3_I.mat sampleMetaOutC

%% MODEL 4
% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(modelTrauma4_I);
model.lb(find(modelTrauma4_I.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(modelTrauma4_I, length(modelTrauma4_I.rxns), [], 8*3600,10000);
save mixedFraction_modelTrauma4_I mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_modelTrauma4_I fba
pointsC = sampleMetaOutC.points;
save points_modelTrauma4_I pointsC
save Model4_I.mat sampleMetaOutC

