function [fitVals,resNorm,exitFlag,resid,J] = ...
    FitModelParams(model,expX,expY,fitParams,paramBounds,fitSpecies,variance,normFlag)
% Vary parameter values to fit simulation output to experimental data.
%
% *** NOTE: Initial guess will be the current parameter value in model ***
%
% ---------------- INPUTS ----------------
%
% model         structure containing all model information (obtained using
%               the LoadEcoliNOmodel function).
%
% expX          X-values from experimental data.
%
% expY          Y-values from experimental data.
%
% fitParams     list of metabolite concentrations or parameter values to
%               vary in order to fit the experimental data.
%
% paramBounds   N x 2 matrix of lower and upper bounds for each varied
%               parameter.
%
% fitSpecies    Species concentration data to which the model is being fit
%
% variance      Variance of experimental Y-values. Data points with greater
%               variance will be weighted less in the optimization.
%
% normFlag      if TRUE, normalizes each dataset by the mean value of the
%               experimental Y-values for that dataset.
%
% Copyright 2016 Jonathan Robinson

% handle undefined inputs
if nargin < 8
    normFlag = false;
    if nargin < 7
        variance = [];
    end
end

% convert fitParams to cell array if it is only a string
if ~iscell(fitParams)
    fitParams = {fitParams};
end

% define model species and parameters (should be the same for all models)
specs = model(1).species;
params = model(1).paramNames;

% combine metabolite concentration and parameter value vectors for easier indexing
XP_list = [specs;params];
for i = 1:numel(model)
    % parameter values may vary between different models
    XP(:,i) = [ model(i).Xi; model(i).P ];
end


for i = 1:length(fitParams)
    % the index of the varied parameters will be the same across all models
    if strcmp(fitParams{i},'dna')
        % this is to change all DNA base concentrations simultaneously
        varParamInd{i,1} = find(strcmp('ds_dad_2',XP_list(:,1)) | strcmp('ds_dgsn',XP_list(:,1)) | strcmp('ds_dcyt',XP_list(:,1)));
    else
        varParamInd{i,1} = find(strcmp(fitParams{i},XP_list(:,1)));
    end
end

% normalize Y-values if option is selected
if ( normFlag )
    for i = 1:size(expY,2)
        normVal(i) = mean(expY(:,i)); 
        expY(:,i) = expY(:,i)./normVal(i);
        variance(:,i)=variance(:,i)/((normVal(i))^2);
    end
else
    normVal = [];
end

% normalize Y-values by variance if provided
if ~isempty( variance )
    for i = 1:size(expY,2)
        % check if one or multiple Xdata columns are being used
        if size(expX,2) > 1
            [~,final_ind] = max(expX(:,i));  % find last point, in case expX vectors differ in length
        else
            final_ind = length(expX);
        end
        expY(1:final_ind,i) = expY(1:final_ind,i) ./ sqrt(variance(1:final_ind,i));
    end
end


% add data to model (primarily to facilitate data transfer into nested function)
for i = 1:numel(model)
    model(i).paramLB = paramBounds(:,1);
    model(i).paramUB = paramBounds(:,2);
    model(i).varParamInd = varParamInd;
    model(i).XP = XP(:,i);
    model(i).XPlist = XP_list;
    model(i).normFlag = normFlag;
    model(i).normVal = normVal;
    model(i).fitSpecies = fitSpecies;
    model(i).variance = variance;
end

x0 = zeros(length(varParamInd),1);
% varied parameters should have the same value across all model types, so
% the initial values will be the same for all models
for i = 1:length(varParamInd)
    % normalize parameter values to a range of 0 to 1
    x0(i,1) = ( XP(varParamInd{i}(1),1) - model(1).paramLB(i) ) / ( model(1).paramUB(i) - model(1).paramLB(i) );
end

if any( x0 < 0 | x0 > 1 )
    error('Initial parameter value does not fall within given bounds.');
end

% options = ['Diagnostics','on'];
options = optimset('DiffMinChange',0.01,'TolX',0.0001,'TolFun',0.001);
% options = [];

LB = zeros(size(x0));
UB = ones(size(x0));
[fitVals,resNorm,resid,exitFlag,~,~,J] = lsqcurvefit(@(x,xdata) bacon(x,xdata,model),x0,expX,expY,LB,UB,options);

% ONLY RETURN ONE COLUMN OF RESIDUALS (if simultaneous opt)
resid = resid(:,1);

% [fitVals,resid,J,covb,mse] = nlinfit(expT,expNO,@(x,xdata) bacon(x,xdata,paramData),x0,options);
% exitFlag = [];
% resNorm = [];


% convert parameter values from normalized values to real values
fitVals = fitVals.*(model(1).paramUB - model(1).paramLB) + model(1).paramLB;



function Y = bacon(x,expX,model)

%############# DISPLAY ALGORITHM PROGRESS (DISABLE ON CLUSTER) ############
% for i = 0:ceil(length(x)/10)-1;
%     if i < ceil(length(x)/10)-1
%         fprintf('%1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f\n', ...
%             x(i*10+1:(i+1)*10));
%     else
%         fprintf('%1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f  %1.3f', ...
%             x(i*10+1:end));
%     end
% end
% fprintf('\n');
%##########################################################################

XP = zeros(numel(model(1).XP),numel(model));  % initialize XP matrix
for i = 1:numel(model)
    XP(:,i) = model(i).XP;
end

pvals = x.*(model(1).paramUB - model(1).paramLB) + model(1).paramLB;  % convert to param value from fraction
varParamInd = model(1).varParamInd;

for i = 1:length(varParamInd)
    XP(varParamInd{i},:) = pvals(i);  % update parameter values
end

% assign species concentrations and parameter values to model(s)
for i = 1:numel(model)
    model(i).Xi = XP(1:length(model(i).Xi),i);
    model(i).P = XP(length(model(i).Xi)+1:end,i);
end

% run simulation(s) and extract species concentration and time data
if strcmp(model(1).fitSpecies,'hmp')  % special case if total Hmp concentration is desired
    [~,y_ind] = ismember({'hmpFe2';'hmpFe3';'hmpFe2h';'hmpFe3h'; ...
        'hmpFe2h2';'hmpFe3h2';'hmpFe2_o2';'hmpFe2_no';'hmpFe3_onoo'; ...
        'hmpFe2h_o2';'hmpFe2h_no';'hmpFe3h_onoo';'hmpFe2h2_o2'; ...
        'hmpFe2h2_no';'hmpFe3h2_onoo'},model(1).species);
else
    [~,y_ind] = ismember(model(1).fitSpecies,model(1).species);
end

for i = 1:numel(model)
    
    % check if one or multiple Xdata columns are being used
    if size(expX,2) > 1
        expXind = i;
        [tfinal,final_ind] = max(expX(:,expXind));
        maxStep = mean(diff(expX(1:final_ind,expXind)));
    else
        expXind = 1;
        tfinal = max(expX);
        maxStep = mean(diff(expX));
    end
    
    % run simulation
    [t_temp,X_temp] = RunEcoliNOmodel(model(i),tfinal,maxStep);
    
    % process simulation output
    t{i} = t_temp;
    X{i} = X_temp;
    if strcmp(model(1).fitSpecies,'hmp')
        y{i} = sum(X{i}(:,y_ind),2);
    else
        y{i} = X{i}(:,y_ind);
    end
    
end


% select data points from simulation to compare with experimental data
% NOTE EXAMPLE LAYOUT:
%                        model 1      model 1      model 2      model 2
%                       species A    species B    species A    species B 
%              expY = [     #      ,     #      ,     #      ,     #    
%                           #      ,     #      ,     #      ,     #     
%                           #      ,     #      ,     #      ,     #     ];
Y = zeros(size(expX,1),numel(model)*numel(y_ind));
for i = 1:numel(model)
    for j = 1:numel(y_ind)
        if size(expX,2) > 1
            expXind = numel(y_ind)*(i-1) + j;
            [~,final_ind] = max(expX(:,expXind));
        else
            expXind = 1;
            final_ind = length(expX);
        end
        % interpolate to obtain simulation values at experimental x-values
        Y(1:final_ind,expXind) = interp1(t{i},y{i}(:,j),expX(1:final_ind,expXind));
    end
end


% normalize Y-values if option is selected
if ( model(1).normFlag )
    for i = 1:size(Y,2)
        Y(:,i) = Y(:,i)./model(1).normVal(i);
    end
end


% normalize Y-values by variance if provided
variance = model(1).variance;
if ~isempty( variance )
    for i = 1:size(Y,2)
        % check if one or multiple Xdata columns are being used
        if size(expX,2) > 1
            [~,final_ind] = max(expX(:,i));  % find last point, in case expX vectors differ in length
        else
            final_ind = length(expX);
        end
        Y(1:final_ind,i) = Y(1:final_ind,i) ./ sqrt(variance(1:final_ind,i));
    end
end


