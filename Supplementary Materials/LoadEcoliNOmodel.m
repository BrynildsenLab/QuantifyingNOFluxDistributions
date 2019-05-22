function EcoliNOmodel = LoadEcoliNOmodel( FileName )
% Reads model data from Excel file and converts to MATLAB structure.
%
%   If FILENAME is not specified, an open file dialog box will prompt user
%   for information.
%
%
% Copyright 2016 Jonathan L. Robinson


%% handle input arguments UNCOMMENT THIS SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!
% if nargin < 1
%     [FileName,FilePath] = uigetfile({'*.xlsx';'*.xls'});
% else
%     FilePath = [pwd,'\'];
% end
% FilePath = 'C:\Users\jlrtwo\Documents\Research\dynamic_RNS_network\Ecoli_NO_model\'; % REMOVE THIS REMOVE THIS REMOVE THIS

% FilePath = 'C:\Users\jlrtwo\Documents\Research\manuscript\clpP\model_files\'; % REMOVE THIS REMOVE THIS
% FileName = 'Ecoli_NO_model-ClpP-04.xlsx';  % REMOVE THIS REMOVE THIS REMOVE THIS

% FilePath = 'C:\Users\jlrtwo\Documents\Research\dynamic_RNS_network\Ecoli_NO_model\'; % REMOVE THIS REMOVE THIS
% FileName = 'Ecoli_NO_model.xlsx';  % REMOVE THIS REMOVE THIS REMOVE THIS

% FilePath = 'C:\Users\jlrtwo\Documents\Research\manuscript\LowO2detox\model_files\'; % REMOVE THIS REMOVE THIS
% FileName = 'Ecoli_NO_model-lowO2-01.xlsx';  % REMOVE THIS REMOVE THIS REMOVE THIS

% FilePath = 'C:\Users\jlrtwo\Documents\Research\manuscript\LowO2detox\model_files\'; % REMOVE THIS REMOVE THIS
% FileName = 'Ecoli_NO_model-lowO2-05_NEWTYPE.xlsx';  % REMOVE THIS REMOVE THIS REMOVE THIS

if nargin < 1
    [FileName,FilePath] = uigetfile({'*.xlsx';'*.xls'});
else
    FilePath = [pwd,'/'];
end

FullPath = [FilePath,FileName];


%% import model data from Excel file
[~,~,MODELDATA] = xlsread(FullPath);

% extract reaction names
[a,b] = find(strcmpi('Reaction Name',MODELDATA));  % locate reaction names
rxnNames = MODELDATA(a+1:end,b);  % extract reaction names
rxnNames(cellfun(@(x) any(isnan(x)),rxnNames)) = [];  % find and remove NaNs

% extract reaction equations
[a,b] = find(strcmpi('Reaction Equation',MODELDATA));  % locate reaction equations
rxnEqns = MODELDATA(a+1:end,b);  % extract reaction equations
rxnEqns(cellfun(@(x) any(isnan(x)),rxnEqns)) = [];  % find and remove NaNs

% extract rate equations
[a,b] = find(strcmpi('Rate Equation',MODELDATA));  % locate reaction rate equations
rateEqns = MODELDATA(a+1:end,b);  % extract reaction rate equations
rateEqns(cellfun(@(x) any(isnan(x)),rateEqns)) = [];  % find and remove NaNs

% extract reaction locations
[a,b] = find(strcmpi('Reaction Location',MODELDATA));  % locate reaction locations
rxnLocs = MODELDATA(a+1:end,b);  % extract reaction locations
rxnLocs(cellfun(@(x) any(isnan(x)),rxnLocs)) = [];  % find and remove NaNs

% extract parameter names
[a,b] = find(strcmpi('Parameter Name',MODELDATA));  % locate parameter names
paramNames = MODELDATA(a+1:end,b);  % extract parameter names
paramNames(cellfun(@(x) any(isnan(x)),paramNames)) = [];  % find and remove NaNs

% extract parameter values
[a,b] = find(strcmpi('Parameter Value',MODELDATA));  % locate parameter values
P = cell2mat(MODELDATA(a+1:end,b));  % extract parameter values and convert to matrix
P(isnan(P)) = [];  % find and remove NaNs

% extract species names
[a,b] = find(strcmpi('Species Name',MODELDATA));  % locate species names
species = MODELDATA(a+1:end,b);  % extract species names
species(cellfun(@(x) any(isnan(x)),species)) = [];  % find and remove NaNs

% extract species locations
[a,b] = find(strcmpi('Species Location',MODELDATA));  % locate species locations
specLocs = MODELDATA(a+1:end,b);  % extract species locations
specLocs(cellfun(@(x) any(isnan(x)),specLocs)) = [];  % find and remove NaNs

% extract initial species concentrations
[a,b] = find(strcmpi('Initial Concentration',MODELDATA));  % locate species concentrations
Xi = cell2mat(MODELDATA(a+1:end,b));  % extract concentrations and convert to matrix
Xi(isnan(Xi)) = [];  % find and remove NaNs

% extract cell volume fraction
[a,b] = find(strcmpi('cell volume fraction',MODELDATA));  % locate cell fraction
cellFrac = cell2mat(MODELDATA(a,b+1));  % extract cell fraction and convert to matrix

% extract media volume fraction
[a,b] = find(strcmpi('media volume fraction',MODELDATA));  % locate media fraction
mediaFrac = cell2mat(MODELDATA(a,b+1));  % extract media fraction and convert to matrix


%% construct stoichiometry matrix (rows = species, columns = reactions)

ignoreSpecies = {};  % some species concentrations do not change (e.g. h2o)
S = zeros(numel(species),numel(rxnEqns));  % initialize stoich matrix

for i = 1:numel(rxnEqns)

    rxnEqn = rxnEqns{i};  % process one reaction at a time
    
    spaceInd = [0,regexp(rxnEqn,'\s'),numel(rxnEqn)+1];  % obtain indices of spaces in equation
    arrowInd = regexp(rxnEqn,'>');  % obtain index of reaction arrow in equation
    Lparenth = regexp(rxnEqn,'(');  % obtain indices of left coefficient parentheses in equation
    Rparenth = regexp(rxnEqn,')');  % obtain indices of right coefficient parentheses in equation
    
    speciesCoeff = 1;  % set species stoichiometric coefficient to default value of 1
    
    for j = 1:numel(spaceInd)-1  % iterate through all spaces in equation
        
        % check if species is on reactant or product side of reaction equation
        if spaceInd(j) < arrowInd
            speciesSign = -1;
        else
            speciesSign = 1;
        end
        
        % check if current space is in front of coefficient bracket
        if ~isempty(Lparenth) && spaceInd(j) + 1 == Lparenth(1)
            
            % read coefficient value
            speciesCoeffText = rxnEqn(Lparenth(1)+1:Rparenth(1)-1);
            speciesCoeff = str2double(speciesCoeffText);
            
            % remove coefficient index from list of coefficient bracket indices
            Lparenth(1) = [];
            Rparenth(1) = [];
            continue  % skip to next iteration of for-loop, saving species coefficient value
        
        end
        
        % check if string between current space and next space (rxnPiece) is a species name
        rxnPiece = rxnEqn(spaceInd(j)+1:spaceInd(j+1)-1);
        [~,speciesInd] = ismember(rxnPiece,species);
        
        if speciesInd > 0
            % add species coefficient to stoichiometry matrix
            S(speciesInd,i) = speciesSign * speciesCoeff;
        elseif ~any(ismember('<>=+',rxnPiece)) && ~isempty(rxnPiece)
            % document ignored/static species (ignore symbols <,>,=,+)
            ignoreSpecies = [ignoreSpecies;rxnPiece];
        end
        
        speciesCoeff = 1;  % reset species coefficient value to default value of 1
        
    end
end

S(end+1,:) = 0;  % add dummy species for indexing purposes in v-vector construction


%% construct compartment-scaling matrix for species (Fspec)
% this matrix is used to convert species concentrations from a
% per-total-volume basis to a basis according to their location 
% (e.g., a per-cell-volume basis)

Fspec = zeros(size(S,1),1);  % initialize Fspec vector

% verify that all of the species locations are valid
if ~all(ismember(specLocs,{'cell','media','all'}))
    error('Species location vector contains invalid entry.');
end

if cellFrac > 0  % skip if cell fraction is zero, to avoid division by zero.
    Fspec(ismember(specLocs,'cell')) = 1/cellFrac;
end
Fspec(ismember(specLocs,'media')) = 1/mediaFrac;
Fspec(ismember(specLocs,'all')) = 1;  % 'all' are species present in all compartments

Fspec = diag(Fspec);  % assign to diagonal of matrix


%% construct compartment-scaling matrix for species (Fspec)
% this matrix is used to convert reaction fluxes to a per-total-volume basis

Frxn = zeros(size(rxnLocs));  % initialize Frxn vector

% verify that all of the reaction locations are valid
if ~all(ismember(rxnLocs,{'cell','media','all'}))
    error('Reaction location vector contains invalid entry.');
end

Frxn(ismember(rxnLocs,'cell')) = cellFrac;
Frxn(ismember(rxnLocs,'media')) = mediaFrac;
Frxn(ismember(rxnLocs,'all')) = 1;  % reactions that are in all compartments

Frxn = diag(Frxn);  % assign to diagonal of matrix


%% construct elementary rate equation matrix and complex rate equation vector
vElement = zeros(numel(rxnEqns),1);  % initialize
vComplex = {};  % initialize 

for i = 1:numel(rateEqns)  % iterate through each rate equation
    currentEqn = rateEqns{i};
    
    % determine locations of brackets (enclose species or parameter names)
    Lbracket = regexp(currentEqn,'[');  
    Rbracket = regexp(currentEqn,']');
    
    for j = 1:numel(Lbracket)  % iterate through each set of brackets
        
        term = currentEqn(Lbracket(j)+1:Rbracket(j)-1);
        [~,speciesInd] = ismember(term,species);
        [~,paramInd] = ismember(term,paramNames);
        
        if speciesInd == 0 && paramInd == 0 && ~isempty(term) % return error if species or parameter is missing
            
            fprintf('\nWARNING! Reaction rate depends on %s, but species or parameter is not in list provided.\n\n',term);
            continue  % skip to next iteration of for-loop
            
        elseif ismember('{',currentEqn)  % rate equation is complex if equation is enclosed in curly braces
            
            % construct complex rate expression
            if speciesInd ~= 0
                newVar = ['X(',num2str(speciesInd),')'];
            elseif paramInd ~= 0
                newVar = ['P(',num2str(paramInd),')'];
            end
            currentEqn = [currentEqn(1:Lbracket(j)-1),newVar,currentEqn(Rbracket(j)+1:end)];
            indOffset = (Rbracket(j) - Lbracket(j) + 1) - length(newVar);
            Lbracket = Lbracket - indOffset;
            Rbracket = Rbracket - indOffset;
        
        else  % rate equation is elementary
            
            vElement(i,j) = speciesInd;
            
        end
        
    end
    
    if ismember('{',currentEqn)
        % remove curly brackets from equation
        currentEqn(1) = []; 
        currentEqn(end) = [];
        vComplex = [vComplex;currentEqn];
    end
    
end

% The bottom entries of the elementary rate vector (vElement) are not used,
% because those rates are accounted for in the complex rate vector (vComplex).
if ~isempty(vComplex)
    vElement((end-numel(vComplex)+1):end,:) = [];
end

Xi(end+1) = 1;  % make dummy species with concentration = 1
vElement(vElement == 0) = numel(species)+1;  % this is to handle indices of zero
species{end+1} = 'DUMMY_SPECIES';

% generate complex rate equation string from vComplex vector
vComplexEquations = 'v2 = [';
for i = 1:length(vComplex)
    vComplexEquations = [vComplexEquations,vComplex{i},';'];
    if i == length(vComplex)
        vComplexEquations(end) = [];
        vComplexEquations = [vComplexEquations,'];'];
    end
end

%% package model contents

EcoliNOmodel.species = species;
EcoliNOmodel.paramNames = paramNames;
EcoliNOmodel.rxnEqns = rxnEqns;
EcoliNOmodel.rxnNames = rxnNames;
EcoliNOmodel.rxnLocs = rxnLocs;
EcoliNOmodel.S = S;
EcoliNOmodel.Fspec = Fspec;
EcoliNOmodel.Frxn = Frxn;
EcoliNOmodel.Xi = Xi;
EcoliNOmodel.vComplex = vComplex;
EcoliNOmodel.vComplexEquations = vComplexEquations;
EcoliNOmodel.vElement = vElement;
EcoliNOmodel.P = P;

fprintf('\nLoaded model data from "%s".\n\n',FileName);
