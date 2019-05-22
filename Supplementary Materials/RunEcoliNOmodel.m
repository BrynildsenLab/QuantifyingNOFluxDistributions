function  [t,X,dX,V] = RunEcoliNOmodel(EcoliNOmodel,tfinal,maxStep)
% Run simulation with EcoliNOmodel.
%
%--------------------------------- INPUTS ---------------------------------
%
% EcoliNOmodel      Ecoli model structure generated using
%                   LoadEcoliNOmodel.m function.
%
% tfinal            Final simulation time, in hours [DEFAULT = 1 hour]. 
%                   Initial time is automatically set as zero.
%
% maxStep           Maximum time step size [hr] allowed during integration.
%                   [DEFAULT = No limit]
%
%--------------------------------- OUTPUTS --------------------------------
%
% t     Vector of simulation time points [h]
%
% X     Matrix of species concentrations at each time point [uM]
%           NOTE: species concentrations are relative to the volume of the
%           compartment in which they exist
%
% dX    Matrix of changes in species concentrations per time at each time 
%       point [uM/h]
%           NOTE: species concentrations are relative to the volume of the
%           compartment in which they exist
%
% V     Matrix of reaction rates calculated at each time point [uM/h]
%           NOTE: reaction rates are relative to the volume of the
%           compartment in which they occur
%
%
% Copyright 2014 Jonathan L. Robinson

if nargin < 3
    maxStep = [];
    if nargin < 2
        tfinal = 1;
    end
end

% unpack EcoliNOmodel
[rxnNames,S,Xi,vElement,P,Fspec,Frxn] = ...
    deal(EcoliNOmodel.rxnNames,  EcoliNOmodel.S, EcoliNOmodel.Xi, ...
         EcoliNOmodel.vElement, EcoliNOmodel.P, EcoliNOmodel.Fspec, ...
         EcoliNOmodel.Frxn);

tspan = [0,tfinal];

% scale stoichiometry matrix
S = Fspec*S*Frxn;


% set ODE parameters. Negative concentrations are not permitted.
if isempty(maxStep)
    options = odeset('NonNegative',1:numel(Xi),'AbsTol',1e-8);
else
    options = odeset('NonNegative',1:numel(Xi),'AbsTol',1e-8,'MaxStep',maxStep);
end
   
% integrate governing mass balances over specified time span.
[t,X] = ode15s(@ODEfunc,tspan,Xi,options,P,S,vElement);

% calculate concentration changes at each time point, if requested
V = zeros(length(t),length(rxnNames));
dX = zeros(size(X));
if ( nargout > 2 )
    for i = 1:length(t)
        [dx,v] = ODEfunc(t(i),X(i,:),P,S,vElement);
        V(i,:) = v;
        dX(i,:) = dx;
    end
end

% remove dummy species from output
X(:,end) = [];
dX(:,end) = [];

% ALPHA ESTIMATION
X(:,5)=X(:,5)*P(11);



function [dX,v] = ODEfunc(~,X,P,S,vElement)

% calculate dX/dt for all elementary-type rate equations
Pk = P(1:size(vElement,1));
v1 = Pk.*prod(X(vElement),2);


%*************** COMPLEX RATE EQUATIONS (vComplexEquations) ***************
%
% NOTE: see model user guide for additional details.
%
% > If the model Excel document is changed (e.g., reaction added, species
%   removed, etc.), re-load the model with "LoadEcoliNOmodel", and paste
%   the updated complex rate equations (vComplexEquations) below.

%ALPHA ESTIMATION
v2 = [((P(1)*X(5)*X(5)*X(10)));P(2)*X(5)*X(7);P(3)*X(3);P(4)*X(3);P(5)*X(7)*X(7);P(6)*X(4);P(7)*X(4);P(8)*X(9);P(9)*X(5);P(10)*(X(11)-X(10))];

%ECOLI MODEL
%v2 = [X(129)*X(19)*X(130)*P(161)/(X(19)*X(130) + P(162)*X(130) + P(163)*X(19));X(129)*X(19)*X(132)*P(161)/(X(19)*X(132) + P(162)*X(132) + P(163)*X(19));X(132)*X(121)*P(164)/(P(165) + X(121));X(134)*X(121)*P(164)/(P(165) + X(121));X(46)*X(39)*P(166)/(P(167) + X(39));X(46)*X(37)*P(168)/(P(169)+ X(37));X(47)*X(38)*P(170)/(P(171) + X(38));X(48)*X(30)*P(172)/(P(173) + X(30));X(48)*X(28)*P(172)/(P(173) + X(28));X(48)*X(29)*P(172)/(P(173) + X(29));P(174)*X(49)*X(24)*X(33)/(P(176)*P(178) + P(175)*X(24) + P(178)*X(33)*(1+P(176)/P(177)) + X(33)*X(24)*(1+P(175)/P(177)));P(174)*X(49)*X(22)*X(31)/(P(176)*P(179) + P(175)*X(22) + P(179)*X(31)*(1+P(176)/P(177)) + X(31)*X(22)*(1+P(175)/P(177)));P(174)*X(49)*X(23)*X(32)/(P(176)*P(180) + P(175)*X(23) + P(180)*X(32)*(1+P(176)/P(177)) + X(32)*X(23)*(1+P(175)/P(177)));X(50)*X(36)*X(92)*P(181)/(X(36)*X(92) + P(183)*X(36) + P(182)*X(92));X(50)*X(34)*X(92)*P(181)/(X(34)*X(92) + P(183)*X(34) + P(182)*X(92));X(50)*X(35)*X(92)*P(181)/(X(35)*X(92) + P(183)*X(35) + P(182)*X(92));P(184)*X(64)*X(55)/(X(55) + P(185));X(63)*(P(186)*X(61)*X(95) + P(189)*X(61)^2*X(95)) / ( P(187)*X(61) + P(188)*X(95) + X(61)*X(95) + P(190)*X(61)^2 + P(191)*X(61)^2*X(95) );P(192)*X(141)*X(139)*X(95)/(X(139)*X(95) + P(193)*X(139) + P(194)*X(95));P(195)*X(140)*X(55)/(P(196) + X(55));P(197)*(X(113)-X(112));P(198)*X(17)*X(148)^2*X(112)/(P(199)*X(112) + P(200)*X(148)*X(112) + P(201)*X(148)^2 + X(148)^2*X(112));P(202)*X(15)*X(148)^2*X(112)/(P(203)*X(112) + P(204)*X(148)*X(112) + P(205)*X(148)^2 + X(148)^2*X(112));P(206)*X(100)*X(17)/(1 + X(112)/P(201)) - P(207)*X(18);P(208)*X(100)*X(15)/(1 + X(112)/P(205)) - P(209)*X(16);((X(117) + X(118))/(X(117) + X(118) + X(119) + X(120) + X(121) + X(122)))*P(210)*X(149)*X(147)*X(93)/(P(211)*P(212) + X(147)*X(93) + P(212)*X(147) +P(213)*X(93));P(214)*X(150)*X(147)*X(93)/(P(215)*P(216) + X(147)*X(93) + P(216)*X(147) +P(217)*X(93));P(218)*X(80)*X(93)/(1+X(100)/P(219));P(220)*X(81)*X(100)/(P(221) + X(100));P(222)*X(82)*X(100)/(X(100)+P(223));P(224) + X(100)/(X(100) + P(226))*(P(225)-P(224));P(227) + X(100)/(X(100) + P(229))*(P(228)-P(227));P(230)*P(232)*X(102)/(P(231)*P(232)+P(232)*X(102)+P(231)*X(112)+X(102)*X(112));P(233)*X(151)*(1 + P(236)*X(112)/(P(237)+X(112)));P(234)*X(152)*(1 + P(236)*X(112)/(P(237)+X(112)));P(235)*X(153)*(1 + P(236)*X(112)/(P(237)+X(112)))];
%**************************************************************************


v = [v1;v2];  % construct v-vector
v(isnan(v)) = 0;  % eliminate erroneous outputs


dX = S*v;



