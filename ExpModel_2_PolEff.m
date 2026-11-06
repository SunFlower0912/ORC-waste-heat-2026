function [out,TS] = ExpModel_2_PolEff(fluid, P_su, h_su, M_dot, P_ex, T_amb, param)

% fluid, P_su, h_su, N_exp, P_ex, T_amb, param % <-- Unconmment to see the inputs of the external call

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 11/05/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% "ExpModel_2_PolEff.m" is a matlab code implementing a polynomial model
% of volumetric expanders (see the Documentation/ExpModel_2_SemiEmp_MatlabDoc)
% Unlike "ExpModel_PolEff.mat" which takes as input N_exp and computes M_dot, 
% this code derives the expander rotational speed based on the inlet mass flow rate.
%
% The model inputs are:
%       - P_su: inlet pressure of the WF                          	[Pa]
%       - h_su: inlet temperature of the WF                        	[J/kg]
%       - P_ex: outlet pressure of the WF                          	[Pa]
%       - fluid: nature of the fluid (string)                       [-]
%       - M_dot: fluid mass flow rate                               [kg/s]
%       - T_amb : ambien temperature                                [K]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' should contain the followings:
%           param.V_s , machine displacement volume                	[m3]
%           param.V, volume of the expander                         [m^3]
%           param.coeffPol_is, polynmial coef for epsilon_is        [-]
%           param.coeffPol_ff, polynmial coef for FF      [-]
%           param.displayResults, flag to display the results or not[1/0]
%           param.h_min, minimum enthalpy of the fluid              [J/kg]
%           param.h_max, maximum enthalpy of the fluid              [J/kg]
%
% The model outputs are:
%       - out: a structure variable which includes
%               - T_ex =  exhaust temperature                    [K]
%               - h_ex =  exhaust enthalpy                       [J/kg]
%               - M_dot = fluid mass flow rate                   [kg/s]
%               - W_dot = mechanical power                       [W]
%               - Q_dot_amb = ambiant losses                     [W]
%               - epsilon_is = isentropic efficiency             [-]
%               - FF = filling factor (volumetric efficiency)    [-]
%               - time = the code computational time             [sec]
%               - flag = simulation flag                         [-1/1]
%               - M = mass of refrigerent inside the expaner     [kg]
%
%       - TS : a stucture variable which contains the vectors of temperature
%              and entropy of the fluid (useful to generate a Ts diagram 
%              when modelling the entire ORC system 
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% DEMONSTRATION CASE -- COMMENT THIS SECTION IF EXTERNAL CALL FOR SPEED IMPROVEMENT
 if nargin == 0
    % Define a demonstration case if ExpanderModel.mat is not executed externally
     fluid = 'R245fa';                           %Nature of the fluid%     
     M_dot = 0.05511;                            %Mass flow rate         [kg/s]
     P_su = 12.77e+05;               %Supply pressure        [Pa]
     h_su = refpropm('H','T',273+112,'P',P_su*0.001,fluid);               %Supply enthalpy        [J/kg]
     P_ex = 4.21e+05;               %Exhaust pressure       [Pa]
     T_amb = 298.1500;                           %Ambient temperature    [K]
     param.coeffPol_is  = [0.640418889478694;0.121516496926235;-0.002040204176348;-0.023473811185544;-2.795162060638594e-04;-2.386060475108433e-05;-8.375764103129725e-05;-9.589286061751416e-09;2.463118782749123e-05;1.299052138854458e-06];
     param.coeffPol_ff = [1.753986206634610;0.093097577652062;-0.005812558544030;-0.007288772042005;-3.214683370757326e-04;4.726357455792485e-05;-2.830506820656102e-04;1.764406891583263e-08;-1.525193171879006e-06;1.759235588693687e-07];
     param.AU_amb = 0.2580;
     param.V_s = 1.279e-05;
     param.V = 1.279e-05*4.05;
     param.h_max =  refpropm('H','T',500,'P',4e3,fluid);
     param.h_min =  refpropm('H','T',253.15,'P',5e1,fluid);
     N_exp = 4420
 end

%% MODELLING SECTION
tstart_exp = tic;

T_su = refpropm('T','P',P_su * 0.001,'H',h_su,fluid);
s_su = refpropm('S','P',P_su * 0.001,'H',h_su,fluid);
rho_su = refpropm('D','P',P_su * 0.001,'H',h_su,fluid);
h_ex_s = refpropm('H','P',P_ex * 0.001,'S',s_su,fluid);

if P_su > P_ex && h_su > refpropm('H','P',P_su * 0.001,'Q',0,fluid)
    
    if length(param.coeffPol_ff) == 6
        % if expander with constant speed -> quadratic function of Rp and rho_su
        FF = max(0.5, min(param.coeffPol_ff(1) + param.coeffPol_ff(2)*(P_su/P_ex) + param.coeffPol_ff(3)*rho_su + param.coeffPol_ff(4)*(P_su/P_ex)^2 + param.coeffPol_ff(5)*(P_su/P_ex)*rho_su + param.coeffPol_ff(6)*(rho_su)^2, 5));
    elseif length(param.coeffPol_ff) == 10
        % if expander with variable speed -> quadratic function of Rp, rho_su and N_exp
        FF = max(0.3,min(param.coeffPol_ff(1) + param.coeffPol_ff(2)*(P_su/P_ex) + param.coeffPol_ff(3)*rho_su + param.coeffPol_ff(4)*(P_su/P_ex)^2 + param.coeffPol_ff(5)*(P_su/P_ex)*rho_su + param.coeffPol_ff(6)*(rho_su)^2 + param.coeffPol_ff(7)*N_exp + param.coeffPol_ff(8)*N_exp^2 + param.coeffPol_ff(9)*N_exp*(P_su/P_ex) + param.coeffPol_ff(10)*N_exp*rho_su,5));
    end
     
    if length(param.coeffPol_is) == 6
        % if expander with constant speed -> quadratic function of Rp and rho_su
        epsilon_is = 0.6*max(0.01,min(param.coeffPol_is(1) + param.coeffPol_is(2)*(P_su/P_ex) + param.coeffPol_is(3)*rho_su + param.coeffPol_is(4)*(P_su/P_ex)^2 + param.coeffPol_is(5)*(P_su/P_ex)*rho_su + param.coeffPol_is(6)*(rho_su)^2,1));
    elseif length(param.coeffPol_is) == 10
        % if expander with variable speed -> quadratic function of Rp, rho_su and N_exp
        epsilon_is = max(0.01,min(param.coeffPol_is(1) + param.coeffPol_is(2)*(P_su/P_ex) + param.coeffPol_is(3)*rho_su + param.coeffPol_is(4)*(P_su/P_ex)^2 + param.coeffPol_is(5)*(P_su/P_ex)*rho_su + param.coeffPol_is(6)*(rho_su)^2 + param.coeffPol_is(7)*N_exp + param.coeffPol_is(8)*N_exp^2 + param.coeffPol_is(9)*N_exp*(P_su/P_ex) + param.coeffPol_is(10)*N_exp*rho_su,1));
    end
   N_exp = 60*M_dot/(param.V_s*FF*rho_su);
    W_dot = M_dot*(h_su - h_ex_s)*epsilon_is;
    Q_dot_amb = max(0,param.AU_amb*(T_su - T_amb));
    h_ex = h_su - (W_dot + Q_dot_amb)/M_dot;
    if h_ex > param.h_min && h_ex < param.h_max
        out.flag = 1;
    else
        out.flag = -1;
    end  
       
else
    out.flag = -2;
end

if out.flag > 0
    out.h_ex = h_ex;
    out.N_exp = N_exp;
    out.W_dot = W_dot;
    out.Q_dot_amb = Q_dot_amb;
    out.epsilon_is = epsilon_is;
    out.FF = FF;
    out.T_ex = refpropm('T','P',P_ex * 0.001,'H',out.h_ex,fluid);
    out.M = (refpropm('D','H',h_su,'P',P_su * 0.001,fluid)+refpropm('D','H',out.h_ex,'P',P_ex * 0.001,fluid))/2*param.V;   
else
    out.N_exp = 60*M_dot/(param.V_s*rho_su);
    out.FF = 1;
    out.W_dot = M_dot*(h_su-h_ex_s);
    out.epsilon_is =1;
    out.Q_dot_amb = 0;
    out.h_ex = h_ex_s;
    out.T_ex = refpropm('T','P',P_ex * 0.001,'H',out.h_ex,fluid);
    out.M = (refpropm('D','H',h_su,'P',P_su * 0.001,fluid)+refpropm('D','H',out.h_ex,'P',P_ex * 0.001,fluid))/2*param.V;
end

TS.T = [T_su out.T_ex];
TS.s = [s_su refpropm('S','H', out.h_ex,'P', P_ex * 0.001, fluid)];
out.time = toc(tstart_exp);

end