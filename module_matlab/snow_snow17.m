function [OUT_EP,OUT_SWE,OUT_MELT,OUT_ADJPR,OUT_INIT] = snow_snow17(IN_SDATE, IN_EDATE, IN_PRCP, IN_TAVG, IN_ELEV, IN_PAR, IN_INIT)
%#codegen
%   Snow process model based on SNOW17 (Anderson, 1976)
%
%   INPUTS
%       IN_SDATE    Start date of simulation
%       IN_EDATE    End date of simulation
%       IN_PRCP     Time series of daily precipitation
%       IN_TAVG     Time series of daily average temperature in Celsius
%       IN_ELEV     Average elevation of the hydrologic response unit
%       IN_PAR      Snow 17 model parameters
%       IN_INIT     Snow storage initial states
%
%   OUTPUTS
%       OUT_EP      Effective precipitation: the amount of precipitation added in the soil
%	    OUT_SWE     State of snow storage: snow water equivalent
%	    OUT_MELT    Amount of snow melt by degree day factor
%	    OUT_ADJPR   Adjusted precipitation used by the snow module
%	    OUT_INIT    Snow storage states at the end
%
%   Model Parameter Description
%	    IN_PAR(1)   SCF     Multiplying factor which adjusts Precipitation (accounts for gage snow catch deficiencies)
%	    IN_PAR(2)   PXTEMP	Temperature that separates rain from snow, deg C
%	    IN_PAR(3)   TTI     Temperature interval for mixture of snow and rain
%	    IN_PAR(4)   MFMAX   Maximum melt factor during non-rain periods - assumed to occur on June 21
%	    IN_PAR(5)   MFMIN   Minimum melt factor during non-rain periods - assumed to occur on Dec 21
%	    IN_PAR(6)   UADJ    Average wind function during rain-on-snow periods
%	    IN_PAR(7)   MBASE   Base temperature for snowmelt computations during non-rain periods, deg C
%	    IN_PAR(8)   TIPM    Antecedent temperature index parameter (0.01 to 1.0)
%	    IN_PAR(9)   PLWHC   Percent liquid water holding capacity (maximum value allowed is 0.4)
%	    IN_PAR(10)  NMF     Maximum negative melt factor
%	    IN_PAR(11)  DAYGM   A constant daily rate of melt at the soil-snow interface
%
%   Initial Snow Storage Condition
%	    IN_INIT(1)  W_i     Accumulated water equivalent of the ice portion of the snow cover (mm)
%       IN_INIT(2)  ATI     Antecedent Temperature Index, deg C
%       IN_INIT(3)  W_q     Liquid water held by the snow (mm)
%       IN_INIT(4)  Deficit Heat Deficit, also known as NEGHS, Negative Heat Storage
% 
%   The directive %#codegen indicates that this function is intended for C/C++ code generation
%   codegen syntax to generate MEX function:
%	    cfg = coder.config('mex');
%       cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%       codegen snow_snow17 -args {zeros(1,3),zeros(1,3),coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),0,zeros(1,11),zeros(1,4)} -config cfg -report
%
% 
%--------------------------------------------------------------------------
% 
%   References
%	    Anderson, E. A. (1976), A point energy and mass balance model of a snow cover, NOAA Tech. Rep. NWS 19, 150 pp., Natl. Oceanic and Atmos.Admin., Silver Spring, Md.
%       Anderson, E. A. (2006), Snow accumulation and ablation model - SNOW17, http://www.nws.noaa.gov/oh/hrl/nwsrfs/users_manual/part2/_pdf/22snow17.pdf
%   
%--------------------------------------------------------------------------
% 
%	Author: Sungwook Wi
%	e-mail: sw2275@cornell.edu
%__________________________________________________________________________ 


%---------------------------
%	Initiailizing Outputs
%---------------------------
L           =  length(IN_PRCP);
OUT_EP      =  nan(L,1);  
OUT_SWE     =  nan(L,1); 
OUT_MELT    =  nan(L,1);  
OUT_ADJPR	=  nan(L,1);    
OUT_INIT    =  nan(4,1);


%------------------------ 
%	Loading Parameters
%------------------------ 
SCF    = IN_PAR(1);
PXTEMP = IN_PAR(2);
TTI    = IN_PAR(3);
MFMAX  = IN_PAR(4);
MFMIN  = IN_PAR(5);
UADJ   = IN_PAR(6);  
MBASE  = IN_PAR(7);
TIPM   = IN_PAR(8); 
PLWHC  = IN_PAR(9); 
NMF    = IN_PAR(10); 
DAYGM  = IN_PAR(11); 


%--------------------------- 
%	Loading Initial State
%--------------------------- 
W_i     = IN_INIT(1); 
ATI     = IN_INIT(2); 
W_q     = IN_INIT(3); 
Deficit = IN_INIT(4);


%--------------------------- 
%	Julian Date Array
%---------------------------
julian  = nan(L,2); % [year juliandate]
yr      = (IN_SDATE(1):IN_EDATE(1));
n       = 1;
for iy = 1:length(yr) 
    
    % Julian date generation
    if iy == 1 
 
        if ( IN_SDATE(2) <= 2 ) % January & February
            year  = IN_SDATE(1) - 1.0;
            month = IN_SDATE(2) + 12.0;
        else
            year = IN_SDATE(1);
            month = IN_SDATE(2);
        end
        day = IN_SDATE(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = IN_SDATE(1) - 1.0;
        month = 1 + 12.0;
        day = 1;
        jd_2 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        str_jdate = jd_1 - jd_2 + 1;
    else
        str_jdate = 1;
    end
    
    if iy == length(yr)
        
        if ( IN_EDATE(2) <= 2 ) % January & February
            year  = IN_EDATE(1) - 1.0;
            month = IN_EDATE(2) + 12.0;
        else
            year = IN_EDATE(1);
            month = IN_EDATE(2);
        end
        day = IN_EDATE(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year = IN_EDATE(1) - 1.0;
        month = 1 + 12.0;
        day = 1;
        jd_2 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        end_jdate = jd_1 - jd_2 + 1;
    else
        if ((mod(yr(iy),4) == 0 && mod(yr(iy),100) ~= 0) || mod(yr(iy),400) == 0)
            end_jdate = 366;
        else
            end_jdate = 365;
        end
    end
    
    jdate       = (str_jdate:end_jdate);
    numday      = length(jdate);
    
    julian(n:n+numday-1,1) = jdate';
    julian(n:n+numday-1,2) = repmat(yr(iy),end_jdate-str_jdate+1,1);
    
    n = n + numday;   
end


%------------------------ 
%	Execute Simulation
%------------------------ 
dtt = 24; % time interval of temperature data
dtp = 24; % time interval of prcipitation data
for t = 1:L
    
    Ta = IN_TAVG(t);     % Air temperature at this time step (deg C)
    Pr = IN_PRCP(t);     % Precipitation at this time step (mm)
    

    %                        FORM OF PRECIPITATION 
    % --------------------------------------------------------------------- 
    if Ta >= (PXTEMP+TTI)   % All rain, no snow
        SNOW = 0; 
        RAIN = Pr;
    elseif Ta <= PXTEMP % All snow, no rain
        SNOW = Pr;  
        RAIN = 0;       
    else  % Linear mixture of snow and rain in interval TTI
        snowfrac = -1/TTI * (Ta-PXTEMP) + 1;
        SNOW = Pr*snowfrac;  
        RAIN = Pr*(1-snowfrac);
    end


    %                    ACCUMULATION OF THE SNOW COVER 
    % ---------------------------------------------------------------------           
    Pn  = SNOW*SCF;   % Water equivalent of new snowfall (mm)
    W_i = W_i + Pn;   % Water equivalent of the ice portion of the snow cover (mm)
    

    %      ENERGY EXCHANGE AT SNOW/AIR SURFACE DURING NON-MELT PERIODS 
    % ---------------------------------------------------------------------

    % Seasonal variation in the non-rain melt factor
    DAYN = julian(t,2);      % Current julian date
    if ((mod(julian(t,1),4) == 0 && mod(julian(t,1),100) ~= 0) || mod(julian(t,1),400) == 0) % Leap year
        days=366;
        N_Mar21=DAYN-81;   % Day of year since March 21 (leap)
    else   % Not a leap year
        days=365;
        N_Mar21=DAYN-80;
    end
    Sv = (0.5*sin((N_Mar21 * 2 * pi)/days)) + 0.5;        % Seasonal variation
    Av = 1.0;                                             % Seasonal variation adjustment, Av=1.0 when lat < 54N
    Mf = (dtt/6) * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN); % Seasonally varying non-rain melt factor
    
    % New snow temperature and heat deficit from new snow
    if Ta < 0
        T_snow_new = Ta;
    else
        T_snow_new = 0;
    end
    
    % Change in the heat deficit due to new snowfall (mm), 80 cal/g: latent
    % heat of fusion, 0.5 cal/g/C: specific heat of ice
    delta_HD_snow = -(T_snow_new*Pn)/(80/0.5);   

    % Change in heat deficit due to a temperature gradient
    delta_HD_T = NMF * (dtp/6) * (Mf/MFMAX) * (ATI - T_snow_new);
    
    % Update ATI(Antecedent Temperature Index)
    if Pn > (1.5*dtp)
        ATI = T_snow_new;      %Antecedent temperature index
    else
        TIPM_dtt = 1.0 - ((1.0 - TIPM)^(dtt/6));
        ATI = ATI + TIPM_dtt * (Ta - ATI);
    end
    ATI = min(ATI,0);
   
    
    %                              SNOW MELT 
    % ---------------------------------------------------------------------                            
    T_rain = max(Ta,0);   % Temperature of rain (deg C), Ta or 0C, whichever greater
    if RAIN > (0.25 * dtp)  % Rain-on-Snow Melt
        stefan = 6.12*(10^(-10));  % Stefan-Boltzman constant (mm/K/hr)
        e_sat  = 2.7489*(10^8)*exp((-4278.63/(Ta+242.792)));  % Saturated vapor pressure at Ta (mb)
        % P_atm: Atmospheric pressure (mb) where elevation is in HUNDREDS
        % of meters (this is incorrectly stated in the manual)
        P_atm  = 33.86*(29.9-(0.335*(IN_ELEV/100))+(0.00022*((IN_ELEV/100)^2.4)));   
        term1 = stefan * dtp * (((Ta+273)^4)-(273^4));
        term2 = 0.0125 * RAIN * T_rain;
        term3 = 8.5 * UADJ * (dtp/6) * ((0.9*e_sat - 6.11) + (0.00057*P_atm*Ta));
        Melt = term1 + term2 + term3;
        Melt = max(Melt,0);
    elseif RAIN <= (0.25 * dtp) && (Ta > MBASE) % Non-Rain Melt
        Melt = (Mf * (Ta - MBASE) * (dtp/dtt)) + (0.0125 * RAIN * T_rain);
        Melt = max(Melt,0);
    else
        Melt = 0;
    end
   
    
    %                    RIPENNESS OF THE SNOW COVER 
    % ---------------------------------------------------------------------
    % E   : Excess liquid water in the snow cover                
    % W_i : water equivalent of the ice portion of the snow cover
    % W_q : liquide water held by the snow
    % W_qx: liquid water storage capacity
    % Qw  : Amount of available water due to melt and rain
    
    Deficit = max(Deficit + delta_HD_snow + delta_HD_T, 0);   % Deficit = heat deficit (mm)

    if Melt < W_i
        
        W_i = W_i-Melt;
        Qw = Melt + RAIN;
        W_qx = PLWHC * W_i;
        
        if Deficit>(0.33*W_i) % limits of heat deficit
            Deficit=0.33*W_i;
        end
        
        if (Qw + W_q) > (Deficit + Deficit*PLWHC + W_qx)
            % THEN the snow is RIPE
            E = Qw + W_q - W_qx - Deficit-(Deficit*PLWHC);  % Excess liquid water (mm)
            W_i = W_i + Deficit;  % W_i increases because water refreezes as heat deficit is decreased
            W_q = W_qx;  % fills liquid water capacity
            Deficit = 0;
            
        elseif (Qw >= Deficit) && ((Qw + W_q) <= ((Deficit*(1+PLWHC)) + W_qx))
            % THEN the snow is NOT yet ripe, but ice is being melted
            E = 0;
            W_i = W_i + Deficit; % W_i increases because water refreezes as heat deficit is decreased 
            W_q = W_q + Qw - Deficit;
            Deficit = 0;
            
        else
            % THEN the snow is NOT yet ripe
            E = 0;
            W_i = W_i + Qw; % W_i increases because water refreezes as heat deficit is decreased
            Deficit = Deficit - Qw;
        end
        
    else % Melt >= W_i
        
        Melt=W_i+W_q;
        W_i=0;
        W_q=0;
        Qw = Melt + RAIN;
        E = Qw;
        
    end
    
    if Deficit == 0
        ATI = 0;
    end
     
    
    %                          CONSTANT RELEASE 
    % ---------------------------------------------------------------------                                  
    % Constant daily amount of melt which takes place at the snow-soil
    % interface whenever there is a snow cover
    if W_i > DAYGM
        gmwlos = (DAYGM/W_i)*W_q;
        gmslos = DAYGM;
        gmro = gmwlos + gmslos;
        W_i = W_i - gmslos;
        W_q = W_q - gmwlos;
        
        E = E + gmro;
        SWE = W_i + W_q;
        
    else
        gmro = W_i+W_q;
        W_i=0;
        W_q=0;
        
        E = E + gmro;
        SWE = 0;
    end
       
    
    %                         UPDATE OUTPUTS
    % ---------------------------------------------------------------------
    OUT_EP(t) =  E;       
    OUT_SWE(t)     =  SWE; 
    OUT_MELT(t)    =  Melt;
    OUT_ADJPR(t)        =  Pn+RAIN;

end
    
% UPDATE INITIAL STATES
% ---------------------
OUT_INIT(1) = W_i;
OUT_INIT(2) = ATI;
OUT_INIT(3) = W_q;
OUT_INIT(4) = Deficit;
      