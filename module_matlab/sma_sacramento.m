function [OUT_SURF,OUT_BASE,OUT_GW,OUT_AET,OUT_UZTWC,OUT_UZFWC,OUT_LZTWC,OUT_LZFPC,OUT_LZFSC,OUT_ADIMC]=sma_sacramento(IN_PET,IN_EP,IN_PAR,IN_INIT)
%#codegen
%   Function for simulating soil moisture accounting based on the
%   Sacramento model written in Fortran originally
%
%   INPUTS
%       IN_PET      Time series of potential evapotranspiration time series (L-by-1 or 1-by-L)     
%       IN_EP       Time series of effective precipitation (L-by-1 or 1-by-L): the amount of precipitation added in the soil
%       IN_PAR      Sacramento model parameters (16-by-1 or 1-by-16)
%       IN_INIT     Storage initial state (6-by-1 or 1-by-6)
%
%   OUTPUTS (L-by-1)
%       OUT_SURF	Surface runoff (direct runoff + overflow + interflow)
%       OUT_BASE	Baseflow
%       OUT_AET     Actual evapotranspiration
%       OUT_UZTWC   Upper zone tension water storage
%       OUT_UZFWC   Upper zone free water storage
%       OUT_LZTWC   Lower zone tension water storage
%       OUT_LZFPC   Lower zone primary free water storage
%       OUT_LZFSC   Lower zone supplementary free water storage
%       OUT_ADIMC   Additional impervious storage
%
%   Model Parameter Description
%       IN_PAR(1)	uztwm	Upper zone tension water storage maximum [mm]
%       IN_PAR(2)	uzfwm	Upper zone free water storage maximum [mm]
%       IN_PAR(3)	lztwm	Lower zone tension water storage maximum [mm]
%       IN_PAR(4)	lzfpm	Lower zone primary free water storage maximum [mm]
%       IN_PAR(5)	lzfsm	Lower zone supplementary free water storage maximum [mm]
%       IN_PAR(6)	uzk     Upper zone free water lateral depletion rate (btw 0 and 1)
%       IN_PAR(7)	lzpk	Lower zone primary free water depletion rate (btw 0 and 1)
%       IN_PAR(8)	lzsk	Lower zone supplementary free water depletion rate (btw 0 and 1)
%       IN_PAR(9)	zperc	Percolation demand scale parameter: multiplier of the percolation equation
%       IN_PAR(10)	rexp	Percolation demand shape parameter: exponent of the percolation equation
%       IN_PAR(11)	pfree	Percolating water split parameter (btw 0 and 1): fraction of water percolating to lower free water storage
%       IN_PAR(12)	pctim	Impervious fraction of the watershed area (btw 0 and 1): permanent impervious area
%       IN_PAR(13)	adimp	Additional impervious areas (btw 0 and 1): temporary impervious area
%       IN_PAR(14)	riva	Fraction of riparian vegetation area (btw 0 and 1)
%       IN_PAR(15)	side	The ratio of deep recharge to river channel base flow
%       IN_PAR(16)	rserv	Fraction of lower zone free water not transferrable to lower zone tension water (btw 0 and 1)
%
%   Initial Storage Condition
%       IN_INIT(1)  uztwc   Initial upper zone tension water content [mm]
%       IN_INIT(2)  uzfwc   Initial upper zone free water content [mm]
%       IN_INIT(3)  lztwc   Initial lower zone tension water content [mm]
%       IN_INIT(4)  lzfpc   Initial lower zone primary free water content [mm]
%       IN_INIT(5)  lzfsc   Initial lower zone supplementary free water content [mm]
%       IN_INIT(6)  adimc   Initial additional impervious area water content [mm]
%
%   The directive %#codegen indicates that this function is intended for C/C++ code generation
%   codegen syntax to generate MEX function:
%       cfg = coder.config('mex');
%       cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%       codegen sma_sacramento -args {coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),zeros(1,16),zeros(1,6)} -config cfg -report
%
% 
%--------------------------------------------------------------------------
% 
%	Author: Sungwook Wi
%	e-mail: sw2275@cornell.edu
%__________________________________________________________________________ 


%---------------------------
%	Initiailizing Outputs
%---------------------------
L           =   length(IN_PET); % length of the input time series
OUT_SURF    =   nan(L,1);  
OUT_BASE    =   nan(L,1);
OUT_GW      =   nan(L,1);
OUT_AET     =   nan(L,1);  
OUT_UZTWC   =   nan(L,1);
OUT_UZFWC   =   nan(L,1);
OUT_LZTWC   =   nan(L,1);
OUT_LZFPC   =   nan(L,1);
OUT_LZFSC   =   nan(L,1);
OUT_ADIMC   =   nan(L,1);


%------------------------ 
%	Loading Parameters
%------------------------ 
% Storage capacity
uztwm   =   IN_PAR(1);
uzfwm   =   IN_PAR(2);
lztwm   =   IN_PAR(3);
lzfpm   =   IN_PAR(4); 
lzfsm   =   IN_PAR(5); 
% Storage release coefficient
uzk     =   IN_PAR(6); 
lzpk    =   IN_PAR(7);  
lzsk    =   IN_PAR(8);  
% Percolation
zperc   =   IN_PAR(9);   
rexp    =   IN_PAR(10);   
pfree   =   IN_PAR(11);   
% Impervious area fraction
pctim   =   IN_PAR(12); 
adimp   =   IN_PAR(13); 
% Others
riva    =   IN_PAR(14);  
side    =   IN_PAR(15);  
rserv   =   IN_PAR(16);  


%--------------------------- 
%	Loading Initial State
%--------------------------- 
uztwc   =   IN_INIT(1);   
uzfwc   =   IN_INIT(2);   
lztwc   =   IN_INIT(3);   
lzfpc   =   IN_INIT(4); 
lzfsc   =   IN_INIT(5);  
adimc   =   IN_INIT(6);  


%------------------------ 
%	Execute Simulation
%------------------------ 
thres = 0.00001;	% if val < thres, val = 0 
parea = 1 - adimp - pctim;    % pervious area  
for i = 1:L
       
    % COMPUTE EVAPOTRANSPIRATION
    % --------------------------   
    edmnd = IN_PET(i);  % ET demand for the time interval
    
    % et1: ET from uztwc
    et1 = edmnd * uztwc/uztwm;
    red = edmnd - et1;  % residual ET demand
    uztwc = uztwc - et1;
    
    % et2: ET from uzfwc
    et2 = 0;
    if uztwc < 0
        et1 = et1 + uztwc; % et1 cannot exceed uztwc
        uztwc = 0;
        red = edmnd - et1;
        
        if uzfwc < red  
            et2 = uzfwc; 
            uzfwc = 0;
            red = red - et2;  
                            
        else 
            et2 = red;  
            uzfwc = uzfwc - et2;
            red = 0;
            
            % Transfer from uzfwc to uztwc
            if (uztwc / uztwm) < (uzfwc / uzfwm)
                uzrat = (uztwc + uzfwc) / (uztwm + uzfwm);
                uztwc = uztwm * uzrat;
                uzfwc = uzfwm * uzrat;
            end
                    
        end
        
    else % ET demand is all consumed at uztwc, therefore NO ET from uzfwc (et2=0)
        
        % Transfer from uzfwc to uztwc
        if (uztwc / uztwm) < (uzfwc / uzfwm) 
            uzrat = (uztwc + uzfwc) / (uztwm + uzfwm);
            uztwc = uztwm * uzrat;
            uzfwc = uzfwm * uzrat;
        end
           
    end
    if uztwc < thres; uztwc = 0; end
    if uzfwc < thres; uzfwc = 0; end
        
    % et3: ET from lztwc
    et3 = red * lztwc / (uztwm + lztwm); 
    lztwc = lztwc - et3;
    if lztwc < 0 % et3 cannot exceed lztwc
        et3 = et3 + lztwc; 
        lztwc = 0;        
    end
    ratlzt = lztwc / lztwm;
        
    % Resupply lztwc from lower zone free water if more water available there
    saved = rserv * (lzfpm + lzfsm);
    ratlz = (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved);
    if ratlzt < ratlz
        del = (ratlz - ratlzt) * lztwm;
        lztwc = lztwc + del;  % Transfer from lzfsc to lztwc
        lzfsc = lzfsc - del;
        if lzfsc < 0  % if tranfer exceeds lzfsc then remainder comes from lzfpc
            lzfpc = lzfpc + lzfsc;
            lzfsc = 0;
        end
    end
    if lztwc < thres; lztwc = 0; end
     
    % et5: ET from the additional impervious area
    et5 = et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm); 
    adimc = adimc - et5;
    if adimc < 0   
        et5 = et5 + adimc; % et5 cannot exceed adimc
        adimc = 0;
    end
    et5 = et5 * adimp;
    
    
    % COMPUTE PERCOLATION & RUNOFF
    % ----------------------------
    
    % Update uztwc & adimc with IN-EP
    % twx: available moisture in excess of uztwm requirements
    twx = IN_EP(i) + uztwc - uztwm;  
    if twx < 0 % all moisture held in uztw - no excess
        uztwc = uztwc + IN_EP(i);
        twx = 0;
    else    
        uztwc = uztwm;
    end
    
    adimc = adimc + IN_EP(i) - twx; % twx is added to adimc later in the for loop below
    

    % Initialize time interval sums
    sbf   = 0;  % Sum of baseflow(from lzfpc and lzfsc)
    ssur  = 0;  % Sum of surface runoff (from uzfwc and adimc)
    sif   = 0;  % Sum of interflow (from uzfwc)
    sperc = 0;  % Sum of percolation
    sdro  = 0;  % Sum of direct runoff (from adimc)
    
    
    % Determine computational time increments for the basic time interval
    ninc = floor(1.0 + 0.2*(uzfwc+twx));  % Number of time increments that the time interval is divided into for further soil-moisture accountng
    dinc = 1.0 / ninc;                    % Length of each increment in days
    pinc = twx / ninc;                    % Amount of available moisture for each increment
    
    % Compute free water depletion fractions for the time increment (basic depletions are for one day)
    duz   = 1 - (1 - uzk)^dinc;
    dlzp  = 1 - (1 - lzpk)^dinc;
    dlzs  = 1 - (1 - lzsk)^dinc;
    
    
    %--- Start incremental for-loop for the time interval
    for n = 1:ninc
        
        adsur = 0; % additional surface runoff from adimp
        
        % Compute direct runoff from adimp area
        ratio = (adimc - uztwc) / lztwm;
        if ratio < 0; ratio = 0; end
        addro = pinc*(ratio^2); % direct runoff from adimp
        
        % Compute baseflow and keep track of time interval sum.
        bf_p = lzfpc * dlzp; % Baseflow from lzfpc
        lzfpc = lzfpc - bf_p;
        if lzfpc <= 0.0001
            bf_p = bf_p + lzfpc;
            lzfpc = 0;
        end
        sbf = sbf + bf_p;
        
        bf_s = lzfsc * dlzs;  % Baseflow from lzfsc
        lzfsc = lzfsc - bf_s;
        if lzfsc <= 0.0001
            bf_s = bf_s + lzfsc;
            lzfsc = 0;
        end
        sbf = sbf + bf_s; % Total Baseflow from lzfpc & lzfsc
        
        % Compute PERCOLATION
        if (pinc + uzfwc) <= 0.01 % if no water available then skip.
            
            uzfwc = uzfwc + pinc;
            
        else
            
            percm = lzfpm * dlzp + lzfsm * dlzs;
            perc = percm * uzfwc / uzfwm;
            
            % defr: lower zone moisture deficiency ratio
            defr = 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm); 
            
            if defr < 0; defr = 0; end
            
            perc = perc * (1.0 + zperc * (defr^rexp));            
            if perc >= uzfwc     % percolation cannot exceeds uzfwc
                perc = uzfwc;
            end
            
            uzfwc = uzfwc - perc;    % update uzfwc with percolation
            
            % Check to see if percolation exceeds lower zone deficiency.
            check = lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm;
            if check > 0    
                perc = perc - check;
                uzfwc = uzfwc + check;
            end
            % sperc: summation of the time interval percolation
            sperc = sperc + perc;
            
            
            % Compute INTERFLOW and keep track of time interval sum
            % Note PINC has not yet been added
            del = uzfwc * duz; % The amount of interflow
            sif = sif + del;
            uzfwc = uzfwc - del;
            
            % Distribute percolated water into the lower zones
            % Tension water must be filled first except for the pfree area
            % perct is percolation to tension water and percf is
            % percolation going to free water
            
            perct = perc * (1.0 - pfree);  
            if (perct + lztwc) <= lztwm
                lztwc = lztwc + perct;
                percf = 0.0; 
            else 
                percf = lztwc + perct - lztwm; % percolation in excess of tension requirements
                lztwc = lztwm;
            end
            
            % Distribute percolation among the free water storages.
            percf = percf + (perc * pfree);
            if percf ~= 0
                % hpl: relative size of the primary storage to total lower zone free water storage
                hpl = lzfpm / (lzfpm + lzfsm); 
                
                % Relative fullness of each storage.
                ratlp = lzfpc / lzfpm;
                ratls = lzfsc / lzfsm;
                
                % fracp: fraction going to primary
                fracp = (hpl * 2.0 * (1.0 - ratlp)) / ((1.0 - ratlp) + (1.0 - ratls)); 
                if fracp > 1.0; fracp = 1.0; end
                
                percp = percf * fracp; % Amount of the excess percolation going to primary
                percs = percf - percp; % Amount of the excess percolation going to supplemental
                lzfsc = lzfsc + percs;
                if lzfsc > lzfsm
                    percs = percs - lzfsc + lzfsm;
                    lzfsc = lzfsm;
                end
                lzfpc = lzfpc + (percf - percs);
                
                % Check to make sure lzfpc does not exceed lzfpm
                if lzfpc >= lzfpm  
                    excess = lzfpc - lzfpm;
                    lztwc = lztwc + excess;
                    lzfpc = lzfpm;
                end
                
            end
            
            
            % Distribute PINC between uzfws and surface runoff
            if pinc ~= 0
                
                % check if pinc exceeds uzfwm
                if (pinc + uzfwc) <= uzfwm  
                    % no surface runoff
                    uzfwc = uzfwc + pinc;  
                else 
                    % compute surface runoff (sur) and keep track of time
                    % interval sum
                    sur = pinc + uzfwc - uzfwm; 
                    uzfwc = uzfwm;
                    ssur = ssur + (sur * parea);
                    
                    % adsur is the amount of surface runoff which comes from
                    % that portion of adimp which is not currently generating
                    % direct runoff. addro/pinc is the fraction of adimp
                    % currently generating direct runoff.
                    adsur = sur * (1.0 - addro / pinc);
                    ssur = ssur + adsur * adimp;
                end
            end
            
            
        end
        
        % adimp area water balance
        adimc = adimc + pinc - addro - adsur;
        if adimc > (uztwm + lztwm)
            addro = addro + adimc - (uztwm + lztwm);
            adimc = uztwm + lztwm;
        end
        
        sdro  = sdro + (addro * adimp); % direct runoff from adimp
        
        if adimc < thres;  adimc = 0; end
        
    end
    %--- END of incremental for loop
    
    
    % Compute sums and adjust runoff amounts by the area over which they are generated.
    
    % eused: the ET from parea (1.0-adimp-pctim)
    eused = et1 + et2 + et3; % not adjusted yet
    sif = sif * parea; % adjusted interflow
    
    % Separate channel component of baseflow from the non-channel component
    tbf = sbf * parea;  % tbf: total baseflow
    bfcc = tbf * (1.0 / (1.0 + side));    % baseflow, channel component
    bfncc = tbf - bfcc; % baseflow, non-channel component (deep groundwater recharge)
    
    
    % Total channel inflow divided into Surface and Base flow
    % tci = roimp + sdro + ssur + sif + bfcc;
    roimp = IN_EP(i) * pctim; % direct runoff from impervious area
    surf = roimp + sdro + ssur + sif;
    base = bfcc; 
    gw = bfncc;
        
    % et4: ET from riparian vegetation.
    et4 = (edmnd - eused) * riva;
       
    % Total inflow to channel for the timestep
    tci = surf + base - et4;
    
    if tci <= 0 % et4, surf, base need to be updated
        et4 = surf + base;
        surf = 0;
        base = 0;
        
    else % surf & base need to be updated
        surf = surf - et4/2;
        base = base - et4/2;
        if surf < 0
            base = base + surf;
            surf = 0;
        end
        if base < 0
            surf = surf + base;
            base = 0;
        end
    end
    
    % Compute total AET
    eused = eused * parea;
    aet = eused + et4 + et5;
    
    % Check if adims >= uztws
    if adimc < uztwc
        adimc = uztwc;
    end
    
    
    % UPDATE OUTPUTS
    % --------------
    OUT_SURF(i)     = surf;
    OUT_BASE(i)     = base;
    OUT_GW(i)       = gw;
    OUT_AET(i)      = aet;
    OUT_UZTWC(i)    = uztwc;
    OUT_UZFWC(i)    = uzfwc;
    OUT_LZTWC(i)    = lztwc;
    OUT_LZFPC(i)    = lzfpc;
    OUT_LZFSC(i)    = lzfsc;
    OUT_ADIMC(i)    = adimc;
    
end