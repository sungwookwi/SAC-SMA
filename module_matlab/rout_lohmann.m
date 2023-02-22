function [OUT_TOTAL,OUT_BASE] = rout_lohmann(IN_DIRECT,IN_BASE,IN_FLOWLEN,IN_PAR,IN_OUTLETIND)
%#codegen
%   Functin for river channel routing model coupled with soil moisture
%   accouting module based on Lohmann routing model
%   - HRU Unit Hydrograph is represented by the Gamma distribution
%   - River channel routing is based on the linearized Saint-Venant Eq
% 
%   INPUTS
%	    IN_DIRECT       Direct runoff from hydrologic response unit (HRU) simulated by SMA module
%	    IN_BASE         Base runoff from hydrologic response unit (HRU) simulated by SMA module
%       IN_FLOWLEN      Travel distance of runoff from HRU outlet to basin outlet
%	    IN_PAR          Unit hydrograph & Lohmann routing parameters
%	    IN_OUTLETIND    Watershed outlet indicator
% 
%   OUTPUTS
%       OUT_TOTAL   Total basin streamflow
%       OUT_BASE    Baseflow component of total basin streamflow
%	
%   Model Parameter Description
%	    IN_PAR(1)   N       HRU unit hydrograph shape parameter
%	    IN_PAR(2)   K       HRU unit hydrograph scale parameter
%	    IN_PAR(3)   VELO    Wave velocity in the linearized Saint-Venant equation(m/s)
%	    IN_PAR(4)   DIFF    Diffusivity in the linearized Saint-Venant equation(m2/s)
% 
%   The directive %#codegen indicates that this function is intended for C/C++ code generation
%   codegen syntax to generate MEX function:
%	cfg = coder.config('mex');
%   cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%   codegen rout_lohmann -args {coder.typeof(0,[100000 1],[1 0]),coder.typeof(0,[100000 1],[1 0]),0,zeros(1,4),0} -config cfg -report
% 
% 
%--------------------------------------------------------------------------
% 
%	Author: Sungwook Wi
%	e-mail: sw2275@cornell.edu
%__________________________________________________________________________


%------------------------ 
%	Loading Parameters
%------------------------ 
N      = IN_PAR(1);
K      = IN_PAR(2);
VELO   = IN_PAR(3);
DIFF   = IN_PAR(4);


%--------------------------------------------------------------------------
%          Unit Hydrograph Base Time for HRU & channel routing 
% -------------------------------------------------------------------------
KE      =  12;          % Base time for HRU UH (day)
UH_DAY  =  96;          % Base time for river routing UH (day)
DT      =  3600;        % Time step in second for solving Saint-Venant equation. This will affect TMAX
TMAX    =  UH_DAY * 24; % Base time of river routing UH in hour because DT is for an hour
LE      =  48*50;       % Base time (hr) for Green function values


% -------------------------------------------------------------------------
%       Derive Daily River Impulse Response Function(Green's function)
% -------------------------------------------------------------------------
UH_river = zeros(1,UH_DAY);
if IN_OUTLETIND == 1 % if the HUR contains the watershed outlet
    
    UH_river(1) =  1;
    
else % if not, calculate Green's function to solve the linearized Saint-Venant Eq
    
    t = 0;
    uhm_grid = zeros(LE,1);
    for k = 1:LE
        t = t + DT;
        
        pot = ((VELO*t-IN_FLOWLEN)^2)/(4*DIFF*t);
        if pot <= 69
            H = IN_FLOWLEN./(2*t.*sqrt(pi*t*DIFF)).*exp(-pot);
        else
            H = 0.0;
        end
        uhm_grid(k) = H;
        
    end
    
    if sum(uhm_grid) == 0
        uhm_grid(1) =  1.0;
    else
        uhm_grid = uhm_grid/sum(uhm_grid);
    end
        
    UHM      = uhm_grid;
    
    FR  = zeros(TMAX,2);
    FR(1:24,1) = 1/24;
    

    for t = 1:TMAX
        for L = 1:TMAX+24
            if t-L > 0
                FR(t,2) = FR(t,2) + FR(t-L,1) * UHM(L);
            end
        end
    end
    
    for t = 1:UH_DAY
        UH_river(t) = sum(FR(24*t-23:24*t,2));
    end
        
end


%--------------------------------------------------------------------------
%                HRU UH represented by Gamma distribution
%--------------------------------------------------------------------------
UH_HRU_direct = zeros(KE,1);
for i = 1:KE
    x = linspace(24*(i-1),24*i,1001);
    hruh_fun = 1/(1/K)/gamma(N)*(x/(1/K)).^(N-1).*exp(-x/(1/K));
    pinteg = sum(hruh_fun)*(x(2)-x(1));
    UH_HRU_direct(i) = pinteg;
end
UH_HRU_base    = zeros(KE,1);
UH_HRU_base(1) = 1;


%--------------------------------------------------------------------------
%         Combined UH for HRU's response at the watershed outlet 
%--------------------------------------------------------------------------
UH_direct = zeros(1,KE+UH_DAY-1);
UH_base   = zeros(1,KE+UH_DAY-1);
for k = 1:KE
    for u = 1:UH_DAY
        UH_direct(k+u-1) = UH_direct(k+u-1) + UH_HRU_direct(k) * UH_river(u);
        UH_base(k+u-1)   = UH_base(k+u-1) + UH_HRU_base(k) * UH_river(u);
    end
end
UH_direct = UH_direct/sum(UH_direct);
UH_base   = UH_base/sum(UH_base);


%--------------------------------------------------------------------------
%             Make Convolution for watershed outlet total flow
%--------------------------------------------------------------------------
directflow = zeros(length(IN_DIRECT),1);
OUT_BASE   = zeros(length(IN_DIRECT),1);
for i = 1:length(IN_DIRECT)
    for j = 1:KE+UH_DAY-1
        if i-j+1 >= 1
            directflow(i) = directflow(i) + UH_direct(j) * IN_DIRECT(i-j+1);
            OUT_BASE(i)   = OUT_BASE(i) + UH_base(j) * IN_BASE(i-j+1);
        end
    end
end
OUT_TOTAL = directflow + OUT_BASE;


