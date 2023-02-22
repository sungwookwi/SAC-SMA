function OUT_PET = pet_hamon(IN_SDATE, IN_EDATE, IN_TAVG, IN_LAT, IN_PAR)
%#codegen
%   Function for calculating potential evapotranspiration based on Hamon method
%
%   INPUTS
%       IN_SDATE    Start Date of Simulation [year, month, day]
%       IN_EDATE    End Date of Simulation [year month day]
%       IN_TAVG     Time series of daily average temperature(C) for the specified simulation period
%       IN_LAT      Latitude of the area
%       IN_PAR      Hamon model parameter
% 
%   OUTPUTS
%       OUT_PET     Potential evapotranspiration
% 
%   Model Parameter Description
%       IN_PAR(1)   coeff   PET adjustment factor
%
%   The directive %#codegen indicates that this function is intended for C/C++ code generation
%   codegen syntax to generate MEX function:
%       cfg = coder.config('mex');
%       cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
%       codegen pet_hamon -args {zeros(1,3),zeros(1,3),coder.typeof(0,[200000 1],[1 0]),0,0} -config cfg -report
%
% 
%--------------------------------------------------------------------------
% 
%   References
%       Lu Jianbiao, Ge Sun, Steven G. McNulty, Devendra M. Amatya (2005), A
%       comparison of six potential evaportranspiration methods for
%       regional use in the southeastern United States. Journal of the
%       American Water Resources Association, Vol. 41, No. 3., pp. 621-633,
%       doi:10.1111/j.1752-1688.2005.tb03759.x

%       Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahla,
%       Hsini Wua, Robert M. Schoolfieldb (1995), A model comparison for
%       daylength as a function of latitude and day of year Ecological
%       Modelling, Volume 80, Issue 1, Pages 87–95
% 
%--------------------------------------------------------------------------
% 
%	Author: Sungwook Wi
%	e-mail: sw2275@cornell.edu
%__________________________________________________________________________


%------------------------ 
%	Loading Parameters
%------------------------ 
coeff   =   IN_PAR;


% -------------------------------------------------------------------------
%                 Daylight hour calculation based on CBM model 
% -------------------------------------------------------------------------
daylighthr_mat = zeros(length(IN_TAVG),1);
yr = (IN_SDATE(1):IN_EDATE(1))';
n  = 1;
for iy = 1:length(yr)
    
    % Julian date generation 
    % ----------------------
    if iy == 1
        if ( IN_SDATE(2) <= 2 ) % January & February
            year  = IN_SDATE(1) - 1.0;
            month = IN_SDATE(2) + 12.0;
        else
            year  = IN_SDATE(1);
            month = IN_SDATE(2);
        end
        day  = IN_SDATE(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year  = IN_SDATE(1) - 1.0;
        month = 1 + 12.0;
        day   = 1;
        jd_2  = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
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
            year  = IN_EDATE(1);
            month = IN_EDATE(2);
        end
        day  = IN_EDATE(3);
        jd_1 = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
            floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5;
        
        year  = IN_EDATE(1) - 1.0;
        month = 1 + 12.0;
        day   = 1;
        jd_2  = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
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
    % ----------------------
    
    var_theta   = 0.2163108 + 2*atan(0.9671396*tan(0.0086*(jdate-186)));
    var_pi      = asin(0.39795*cos(var_theta));
    daylighthr  = 24-24/pi*acos((sin(0.8333*pi/180)+sin(IN_LAT*pi/180)*sin(var_pi))./(cos(IN_LAT*pi/180)*cos(var_pi)));
    daylighthr_mat(n:n+length(jdate)-1) = daylighthr';
    
    n = n + length(jdate);
end


% -------------------------------------------------------------------------
%       Computing Potential Evapotranspiration based on Hamon euqation
% -------------------------------------------------------------------------
esat = 0.611*exp(17.27*IN_TAVG./(237.3+IN_TAVG));
OUT_PET  = coeff*29.8*daylighthr_mat.*(esat./(IN_TAVG+273.2));

