Information = {
	Name = "Elliptical MT",
	Type = "Symmetry",
	NLP = 1,
	MinLayers = 3,
	MaxLayers = 3,
};

function Populate(p, nlayers)	 
	if (p == nil or nlayers ~= 3 or table.getn(p[1]) ~= 1) then				
		error("Parameter matrix must be 7x1, it is " .. nlayers .. "x" .. table.getn(p[1]));
	end
	
	a 			= p[1][1];
	b 			= p[2][1];
	npoints		= 1000;
	npoints_new = p[3][1];
	res = {};
	ind = 1;	
	
delta_theta=2.0*math.pi/npoints;

theta = {};
table.insert(theta,0.0);
delta_s = {};
table.insert(delta_s,0.0);
integ_delta_s = {};
table.insert(integ_delta_s,0.0);
integ_delta_s_val = 0;--{};
--table.insert(integ_delta_s_val,0.0); --# integrated probability density

for iTheta = 0,(npoints) do
    --# ds/d(theta):
    delta_s_val=math.sqrt(math.pow(a,2) * math.pow(math.sin(iTheta*delta_theta),2) + math.pow(b,2)*math.pow(math.cos(iTheta*delta_theta),2));
    table.insert(theta,iTheta*delta_theta); --theta.append(iTheta*delta_theta)
    --delta_s.append(delta_s_val)
    --# do integral
    integ_delta_s_val = integ_delta_s_val+delta_s_val*delta_theta;
    table.insert(integ_delta_s,integ_delta_s_val); --integ_delta_s.append(integ_delta_s_val)
end
    
--# normalize integrated ds/d(theta) to make into a scaled CDF (scaled to 2*pi)
integ_delta_s_norm = {}
for iEntry = 1,(npoints+1) do
    table.insert(integ_delta_s_norm,integ_delta_s[iEntry]/integ_delta_s[npoints]*2.0*math.pi);--integ_delta_s_norm.append(iEntry/integ_delta_s[-1]*2.0*math.pi)   
end	

--# Create corrected ellipse using lookup function
ellip_x_prime={};
ellip_y_prime={};

delta_theta_new=2*math.pi/npoints_new;

for theta_index = 0,(npoints_new-1) do -- do in range(npoints_new):
    theta_val = theta_index*delta_theta_new
    
--# Do lookup:
    for lookup_index = 1,(npoints+1) do --in range(len(integ_delta_s_norm)):
        if ( (theta_val >= integ_delta_s_norm[lookup_index]) and (theta_val < integ_delta_s_norm[lookup_index+1]) ) then
            theta_prime = theta[lookup_index];
			break;
        end
	end
	res[ind] = {a*math.cos(theta_prime),b*math.sin(theta_prime),0,0,0,0};
	ind = ind + 1;
end

	return res;
end

-----------------------------------------------------

-- Optional display parameters
function GetLayerName(index)
	if index == 0 then
		return "L_x";
	elseif index == 1 then
		return "L_y";
	elseif index == 2 then
		return "# Points";
	else	
		return "N/A";
	end
end

function GetLayerParameterName(index)
	if index == 0 then
		return "Parameter";
	else
		return "N/A"
	end
end
	
function IsParamApplicable(layer, layerParam)
	return true;
end

function GetDefaultValue(layer, layerParam)
	if layer == 0 then
		return 15.585;
	elseif layer == 1 then
		return 7.5;
	elseif layer == 2 then
		return 50;
	end
end

	
------ Function iterators to use in "for" loops
function list_iter (t)
    local i = 0
    local n = table.getn(t)
    return function ()
        i = i + 1
        if i <= n then return t[i] end
    end
end





