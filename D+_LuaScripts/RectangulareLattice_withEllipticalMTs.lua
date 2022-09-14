Information = {
	Name = "Rectangular Close-Packed Elliptical MTs",
	Type = "Symmetry",
	NLP = 1,
	MinLayers = 8,
	MaxLayers = 8,
};

function Populate(p, nlayers)	 
	if (p == nil or nlayers ~= 8 or table.getn(p[1]) ~= 1) then				
		error("Parameter matrix must be 7x1, it is " .. nlayers .. "x" .. table.getn(p[1]));
	end
	
	--a 			= p[1][1]; --major axis
	b 			= p[1][1]; --minor axis
	npoints_new = p[2][1];
	x_repeats   = p[3][1];
	y_repeats   = p[4][1];
	extra_row   = p[5][1]; 
	noise 		= p[6][1];
	square 		= p[7][1];
	L_wall 		= p[8][1];
		
	res 		= {};
	ind 		= 1;	
	b_new 		= b;

	--
	
--Create the Array
x_shift = -(x_repeats-1)/2;
for rep_x = 0,(x_repeats-1) do	
	y_shift = -(y_repeats-1)/2;
	if square == 0 then	y_shift = y_shift + (rep_x % 2)*(b+L_wall) end
	for rep_y = 0,(y_repeats-1) do	
		y_shift = y_shift + (b_new+L_wall);
		noise_minor = noise*(math.random()-0.5);
		b_new = b+noise_minor
		ellip_x, ellip_y, a = CreateEllipse(b_new,npoints_new);
		for ii = 1,(npoints_new) do
			res[ind] = {x_shift+ellip_x[ii],y_shift+ellip_y[ii],0,0,0,0};
			ind = ind + 1;
		end
		y_shift = y_shift + (b_new+L_wall);
	end
	if square == 0 then
		x_shift = x_shift + 2*(a+L_wall)*math.sqrt(3)/2;
	else
		x_shift = x_shift + 2*(a+L_wall);
	end
end
--One last loop for our extra row	
y_shift = -(y_repeats-1)/2;
if square == 0 then	y_shift = y_shift + ((x_repeats) % 2)*(b+L_wall) end
for rep_y = 0,(extra_row-1) do	
	y_shift = y_shift + (b_new+L_wall);
	noise_minor = noise*(math.random()-0.5);
	b_new = b+noise_minor
	ellip_x, ellip_y, a = CreateEllipse(b_new,npoints_new);
	for ii = 1,(npoints_new) do
	    res[ind] = {x_shift+ellip_x[ii],y_shift+ellip_y[ii],0,0,0,0};
		ind = ind + 1;
	end
	y_shift = y_shift + (b_new+L_wall);
end	

	return res;
end

-- Optional display parameters
function GetLayerName(index)
	if index == 0 then
		return "Minor Axis";
	elseif index == 1 then
		return "# Points";
	elseif index == 2 then
		return "X Repeats";
	elseif index == 3 then
		return "Y Repeats";
	elseif index == 4 then
		return "Extra Row";
	elseif index == 5 then
		return "Noise";
	elseif index == 6 then
		return "Square?";
	elseif index == 7 then
		return "L_Wall";
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
		return 7.5;
	elseif layer == 1 then
		return 13;
	elseif layer == 2 then
		return 1;
	elseif layer == 3 then
		return 1;
	elseif layer == 4 then
		return 4;
	elseif layer == 5 then
		return 0;
	elseif layer == 6 then
		return 0;
	elseif layer == 7 then
		return 2.45;
	end
end

function CreateEllipse(b,npoints_new)
a = (1/30) *(math.sqrt(-500*math.pow(b,2)+7140*b+42483) - 20*b + 357 ); --Ramanujan simple approximation of an ellipse with same perimeter as a circle of radius 11.9
npoints		= 1000;

delta_theta=2.0*math.pi/npoints;
theta = {}; table.insert(theta,0.0);
delta_s = {}; table.insert(delta_s,0.0);
integ_delta_s = {}; table.insert(integ_delta_s,0.0);
integ_delta_s_val = 0;--{};
for iTheta = 0,(npoints) do
    delta_s_val=math.sqrt(math.pow(a,2) * math.pow(math.sin(iTheta*delta_theta),2) + math.pow(b,2)*math.pow(math.cos(iTheta*delta_theta),2));
    table.insert(theta,iTheta*delta_theta);
    integ_delta_s_val = integ_delta_s_val+delta_s_val*delta_theta;
    table.insert(integ_delta_s,integ_delta_s_val);
end
integ_delta_s_norm = {}
for iEntry = 1,(npoints+1) do
    table.insert(integ_delta_s_norm,integ_delta_s[iEntry]/integ_delta_s[npoints]*2.0*math.pi);--integ_delta_s_norm.append(iEntry/integ_delta_s[-1]*2.0*math.pi)   
end	
ellip_x={}; ellip_y={};
delta_theta_new=2*math.pi/npoints_new;
for theta_index = 0,(npoints_new-1) do -- do in range(npoints_new):
    theta_val = theta_index*delta_theta_new
    for lookup_index = 1,(npoints+1) do --in range(len(integ_delta_s_norm)):
        if ( (theta_val >= integ_delta_s_norm[lookup_index]) and (theta_val < integ_delta_s_norm[lookup_index+1]) ) then
            theta_prime = theta[lookup_index];
			break;
        end
	end
	ellip_x[theta_index+1] = a*math.cos(theta_prime);
	ellip_y[theta_index+1] = b*math.sin(theta_prime);
end

return ellip_x, ellip_y, a

end


