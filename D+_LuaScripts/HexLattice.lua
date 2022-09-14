Information = {
	Name = "Hex Lattice of MTs with noise",
	Type = "Symmetry",
	NLP = 1,
	MinLayers = 6,
	MaxLayers = 6,
};

function Populate(p, nlayers)	 
	if (p == nil or nlayers ~= 6 or table.getn(p[1]) ~= 1) then				
		error("Parameter matrix must be 2x1, it is " .. nlayers .. "x" .. table.getn(p[1]));
	end
	
	L_a = p[1][1];
	L_extra = p[2][1];
	repeats = p[3][1];
	protofilaments		= p[4][1];
	noise = p[5][1];
	roughen_percent = p[6][1];
	
	angle_spacing		= 2.0 * math.pi / protofilaments;
	
	res = {};
	ind = 1;	
	L_tot = L_a + L_extra; -- this is your center to center distance

	
		y_shift = -repeats*(L_tot)*math.sqrt(3)/2;
		for rep_y = -repeats,(repeats) do		
		x_shift = -repeats*(L_tot);
		x_shift = x_shift + (rep_y % 2)*(L_tot)/2;
		for rep_x = -repeats,(repeats) do	
		
			if (((rep_y == -repeats) or (rep_y == repeats) or (rep_x == -repeats) or (rep_x == repeats)) and (math.random() < roughen_percent)) then
			--do nothin
			else
				for theta = 1, protofilaments do 
					
					x		= L_a/2 * math.sin(theta*angle_spacing);
					y		= L_a/2 * math.cos(theta*angle_spacing);

					res[ind] = {x_shift+x,y_shift+y,0,0,0,0};
					ind = ind + 1;
				end
			end
			
		x_shift = x_shift + (L_tot) + noise*(math.random()-0.5);
		end

		y_shift = y_shift + (L_tot)*math.sqrt(3)/2 + noise*(math.random()-0.5);
		end
	
	return res;


end

-----------------------------------------------------
-- UI
function GetLayerName(index)
	if index == 0 then
		return "MT Diameter";
	elseif index == 1 then
		return "Extra space between MTs";
	elseif index == 2 then
		return "Lattice repeats";
	elseif index == 3 then
		return "# protofilaments";
	elseif index == 4 then
		return "Spatial Noise Max Displacement";
	elseif index == 5 then
		return "Percent Missing from Outer Layer";
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
		return 11.9;
	elseif layer == 1 then
		return 4.9;--6;
	elseif layer == 2 then
		return 2;
	elseif layer == 3 then
		return 50;
	elseif layer == 4 then
		return 0;
	elseif layer == 5 then
		return 0;
	end
end

	