function down_cpf_lun
data_and_reflector=dir('*.a15');
if length(data_and_reflector) == 0
    data_and_reflector=dir('*.a14');
elseif length(data_and_reflector) == 0
    data_and_reflector=dir('*.a11');
elseif length(data_and_reflector) == 0
    data_and_reflector=dir('*.l17'); 
elseif length(data_and_reflector) == 0
    data_and_reflector=dir('*.l21');     
end

temp = data_and_reflector.name;
YR = num2str(temp(1:4));
MONTH = num2str(temp(5:6));
DATE = num2str(temp(7:8));
abb_name = num2str(temp(end-2:end));
if abb_name == 'a15'
    target = 'apollo15';
elseif abb_name == 'a14'
    target = 'apollo14';
elseif abb_name == 'a11'
    target = 'apollo11';
elseif abb_name == 'l17'
    target = 'luna17';
elseif abb_name == 'l21'
    target = 'luna21';
end
fullURL = ['ftp://edc.dgfi.tum.de/pub/slr/cpf_predicts','/',YR,'/',target,'/','*',YR(3:4), MONTH,DATE,'*'];
filename = 'cpf';
filename = urlwrite(fullURL,filename);
opts = delimitedTextImportOptions("NumVariables", 22);
opts.DataLines = [4, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["H1", "CPF", "VarName3", "OPA", "VarName5", "VarName6", "VarName7", "VarName8"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
cpf = readtable("cpf", opts);
cpf = table2array(cpf);
clear opts
end

