function res = cinput(s,value)
% 
% res = cinput(s,value)

text = [s  ' (default '  num2str(value)  ')= ']; 
res = input(text); 
if isempty(res)
    res = value; 
end

